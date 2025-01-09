import numpy as np
import matplotlib.pyplot as plt
import random

# Charger les données depuis le fichier marche.txt
data = np.loadtxt("marche.txt")
temps = data[:, 0]  # Première colonne : temps
position = data[:, 1]  # Deuxième colonne : position

# Calculer la vitesse et l'accélération
vitesse = np.gradient(position, temps[1] - temps[0])
acceleration = np.gradient(vitesse, temps[1] - temps[0])

# Constantes
Vsst = 790
Rsst = 33e-3
rho_lac = 95e-6
rho_rail = 10e-6
M = 70000
A0 = 780
A1 = 6.4e-3
B0 = 0
B1 = 0.14e-3 / 3.6
C0 = 0.3634 / (3.6**2)
C1 = 0
g = 9.81
rendement = 0.8
conso = 35000

#Valeurs variables
#----------------Partie Mécanique----------------
R_lac1=rho_lac*position
R_rail1=rho_rail*position
R_lac2=rho_lac*(position[-1]-position)
R_rail2=rho_rail*(position[-1]-position)
F_resistive=(A0+A1*M)+(B0+B1*M)*vitesse+(C0/(10**3)+C1*M)*(vitesse**2)
F_motrice=(M*acceleration)+F_resistive
P_meca=F_motrice*vitesse
#----------------Partie Electrique----------------
Req1=(Rsst+R_lac1+R_rail1)
Req2=(Rsst+R_lac2+R_rail2)
Req=(Req1*Req2)/(Req1+Req2)
I1=Vsst/Req1
I2=Vsst/Req2
P_lac_max=min((Vsst**2)/(4*Req))
P_lac=(Vsst**2)/(4*Req)
P_train = np.copy(P_meca)

for i in range(len(P_meca)):
    if P_meca[i]>0:
        P_train[i]=(P_meca[i]/rendement) + conso
    else:
        P_train[i]=(P_meca[i]*rendement) + conso

P_train_copy=P_train.copy()
V_train=0.5*(Vsst+np.sqrt(Vsst**2-4*P_train*Req))
I_train=(Vsst-V_train)/Req

#-----------------Ajout d'une batterie-----------------
def calculer_P_train_avec_batterie(capacite_batterie, P_train, temps, vitesse, seuil):
    rendement_batt = 0.9 # Rendement de la batterie
    capacite = capacite_batterie  # Capacité maximale de la batterie (kWh)
    charge_batterie = np.zeros(len(P_train))
    P_rheo = np.zeros(len(P_train))
    P_train_copy = np.copy(P_train)

    for i in range(1, len(P_train)):
        dt = temps[i] - temps[i - 1]

        if P_train_copy[i] < 0 and charge_batterie[i - 1] < capacite:  # Freinage
            energie_stockee = abs(P_train_copy[i]) * dt / 3600 * rendement_batt
            charge_batterie[i] = min(charge_batterie[i - 1] + energie_stockee, capacite)
            if charge_batterie[i] == capacite:
                P_rheo[i] = P_train_copy[i]  # Dissiper l'excédent d'énergie

        elif P_train_copy[i] >= seuil and charge_batterie[i - 1] > 0:  # Demande énergétique élevée
            energie_dispo = charge_batterie[i - 1] * 3600 / dt * rendement_batt
            if P_train_copy[i] > energie_dispo:
                P_train_copy[i] -= energie_dispo
                charge_batterie[i] = 0
            else:
                charge_batterie[i] = charge_batterie[i - 1] - P_train_copy[i] * dt / 3600 / rendement_batt
                P_train_copy[i] = 0

        elif vitesse[i] == 0 and charge_batterie[i-1] < capacite:  # Si le train est à l'arrêt et que la batterie n'est pas pleine : utiliser la LAC pour recharger
            energie_chargee = P_lac[i] * dt / 3600 * rendement_batt
            charge_batterie[i] = min(charge_batterie[i - 1] + energie_chargee, capacite)
            if charge_batterie[i] > capacite:  # Assurer que la batterie ne dépasse pas sa capacité maximale
                charge_batterie[i] = capacite
                P_rheo[i] = P_rheo[i-1]+P_train_copy[i] + (charge_batterie[i] - capacite) * 3600 / dt

        else:  # Autres cas (charge stable)
            charge_batterie[i] = charge_batterie[i - 1]

    V_train = 0.5 * (Vsst + np.sqrt(Vsst**2 - 4 * P_train_copy * Req))
    return V_train, P_train_copy, P_rheo, charge_batterie

#-----------------Monte Carlo et Front de Pareto-----------------
capacites = np.linspace(0, 14000, 5000)  # Capacité batterie en kWh
seuils = np.linspace(0, 1e6, 5000)  # Seuil de puissance (W)

# Tirer x couples aléatoires pour l'espace des solutions
couples_aleatoires = [(random.choice(capacites), random.choice(seuils)) for _ in range(3000)]
solutions = []

for capacite, seuil in couples_aleatoires:
    V_train, _, _, _ = calculer_P_train_avec_batterie(capacite, P_train, temps, vitesse, seuil)
    chute_tension_max = Vsst - min(V_train)  # Chute de tension maximale
    solutions.append((capacite, seuil, chute_tension_max))

# Extraction du front de Pareto
pareto_solutions = []
for sol in solutions:
    dominated = False
    for other in solutions:
        if (other[0] <= sol[0] and other[2] <= sol[2] and other != sol):    # Si un autre point est meilleur sur les deux critères
            dominated = True
            break
    if not dominated:
        pareto_solutions.append(sol)

# Graphique : Front de Pareto
pareto_capacites = [sol[0] for sol in pareto_solutions]
pareto_chutes = [sol[2] for sol in pareto_solutions]

def choisir_meilleur_point(pareto_front, poids_capacite, poids_chute):
    scores = []
    for point in pareto_front:
        capacite, _ , chute = point
        score = poids_capacite * capacite + poids_chute * chute
        scores.append(score)
    
    meilleur_indice = np.argmin(scores)
    return pareto_front[meilleur_indice]

best_capacite,best_seuil,best_chute = choisir_meilleur_point(pareto_solutions, 1, 17)
print(f"Meilleure solution : Capacité = {best_capacite} Wh, Chute de tension maximale = {best_chute} V")
print("Marge de sécurité : ", Vsst - best_chute - 500, "V")

plt.figure(figsize=(10, 6))
plt.scatter([sol[0] for sol in solutions], [sol[2] for sol in solutions], alpha=0.3, label="Solutions Monte Carlo")
plt.scatter(pareto_capacites, pareto_chutes, color="red", label="Front de Pareto")
plt.scatter(best_capacite, best_chute, color="green", label="Meilleure solution")
plt.xlabel("Capacité de la batterie (kWh)")
plt.ylabel("Chute de tension maximale (V)")
plt.title("Optimisation bi-objectifs : Capacité vs Chute de tension maximale")
plt.legend()
plt.grid()


# Calculer les graphes pour la meilleure capacité et le meilleur seuil trouvés
V_train, P_train_opt, P_rheo, charge_batterie = calculer_P_train_avec_batterie(best_capacite, P_train, temps, vitesse, best_seuil)

# Tracer les graphes
plt.figure(figsize=(14, 10))

# Graphique 1 : Puissance
plt.subplot(3, 1, 1)
plt.plot(temps, P_train_opt, label="Puissance effective (W)")
plt.plot(temps, P_rheo, label="Puissance dissipée (rhéostatique, W)")
plt.xlabel("Temps (s)")
plt.ylabel("Puissance (W)")
plt.title("Puissance en fonction du temps")
plt.legend()
plt.grid()

# Graphique 2 : Tension
plt.subplot(3, 1, 2)
plt.plot(temps, V_train, label="Tension (V)")
plt.axhline(790, color="red", linestyle="--", label="Tension nominale (790V)")
plt.xlabel("Temps (s)")
plt.ylabel("Tension (V)")
plt.title("Tension en fonction du temps")
plt.legend()
plt.grid()

# Graphique 3 : Charge de la batterie
plt.subplot(3, 1, 3)
plt.plot(temps, charge_batterie, label="Charge de la batterie (kWh)")
plt.xlabel("Temps (s)")
plt.ylabel("Charge (Wh)")
plt.title("Charge de la batterie en fonction du temps")
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()
