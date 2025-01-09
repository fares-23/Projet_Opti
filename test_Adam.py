import numpy as np
import matplotlib.pyplot as plt
import random


# IMPORTANT: Les unités sont toutes en SI, sauf indication contraire explicite.

def Simulation(Bat_cap, Train_Seuil):
    """
    Simulation d'un cas particulier, pour une capacité de batterie Bat_cap connue et un seuil de demande d'énergie à la batterie Train_Seuil donné également.
    """

    # Constantes et paramètres
    M = 70000  # Masse du train (kg)
    A_0 = 780  # Coefficient de résistance constant (N)
    A_1 = 6.4 * 10**-3  # Coefficient de résistance proportionnel à la masse (N/kg)
    B_1 = 0.14 * 10**-3 / 3.6  # Coefficient de résistance proportionnel à la vitesse (N/(m/s))
    C_0 = 0.3634 / (3.6**2)  # Coefficient de résistance proportionnel au carré de la vitesse (N/(m/s)^2)

    V_SST = 790  # Tension de la sous-station (V)
    R_SST = 33 * 10**-3  # Résistance interne de la sous-station (Ω)
    rho_LAC = 95 * 10**-6  # Résistance linéique de la LAC (Ω/m)
    rho_rail = 10 * 10**-6  # Résistance linéique des rails (Ω/m)
    Bat_rend = 0.9  # Rendement de la batterie

    # Chargement des données
    marche = np.loadtxt("marche.txt")
    t = marche[:, 0]  # Temps (s)
    x = marche[:, 1]  # Position (m)

    # Calculs de vitesse et d'accélération
    v = np.gradient(x, t[1] - t[0])
    a = np.gradient(v, t[1] - t[0])

    # Calcul des forces et puissances
    F_resistive = A_0 + A_1 * M + B_1 * M * v + C_0 * (v**2)
    F_motrice = M * a + F_resistive
    P_mecanique = F_motrice * v

    P_train = np.zeros(len(P_mecanique))
    for i in range(len(P_mecanique)):
        if P_mecanique[i] >= 0:
            P_train[i] = (P_mecanique[i] / 0.8) + 35000  # Consommation d'énergie
        else:
            P_train[i] = (P_mecanique[i] * 0.8) + 35000  # Production d'énergie

    # Résistances équivalentes
    R_LAC1 = rho_LAC * x
    R_LAC2 = rho_LAC * (x[-1] - x)
    R_rail1 = rho_rail * x
    R_rail2 = rho_rail * (x[-1] - x)

    R_eq1 = R_SST + R_LAC1 + R_rail1
    R_eq2 = R_SST + R_LAC2 + R_rail2
    R_eq = (R_eq1 * R_eq2) / (R_eq1 + R_eq2)

    # Initialisation des vecteurs
    Capabatterie = Bat_cap * 1e3 * 3600  # Conversion en Wh
    E_bat = np.zeros(len(t))
    E_bat[0] = Capabatterie

    P_bat = np.zeros(len(t))
    P_LAC = np.zeros(len(t))
    V_train = np.zeros(len(t))

    P_seuil = Train_Seuil * 10**3  # Conversion en W

    # Simulation
    for i in range(1, len(t)):
        if P_train[i] > P_seuil:  # Cas où la puissance du train dépasse le seuil
            P_LAC[i] = P_seuil
            P_bat[i] = (P_train[i] - P_LAC[i]) / Bat_rend
            E_bat[i] = E_bat[i - 1] - (P_bat[i] + P_bat[i - 1]) * 0.5

            if E_bat[i] <= 0:  # Batterie vide
                E_bat[i] = 0
                P_LAC[i] = P_train[i]
                P_bat[i] = 0

        elif P_train[i] <= P_seuil and P_train[i] > 0:  # Cas où la puissance est en dessous du seuil
            P_LAC[i] = P_train[i]
            P_bat[i] = 0
            E_bat[i] = E_bat[i - 1]

        elif P_train[i] <= 0:  # Cas de freinage
            P_LAC[i] = 0
            if E_bat[i - 1] >= Capabatterie:
                P_bat[i] = 0
                E_bat[i] = Capabatterie
            else:
                P_bat[i] = P_train[i] * Bat_rend
                E_bat[i] = E_bat[i - 1] - (P_bat[i] + P_bat[i - 1]) * 0.5

                if E_bat[i] > Capabatterie:
                    E_bat[i] = Capabatterie

    # Calcul de la tension aux bornes du train
    for i in range(len(t)):
        if (V_SST**2) - (4 * R_eq[i] * P_LAC[i]) > 0:
            V_train[i] = 0.5 * (V_SST + np.sqrt((V_SST**2) - (4 * R_eq[i] * P_LAC[i])))
        else:
            V_train[i] = 0

    # Retour de la chute de tension maximale
    return max(V_SST - V_train)


N = 1000  # Nombre de points
Capabatterie = np.random.uniform(0, 14, N)  # Capacité entre 0.5 et 14 kWh
Pseuil = np.random.uniform(0, 1000, N)  # Puissance seuil entre 0 et 1000 kW
deltaV =  []
for i in range(0,len(Capabatterie)):
    
    deltaV.append(Simulation(Capabatterie[i],Pseuil[i]))
        
deltaV = np.array(deltaV)




def non_dominated_sort(capacite, puissance):
    pareto_points = []
    for i in range(len(capacite)):
        is_dominated = False
        for j in range(len(capacite)):
            if (capacite[j] < capacite[i] and 
                puissance[j] < puissance[i] and 
                (capacite[j] < capacite[i] or puissance[j] > puissance[i])):
                is_dominated = True
                break
        if not is_dominated:
            pareto_points.append(i)
    return np.array(pareto_points)

pareto_indices = non_dominated_sort(Capabatterie, deltaV)

#4. Visualisation des résultats

plt.figure()
# Espace des solutions

plt.scatter(Capabatterie, Pseuil/1000)
plt.scatter(Capabatterie[pareto_indices], Pseuil[pareto_indices]/1000, color='orange')
plt.xlabel("Capacité (kWh)")
plt.ylabel("Puissance seuil (MW)")
plt.title("Espace des solutions")
plt.grid(visible=True)



# Espace des solutions
plt.figure()
plt.scatter(Capabatterie, deltaV)
plt.scatter(Capabatterie[pareto_indices], deltaV[pareto_indices], color='orange')
plt.xlabel("Capacité [kWh]")
plt.ylabel("Différence de tension [V]")
plt.title("Espace des objectifs")
plt.grid(visible=True)

plt.show()

print(Simulation(100000, 350000))