import numpy as np
import matplotlib.pyplot as plt
import random


# IMPORTANT: Les unités sont toutes en SI, sauf indication contraire explicite.


def Simulation(Bat_cap, Train_Seuil):
    """
        Simulation d'un cas particulier, pour une capacité de batterie Bat_cap connue et un seuil de demande d'énergie à la batterie Train_Seuil donnée également.
    """

    """
        VARIABLES
        =========
    """

    M = 70000 # Masse du train.
    A_0 = 780
    A_1 = 6.4*10**(-3)
    B_0 = 0
    B_1 = 0.14*3.6*10**(-3)
    C_0 = 0.3634*3.6**2
    C_1 = 0

    alpha = 0 # Pente que le train doit gravir (compté positivement).
    g = 9.81 # Accélération de la pesanteur.
    
    SysBor = 35000 # Consommation du système de bord.
    rend = 0.8 # Rendement du moteur+convertisseurs.

    """
        GRANDEURS MESURÉES
        ==================
    """

    marche = np.loadtxt("marche.txt")
    t = marche[:, 0] # Temps.
    x = marche[:, 1] # Distance parcourue par le train.

    """
        TRAIN
        =====
    """

    # Vitesse v.
    v = np.gradient(x, t[1]-t[0])

    # Accélération a.
    a = np.gradient(v, t[1]-t[0])

    # F_resistive.
    F_resistive = (A_0+A_1*M)+(B_0+B_1*M)*v+(C_0+C_1*M)*(v**2)

    # F_motrice.
    F_motrice = M*a + F_resistive + M*g*np.sin(alpha)

    # P_mecanique.
    P_mecanique = F_motrice*v

    # Puissance totale consommée par le train: P_train.
    P_train = np.zeros(len(P_mecanique))
    
    for k in range(len(P_mecanique)):
        if P_mecanique[k] > 0:
            P_train[k] = P_mecanique[k]/rend+SysBor # On consomme l'énergie.
        else:
            P_train[k] = P_mecanique[k]*rend+SysBor # On produit l'énergie.

    # Puissance et tension du train SANS batterie
    Vsst = 790 # V
    Rsst = 33*10**(-3) # Ohm
    pLac = 95*10**(-6) # Ohm/m
    pRail = 10*10**(-6) # Ohm/m
    
    
    return 0
    

print(Simulation(100000, 350000)) # Test de la fonction Simulation.

N = 1000  # Nombre de points
Capabatterie = np.random.uniform(0, 14, N)  # Capacité entre 0.5 et 14 kWh
Pseuil = np.random.uniform(0, 1000, N)  # Puissance seuil entre 0 et 1000 kW
deltaV =  []
for i in range(0,len(Capabatterie)):
    
    deltaV.append(Simulation(Capabatterie[i],Pseuil[i]))
        
deltaV = np.array(deltaV)


def nom_dominated_sort(capacite, puissance):
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

pareto_indices = nom_dominated_sort(Capabatterie, deltaV)

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