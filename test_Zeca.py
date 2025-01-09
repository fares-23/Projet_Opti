import numpy as np
import matplotlib.pyplot as plt
import random

def Simulation(Bat_cap,seuil):
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
    V_SST = 790 # Tension délivrée par la sous-station (tension nominale).
    R_SST = 33*10**(-3) # Résistance interne de la sous-station.
    rho_LAC = 95*10**(-6) # Résistance linéique de la LAC.
    rho_rail = 10*10**(-6) # Résistance linéique des rails.
    x_Tot = 5000 # Distance totale que le train doit parcourir durant son trajet.
    SysBor = 35000 # Consommation du système de bord.
    rend = 0.8 # Rendement du moteur+convertisseurs

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
    F_resistive = (A_0+A_1*M)+(B_0+B_1*M)*v+(C_0+C_1*M)*v**2

    # F_motrice.
    F_motrice = M*a+M*g*np.sin(alpha)+F_resistive

    # P_mecanique.
    P_mecanique = F_motrice*v

    # Puissance totale consommée par le train: P_train.
    P_train = np.zeros(len(P_mecanique))
    for k in range(len(P_mecanique)):
        if P_mecanique[k] > 0:
            P_train[k] = P_mecanique[k]/rend+SysBor # On consomme l'énergie.
        else:
            P_train[k] = P_mecanique[k]*rend+SysBor # On produit l'énergie.

    """
        BATTERIE
        ========
    """
    # Bat_cap = 100000 # Capacité de la batterie (Wh).
    Bat_E = np.zeros(len(P_train)) # Énergie contenue dans la batterie (à un instant t).
    Bat_Charge = np.zeros(len(P_train)) # Charge de la batterie, au cours du temps.
    Bat_Charge[0] = 0 # Charge de départ de la batterie.
    if Bat_Charge[0] > Bat_cap*3600 and Bat_Charge < 0: # On vérifie qu'on ne dépasse pas le capacité de la batterie.
        raise ValueError("On dépasse les capacités de la batterie!")
    Train_Seuil = seuil # Seuil de puissance (négatif, en joule) auquel le train demande de l'énergie.
    Rheo_P = np.zeros(len(P_train)) # Puissance dissipée par le rhéostat.
    Rheo_P[0] = 0 # Le rhéostat ne dissipe rien au départ.

    for k in range(1, len(t)):
        if P_train[k] < 0: # Le train freine.
            Bat_E[k] = Bat_E[k-1]-P_train[k]
            P_train[k] = 0
            if Bat_E[k] > 3600*Bat_cap: # Dépasse-t-on les limites de la batterie?
                Diff = Bat_E[k] - Bat_cap*3600
                Bat_E[k] = Bat_cap*3600
                Rheo_P[k] = Diff # On dissipe dans le rhéostat.
        elif P_train[k] >= Train_Seuil: # Si le train demande de l'énergie.
            Bat_E[k] = Bat_E[k-1] - (P_train[k]-Train_Seuil)
            P_train[k] = Train_Seuil
            if Bat_E[k] < 0: # Dépasse-t-on les limites de la batterie?
                Diff = -Bat_E[k]
                Bat_E[k] = 0
                P_train[k] += Diff
        else:
            Bat_E[k] = Bat_E[k-1] # La batterie conserve son énergie.

    """
        RÉSEAU FERROVIAIRE
        ==================
    """

    # Résistances du circuit.
    R_LAC1 = rho_LAC*x
    R_LAC2 = (x_Tot-x)*rho_LAC
    R_rail1 = x*rho_rail
    R_rail2 = (x_Tot-x)*rho_rail
    R_eq = ((R_SST+R_LAC1+R_rail1)*(R_SST+R_LAC2+R_rail2))/(2*R_SST+R_LAC1+R_LAC2+R_rail1+R_rail2)

    P_LAC = np.zeros(len(P_train))
    for k in range(len(P_train)):
        if P_train[k] > (V_SST**2)/(4*R_eq[k]):
            P_LAC[k] = (V_SST**2)/(4*R_eq[k])-SysBor
        else:
            P_LAC[k] = P_train[k]

    # Tensions aux bornes du train.
    V_train = 0.5*(V_SST+np.sqrt(V_SST**2-4*R_eq*P_LAC))

    # Courant tranversant le train.
    I_train = P_train/V_train

    # Courants dans chaque branche.
    I_1 = (V_SST-V_train)/(R_SST+R_LAC1+R_rail1)
    I_2 = (V_SST-V_train)/(R_SST+R_LAC2+R_rail2)

    # Puissances fournies par les sous-stations
    P_SST1 = I_1*V_SST
    P_SST2 = I_2*V_SST

    # Puissance perdue par sous-station.
    P_SST1_loss = (R_SST+R_LAC1+R_rail1)*I_1**2
    P_SST2_loss = (R_SST+R_LAC2+R_rail2)*I_2**2

    # Valeur maximale entre V_SST et la tension aux bornes du train (indicateur de qualité).
    Ind_qual = V_SST - np.min(V_train)
    IndCrit = V_SST - 500 # Valeur critique du Ind_qual, à ne pas dépasser.

    return Ind_qual # Indicateur de qualité = Chute de tension maximale.

print(Simulation(10000,35000))

def Paires_P_seuil_CapaBat():
    capacites = np.linspace(0, 14000, 5000)  # Capacité batterie en kWh
    seuils = np.linspace(0, 1e6, 5000)  # Seuil de puissance (W)

    # Tirer x couples aléatoires pour l'espace des solutions
    couples_aleatoires = [(random.choice(capacites), random.choice(seuils)) for _ in range(1000)]
    solutions = []

    for capacite, seuil in couples_aleatoires:
        chute_tension_max= Simulation(capacite,seuil)
        solutions.append((capacite, seuil, chute_tension_max))

    capacites, seuils, chute_tension_max_values = zip(*solutions)
    plt.figure(figsize=(6, 4))
    plt.scatter(capacites, chute_tension_max_values)  # Scatter plot de l'espace des solutions
    plt.title("Espace des solutions", fontsize=14)
    plt.xlabel("Capacité de la batterie (kWh)", fontsize=12)
    plt.ylabel("Chute tension (V)", fontsize=12)
    plt.grid(True)
    plt.show()
    
    return solutions

solutions = Paires_P_seuil_CapaBat()

pareto_solutions = []
for sol in solutions:
    dominé = False
    for options in solutions:
        if (options[0] <= sol[0] and options[0] <= sol[2] and options != sol):
            dominated = True 
            break
    if not dominé:
        pareto_solutions.append(sol)
    

# Graphique : Front de Pareto
pareto_capacites = [sol[0] for sol in pareto_solutions]
pareto_chutes = [sol[2] for sol in pareto_solutions]
    

def MonteCarlo(CapaBat_max, CapaBat_min, CapaBat_step):
    """
        Teste de Monte-Carlo pour une liste de valeurs de capacités de batterie ainsi que de seuils de demandes d'énergie lorsque le train roule (ma poule).
    """
    CapaBat = [] # Création des valeurs de la capacité de la batterie.
    for i in range(CapaBat_max, CapaBat_min, CapaBat_step):
        CapaBat.append(i)

    V_SST = 790 # Tension délivrée par la sous-station (tension nominale).
    IndCrit = V_SST - 500 # Valeur critique du Ind_qual, à ne pas dépasser.
    Sim = [] # Liste des paramètres de chaque simulation.
    Qual = [] # Liste des résultats de chaque simulation.
    for capa in CapaBat:
        para = f"{capa}" # Paramètres de la simulation.
        Sim.append(para)
        Qual.append(Simulation(capa))
        print(Simulation(capa))

    # # Affichage des résultats.
    # plt.figure("Résultats des Simulations de Monte-Carlo")
    # plt.axhline(y=IndCrit, color="red", linestyle="-")
    # plt.bar(Sim, Qual, color="grey", edgecolor="black")
    # plt.title("Résultats des Simulations de Monte-Carlo")
    # plt.xlabel("Capacité de la Batterie [Wh]")
    # plt.ylabel("Différence de Tension [V]")
    # plt.grid()
    # plt.legend()
    # plt.show()


"""
    MONTE-CARLO
    ===========
"""

#MonteCarlo(0, 2000, 10)

"""
        Simulation d'un cas particulier, pour une capacité de batterie Bat_cap connue et un P_seuil de demande d'énergie à la batterie Train_Seuil donnée également.
    """

"""
        VARIABLES
        =========
    """

    # M = 70000 # Masse du train.
    # A_0 = 780
    # A_1 = 6.4*10**(-3)
    # B_0 = 0
    # B_1 = 0.14*3.6*10**(-3)
    # C_0 = 0.3634*3.6**2
    # C_1 = 0
    # alpha = 0 # Pente que le train doit gravir (compté positivement).
    # g = 9.81 # Accélération de la pesanteur.
    # V_SST = 790 # Tension délivrée par la sous-station (tension nominale).
    # R_SST = 33*10**(-3) # Résistance interne de la sous-station.
    # rho_LAC = 131*10**(-6) # Résistance linéique de la LAC.
    # rho_rail = 18*10**(-6) # Résistance linéique des rails.
    # x_Tot = 5000 # Distance totale que le train doit parcourir durant son trajet.
    # SysBor = 35000 # Consommation du système de bord.
    # rend = 0.8 # Rendement du moteur+convertisseurs.



    # """
    #     GRANDEURS MESURÉES
    #     ==================
    # """

    # marche = np.loadtxt("marche.txt")
    # t = marche[:, 0] # Temps.
    # x = marche[:, 1] # Distance parcourue par le train.



    # """
    #     TRAIN
    #     =====
    # """

    # # Vitesse v.
    # v = np.gradient(x, t[1]-t[0])


    # # Accélération a.
    # a = np.gradient(v, t[1]-t[0])

    # # F_resistive.
    # F_resistive = (A_0+A_1*M)+(B_0+B_1*M)*v+(C_0+C_1*M)*v**2

    # # F_motrice.
    # F_motrice = M*a+M*g*np.sin(alpha)+F_resistive

    # # P_mecanique.
    # P_mecanique = F_motrice*v

    # # Puissance totale consommée par le train: P_train.
    # P_train = np.zeros(len(P_mecanique))
    # for k in range(len(P_mecanique)):
    #     if P_mecanique[k] > 0:
    #         P_train[k] = P_mecanique[k]/rend+SysBor # On consomme l'énergie.
    #     else:
    #         P_train[k] = P_mecanique[k]*rend+SysBor # On produit l'énergie.

    # """
    #     RÉSEAU FERROVIAIRE
    #     ==================
    # """

    # # Résistances du circuit.
    # R_LAC1 = rho_LAC*x
    # R_LAC2 = (x_Tot-x)*rho_LAC
    # R_rail1 = x*rho_rail
    # R_rail2 = (x_Tot-x)*rho_rail
    # R_eq = ((R_SST+R_LAC1+R_rail1)*(R_SST+R_LAC2+R_rail2))/(2*R_SST+R_LAC1+R_LAC2+R_rail1+R_rail2)

    # P_LAC = np.zeros(len(P_train))
    # for k in range(len(P_train)):
    #     if P_train[k] > (V_SST**2)/(4*R_eq[k]):
    #         P_LAC[k] = (V_SST**2)/(4*R_eq[k])-SysBor
    #     else:
    #         P_LAC[k] = P_train[k]

    # # Tensions aux bornes du train.
    # V_train = 0.5*(V_SST+np.sqrt(V_SST**2-4*R_eq*P_LAC))

    # # Courant tranversant le train.
    # I_train = P_train/V_train

    # # Courants dans chaque branche.
    # I_1 = (V_SST-V_train)/(R_SST+R_LAC1+R_rail1)
    # I_2 = (V_SST-V_train)/(R_SST+R_LAC2+R_rail2)

    # # Puissances fournies par les sous-stations
    # P_SST1 = I_1*V_SST
    # P_SST2 = I_2*V_SST

    # # Puissance perdue par sous-station.
    # P_SST1_loss = (R_SST+R_LAC1+R_rail1)*I_1**2
    # P_SST2_loss = (R_SST+R_LAC2+R_rail2)*I_2**2

    # # Valeur maximale entre V_SST et la tension aux bornes du train (indicateur de qualité).
    # Ind_qual = V_SST - np.min(V_train)
    # IndCrit = V_SST - 500 # Valeur critique du Ind_qual, à ne pas dépasser.


    # """
    #     BATTERIE
    #     ========
    # """
