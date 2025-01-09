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
    V_SST = 790 # Tension délivrée par la sous-station (tension nominale).
    R_SST = 33*10**(-3) # Résistance interne de la sous-station.
    rho_LAC = 95*10**(-6) # Résistance linéique de la LAC.
    rho_rail = 10*10**(-6) # Résistance linéique des rails.
    x_Tot = 5000 # Distance totale que le train doit parcourir durant son trajet.
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
    Bat_rend = 0.9 # Rendement de la batterie.
    Bat_E = np.zeros(len(P_train)) # Énergie contenue dans la batterie (à un instant t).
    Bat_Charge = np.zeros(len(P_train)) # Charge de la batterie, au cours du temps.
    Bat_Charge[0] = 0 # Charge de départ de la batterie.
    if Bat_Charge[0] > Bat_cap*3600 and Bat_Charge < 0: # On vérifie qu'on ne dépasse pas le capacité de la batterie.
        raise ValueError("On dépasse les capacités de la batterie!")
    # Train_Seuil = 350000 # Seuil de puissance (négatif, en joule) auquel le train demande de l'énergie.
    Rheo_P = np.zeros(len(P_train)) # Puissance dissipée par le rhéostat.
    Rheo_P[0] = 0 # Le rhéostat ne dissipe rien au départ.

    for k in range(1, len(t)):
        if P_train[k] < 0: # Le train freine.
            Bat_E[k] = Bat_E[k-1]-P_train[k]*Bat_rend
            P_train[k] = 0
            if Bat_E[k] > 3600*Bat_cap: # Dépasse-t-on les limites de la batterie?
                Diff = Bat_E[k] - Bat_cap*3600
                Bat_E[k] = Bat_cap*3600
                Rheo_P[k] = Diff # On dissipe dans le rhéostat.
        elif P_train[k] >= Train_Seuil: # Si le train demande de l'énergie.
            Bat_E[k] = Bat_E[k-1] - (P_train[k]-Train_Seuil)
            P_train[k] = P_train[k] - (P_train[k]-Train_Seuil) * Bat_rend
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
            P_LAC[k] = (V_SST**2)/(4*R_eq[k])
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



    """
        AFFICHAGE
        =========
    """

    # plt.figure("Position, vitesse, accélération du train en fonction du temps")
    # # Affichage position :
    # plt.subplot(3, 1, 1)
    # plt.plot(t, x/1000, "-k", label="Position du train") # Position normalisée en km.
    # plt.title("Position, vitesse, accélération du train en fonction du temps")
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Longueur [km]")
    # plt.grid()
    # plt.legend()
    #
    # # Affichage vitesse :
    # plt.subplot(3, 1, 2)
    # plt.plot(t, v/1000, "-k", label="Vitesse du train") # Vitesse normalisée en km/s.
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Vitesse [km/s]")
    # plt.grid()
    # plt.legend()
    # # Affichage accélération :
    # plt.subplot(3, 1, 3)
    # plt.plot(t, a/g, "-k", label="Accélération du train") # Accélération normalisée en g.
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Accélération [g]")
    # plt.grid()
    # plt.legend()
    #
    #
    # plt.figure("Puissance, tension et courant dans le train")
    # # Affichage de la puissance :
    # plt.subplot(3, 1, 1)
    # plt.plot(t, P_train/1000000, "-k", label="Puissance consommée par le train") # P_train normalisée en MW.
    # plt.title("Puissance, tension et courant dans le train")
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Puissance [MW]")
    # plt.legend()
    # plt.grid()
    #
    # # Affichage de la tension:
    # plt.subplot(3, 1, 2)
    # plt.plot(t, V_train, "-k", label="Tensions aux bornes de la locomotive")
    # plt.legend()
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Tension [V]")
    # plt.grid()
    #
    # # Affichage du courant:
    # plt.subplot(3, 1, 3)
    # plt.plot(t, I_train, "-k", label="Courant traversant le train")
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Courant [A]")
    # plt.legend()
    # plt.grid()
    #
    #
    # plt.figure("Courants dans les branches du circuit")
    # # Affichage de I_1, I_2 et I_train:
    # plt.plot(t, I_1, "-b", label="I_1")
    # plt.plot(t, I_2, "-g", label="I_2")
    # plt.plot(t, I_train, "-k", label="I_train")
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Courant [A]")
    # plt.title("Courants dans les branches du circuit")
    # plt.grid()
    # plt.legend()
    #
    #
    # plt.figure("Puissances mises en jeu")
    # # Affichage Puissance des stations et du train.
    # plt.subplot(3, 1, 1)
    # plt.plot(t, P_SST1/1000000, "-b", label="Sous-station 1") # normalisé en MW
    # plt.plot(t, P_SST2/1000000, "-g", label="Sous-station 2") # normalisé en MW
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Puissance [MW]")
    # plt.title("Puissances mises en jeu")
    # plt.legend()
    # plt.grid()
    #
    # # Affichage de la puissance des stations avec pertes.
    # plt.subplot(3, 1, 2)
    # plt.plot(t, P_SST1_loss/1000000, "-m", label="Perte 1") # normalisé en MW
    # plt.plot(t, P_SST2_loss/1000000, "-c", label="Perte 2") # normalisé en MW
    # plt.grid()
    # plt.legend()
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Puissance [MW]")
    #
    # # Affichage de la puissance des stations avec pertes.
    # plt.subplot(3, 1, 3)
    # plt.plot(t, (P_SST1+P_SST2-P_SST1_loss-P_SST2_loss)/1000000, "-r", label="Puissance des stations, avec pertes") # normalisé en MW
    # plt.plot(t, P_train/1000000, "-k", label="Puissance consommée par le train") # P_train normalisé en MW
    # plt.grid()
    # plt.legend()
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Puissance [MW]")
    #
    #
    # plt.figure("Gestion de la Batterie")
    # # Affichage de l'énergie de la batterie, de la consommation du rhéostat et de l'énergie du train ainsi que de sa tension.
    # plt.subplot(3, 1, 1)
    # plt.plot(t, Bat_E/1000000, "-b", label="Énergie dans la batterie")
    # plt.plot(t, Rheo_P/1000000, "-g", label="Énergie perdue dans le rhéostat")
    # plt.title("Gestion de la Batterie")
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Énergie [MJ]")
    # plt.grid()
    # plt.legend()
    #
    # plt.subplot(3, 1, 2)
    # plt.plot(t, P_train/1000000, "-k", label="Puissance consommée par le train")
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Puissance [MW]")
    # plt.grid()
    # plt.legend()
    #
    # plt.subplot(3, 1, 3)
    # plt.plot(t, V_train, "-k", label="Tension aux bornes du train")
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Tension [V]")
    # plt.grid()
    # plt.legend()
    #
    #
    # plt.figure("Indication de la Qualité du Système")
    # # Comparaison, dans un histogramme, de Ind_qual et IndCrit.
    # plt.bar(["Indicateur Critique\n(à ne pas dépasser)", "Indicateur Actuelle"], [IndCrit, Ind_qual], color="grey", edgecolor="black")
    # plt.title("Indication de la Qualité du Système")
    # plt.ylabel("Différence de Tension [V]")
    # plt.grid()
    # plt.legend()

    return Ind_qual # Indicateur de qualité = Chute de tension maximale.

print(Simulation(100000, 350000))
def non_dominant_sort(pop):
    """
        Renvoie les rangs, contenant différents individus suivant la dominance de ces derniers (1er rang = les non-dominées, 2e rang = les autres non-dominés, etc...).
    """
    fronts = []
    front = []
    for individual in pop:
        dominated = False # Par défaut, on considère l'individu non-dominé.
        for other in pop:
            if (other[0] < individual[0] and other[1] <= individual[1] and other[2] <= individual[2]) or (other[0] <= individual[0] and other[1] < individual[1] and other[2] <= individual[2]) or (other[0] <= individual[0] and other[1] <= individual[1] and other[2] < individual[2]):
                dominated = True # L'individu est dominé.
                break
        if not dominated:
            front.append(individual) # On n'y met que les individus non-dominés.
    fronts.append(front)
    while len(front) > 0:
        new_front = []
        for individual in front:
            for other in pop:
                if (other[0] < individual[0] and other[1] <= individual[1] and other[2] <= individual[2]) or (other[0] <= individual[0] and other[1] < individual[1] and other[2] <= individual[2]) or (other[0] <= individual[0] and other[1] <= individual[1] and other[2] < individual[2]):
                    new_front.append(other)
        fronts.append(new_front)
        front = new_front
    return fronts



def MonteCarlo(CapaBat_min, CapaBat_max, CapaBat_step, Seuil_min, Seuil_max, Seuil_step):
    """
        Teste de Monte-Carlo pour une liste de valeurs de capacités de batterie ainsi que de seuils de demandes d'énergie lorsque le train roule (ma poule).
    """
    population = []
    CapaBat = [] # Création des valeurs de la capacité de la batterie.
    for i in range(CapaBat_min, CapaBat_max, CapaBat_step):
        CapaBat.append(i)
    TrainSeuil = [] # Création des valeurs du seuil.
    for i in range(Seuil_min, Seuil_max, Seuil_step):
        TrainSeuil.append(i)

    for capa in CapaBat:
        for seuil in TrainSeuil:
            population.append([capa, seuil, Simulation(capa, seuil)]) # Création des individus.
    fronts = non_dominant_sort(population) # Rangement des individus.

    # Affichage des résultats.
    plt.figure("Espace des Solutions / Espace des Objectifs")
    plt.subplot(2, 1, 1)
    for fronts in fronts[:0]:
        for individual in fronts:
            plt.plot(individual[0], individual[1], "+k")
    for individual in fronts[0]:
        plt.plot(individual[0], individual[1], "+k")
    plt.xlabel("Capacités de la Batterie [Wh]")
    plt.ylabel("Seuils [J]")
    plt.title("Espace des Solutions / Espace des Objectifs")
    plt.grid()

    plt.subplot(2, 1, 2)
    for front in fronts[:0]:
        for individual in front:
            plt.plot(fronts[0], Simulation(individual[0], individual[1]), "+k")
    for individual in fronts[0]:
        plt.plot(individual[0], Simulation(individual[0], individual[1]), "+b") # Dernière génération.
    plt.xlabel("Capacités de la Batterie [Wh]")
    plt.ylabel("Chutes de Tension [V]")
    plt.grid()



def NSGA2(Bounds_Capa, Step_Capa, Bounds_Seuil, Step_Seuil, PopSize, NumGen, CrossProba=0.5, MutationProba=0.25):
    """
        Utilise l'algorithme NSGA-II pour trouver les meilleures solutions entre la capacité de la batterie et la chute de tension.

        - PopSize: Taille de la population.
        - NumGen: Nombre de générations à générer.
        - Bounds_Capa: Intervalle des valeurs de la capacité de la batterie.
        - Step_Capa: Écart entre deux valeurs adjacentes de la capacité de la batterie.
        - Bouds_Seuil: Intervalle des valeurs du seuil.
        - Step_Seuil: Étape entre chaque valeur du seuil.
        - CrossProba: Probabilité d'avoir une crossover (défaut: 0.5).
        - MutationProba: Probabilité d'avoir une mutation (défaut: 0.1).
    """
    process = [] # Sauvegarde des populations.
    population = [[random.choice(np.arange(Bounds_Capa[0], Bounds_Capa[1], Step_Capa)), random.choice(np.arange(Bounds_Seuil[0], Bounds_Seuil[1], Step_Seuil))] for _ in range(PopSize)]
    process.append(population)
    V_SST = 790 # Tension délivrée par la sous-station (tension nominale).

    # Évolution de la populace.
    for _ in range(NumGen):
        fronts = non_dominant_sort([[individual[0], individual[1], Simulation(individual[0], individual[1])] for individual in population]) # Classement des individus.
        best = [] # Liste des meilleurs individus.
        out = False
        for front in fronts:
            for individual in front:
                best.append(individual) # On rajoute le meilleur.
                if len(best) >= 0.5 * PopSize:
                    out = True # On a sélectionné 50 % des meilleurs.
                    break
            if out:
                break
        offspring = [] # Liste des enfants.
        if len(best) < 2: # Ensuring best is big enough.
            best = population
        for _ in range(PopSize): # On crée les enfants.
            parent1, parent2 = random.sample(best, 2) # Sélection des parents.
            child = [0.5 * (parent1[0] + parent2[0]), 0.5 * (parent1[1] + parent2[1])] # Croisement.
            if random.random() < MutationProba: # Mutation de la batterie?
                child[0] += random.uniform(-Step_Capa, Step_Capa)
                child[0] = max(Bounds_Capa[0], min(Bounds_Capa[1], child[0]))
            if random.random() < MutationProba: # Mutation du seuil?
                child[1] += random.uniform(-Step_Seuil, Step_Seuil)
                child[1] = max(Bounds_Seuil[0], min(Bounds_Seuil[1], child[0]))
            offspring.append(child)
        population = offspring # Nouvelle population.
        process.append(population)

    # Affichage.
    plt.figure("Espace des Solutions / Espace des Objectifs")
    plt.subplot(2, 1, 1)
    for pop in process[:-1]:
        for individual in pop:
            plt.plot(individual[0], individual[1], "+k")
    for individual in process[-1]:
        plt.plot(individual[0], individual[1], "+b") # Dernière génération.
    plt.xlabel("Capacités de la Batterie [Wh]")
    plt.ylabel("Seuils [J]")
    plt.title("Espace des Solutions / Espace des Objectifs")
    plt.grid()

    plt.subplot(2, 1, 2)
    for pop in process[:-1]:
        for individual in pop:
            plt.plot(individual[0], Simulation(individual[0], individual[1]), "+k")
    for individual in process[-1]:
        plt.plot(individual[0], Simulation(individual[0], individual[1]), "+b") # Dernière génération.
    plt.xlabel("Capacités de la Batterie [Wh]")
    plt.ylabel("Chutes de Tension [V]")
    plt.grid()




"""
    Teste de SIMULATION
    ===================
"""

Simulation(10000, 300000)



"""
    MONTE-CARLO
    ===========
"""

# MonteCarlo(1000, 10000, 1000, 0, 1000000, 100000)


"""
    NSGA-II
    =======
"""

# NSGA2([1000, 10000], 1000, [0, 100000], 100, 50, 50)

plt.show()
