import numpy as np
import matplotlib.pyplot as plt
import random


# IMPORTANT: Les unités sont toutes en SI, sauf indication contraire explicite.


def Simulation(Bat_cap,seuil):
    """
        Simulation d'un cas particulier, pour une capacité de batterie Bat_cap connue et un seuil de demande d'énergie à la batterie seuil donnée également.
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
    Bat_rend = 0.9 # Rendement de la batterie.
    Bat_E = np.zeros(len(P_train)) # Énergie contenue dans la batterie (à un instant t).
    Bat_Charge = np.zeros(len(P_train)) # Charge de la batterie, au cours du temps.
    Bat_Charge[0] = 0 # Charge de départ de la batterie.
    if Bat_Charge[0] > Bat_cap*3600 and Bat_Charge < 0: # On vérifie qu'on ne dépasse pas le capacité de la batterie.
        raise ValueError("On dépasse les capacités de la batterie!")
    Rheo_P = np.zeros(len(P_train)) # Puissance dissipée par le rhéostat.
    Rheo_P[0] = 0 # Le rhéostat ne dissipe rien au départ.

    for k in range(1, len(t)):
        if P_train[k] < 0: # Le train freine.
            Bat_E[k] = Bat_E[k-1]-P_train[k] * Bat_rend
            P_train[k] = 0
            if Bat_E[k] > 3600*Bat_cap: # Dépasse-t-on les limites de la batterie?
                Diff = Bat_E[k] - Bat_cap*3600
                Bat_E[k] = Bat_cap*3600
                Rheo_P[k] = Diff # On dissipe dans le rhéostat.
        elif P_train[k] >= seuil: # Si le train demande de l'énergie.
            Bat_E[k] = Bat_E[k-1] - (P_train[k]-seuil)
            P_train[k] = P_train[k] - seuil * Bat_rend
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



def non_dominant_sort(pop):
    """
        Renvoie les rangs, contenant différents individus suivant la dominance de ces derniers (1er rang = les non-dominées, 2e rang = les autres non-dominés, etc).
    """
    fronts = []
    front = [] # Premier rang.
    population = pop[:] # Initialisation de la population.
    while len(population) > 0:
        for individual in population:
            dominated = False # On suppose que l'individu n'est pas dominé.
            for other in population:
                if (other[0] < individual[0] and other[1] <= individual[1] and other[2] <= individual[2]) or (other[0] <= individual[0] and other[1] < individual[1] and other[2] <= individual[2]) or (other[0] <= individual[0] and other[1] <= individual[1] and other[2] < individual[2]):
                    dominated = True # L'individu est dominé.
                    break
            if not dominated:
                front.append(individual[:])
        if len(front) == 0:
            front = population[:]
        fronts.append(front[:])
        for individual in front:
            population.remove(individual) # On retire les individus déjà classés.
        front = [] # On réinitialise le premier rang.
    return fronts



def MonteCarlo():
    capacites = np.linspace(0, 14000, 1000)  # Capacité batterie en kWh
    seuils = np.linspace(0, 1e6, 1000)  # Seuil de puissance (W)

    # Tirer x couples aléatoires pour l'espace des solutions
    couples_aleatoires = [(random.choice(capacites), random.choice(seuils)) for _ in range(1000)]
    solutions = []

    for capacite, seuil in couples_aleatoires:
        chute_tension_max= Simulation(capacite,seuil)
        solutions.append((capacite, seuil, chute_tension_max))

    # capacites, seuils, chute_tension_max_values = zip(*solutions)
    # plt.figure(figsize=(6, 4))
    # plt.scatter(capacites, chute_tension_max_values)  # Scatter plot de l'espace des solutions
    # plt.title("Espace des solutions", fontsize=14)
    # plt.xlabel("Capacité de la batterie (kWh)", fontsize=12)
    # plt.ylabel("Chute tension (V)", fontsize=12)
    # plt.grid(True)
    # plt.show()

    return solutions



def non_dominant_sort(pop):
    """
        Renvoie les rangs, contenant différents individus suivant la dominance de ces derniers (1er rang = les non-dominées, 2e rang = les autres non-dominés, etc).
    """
    fronts = []
    front = [] # Premier rang.
    population = pop[:] # Initialisation de la population.
    while len(population) > 0:
        for individual in population:
            dominated = False # On suppose que l'individu n'est pas dominé.
            for other in population:
                if (other[0] <= individual[0] and other[2] <= individual[2] and other != individual):
                    dominated = True # L'individu est dominé.
                    break
            if not dominated:
                front.append(individual[:])
        if len(front) == 0:
            front = population[:]
        fronts.append(front[:])
        for individual in front:
            population.remove(individual) # On retire les individus déjà classés.
        front = [] # On réinitialise le premier rang.
    return fronts


def select_half_best(fronts, popSize):
    """
        Retourne une liste des 50% meilleurs individus de rangs.
    """
    best = []
    for front in fronts:
        for individual in front:
            best.append(individual)
            if len(best) >= popSize // 2:
                return best
    return best


def NSGA2(CapaLim, CapaStep, SeuilLim, SeuilStep, PopSize, N, mutant=0.25):
    """
        Algorithme NSGA-II (simpififé).

        - CapaLim: [Capacité Minimale, Capacité Maximale]
        - CapaStep: Pas pour la capacité.
        - SeuilLim: [Seuil_Mini, Seuil_Maxi]
        - SeuilStep: Pas du seuil.
        - PopSize: Taille de la population.
        - N: Nombre d'itérations.
        - mutant: Probabilité d'une mutation.
    """
    process = [] # Ensemble des populations.
    population = [] # Création de la population.
    for _ in range(PopSize):
        capa = random.choice(np.arange(CapaLim[0], CapaLim[1], CapaStep))
        seuil = random.choice(np.arange(SeuilLim[0], SeuilLim[1], SeuilStep))
        population.append([capa, seuil, Simulation(capa, seuil)])

    for _ in range(N):
        fronts = non_dominant_sort(population)
        best = select_half_best(fronts, PopSize)
        new_pop = [] # Nouvelle génération.
        for _ in range(PopSize):
            parent1, parent2 = random.sample(best, 2)
            capa = (parent1[0] + parent2[0]) * 0.5
            seuil = (parent1[1] + parent2[1]) * 0.5
            if random.random() <= mutant: # Mutation?
                capa += random.uniform(-CapaStep, CapaStep)
                seuil += random.uniform(-SeuilStep, SeuilStep)
            if capa < CapaLim[0]:
                capa = CapaLim[0]
            if capa > CapaLim[1]:
                capa = CapaLim[1]
            if seuil < SeuilLim[0]:
                seuil = SeuilLim[0]
            if seuil > SeuilLim[1]:
                seuil = SeuilLim[1]
            child = [capa, seuil, Simulation(capa, seuil)] # Création d'un nouvel individu.
            new_pop.append(child[:])
        process.append(population[:])
        population = new_pop[:] # MÀJ.

    # Affichage:
    plt.figure("NSGA-II")
    plt.subplot(2, 1, 1)
    for pop in process[:-1]:
        for individual in pop:
            plt.plot(individual[0], individual[1], "+k")
    for individual in process[-1]:
        plt.plot(individual[0], individual[1], "+b")
    plt.title("Espace des Solutions / Espace des Objectifs")
    plt.xlabel("Capacités de la Batterie [Wh]")
    plt.ylabel("Seuils [W]")
    plt.subplot(2, 1, 2)
    for pop in process[:-1]:
        for individual in pop:
            plt.plot(individual[0], individual[2], "+k")
    for individual in process[-1]:
        plt.plot(individual[0], individual[2], "+b")
    plt.xlabel("Capacités de la Batterie [Wh]")
    plt.ylabel("Chutes de Tension [V]")
    plt.grid()





"""
    Teste de SIMULATION
    ===================
"""

# Simulation(10000, 300000)



"""
    MONTE-CARLO
    ===========
"""

# MonteCarlo(1000, 10000, 100, 0, 1000000, 10000)


"""
    NSGA-II
    =======
"""

NSGA2([1000, 10000], 1000, [0, 100000], 10000, 10, 20)

plt.show()
