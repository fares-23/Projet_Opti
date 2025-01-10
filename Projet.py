import numpy as np
import matplotlib.pyplot as plt
import random


# IMPORTANT: Les unités sont toutes en SI, sauf indication contraire explicite.



def Simulation(Bat_cap,seuil, PlotSim=False):
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

    """
        AFFICHAGE
        =========
    """
    if PlotSim:
        plt.figure("Position, vitesse, accélération du train en fonction du temps")
        # Affichage position :
        plt.subplot(3, 1, 1)
        plt.plot(t, x/1000, "-k", label="Position du train") # Position normalisée en km.
        plt.title("Position, vitesse, accélération du train en fonction du temps")
        plt.xlabel("Temps [s]")
        plt.ylabel("Longueur [km]")
        plt.grid()
        plt.legend()

        # Affichage vitesse :
        plt.subplot(3, 1, 2)
        plt.plot(t, v/1000, "-k", label="Vitesse du train") # Vitesse normalisée en km/s.
        plt.xlabel("Temps [s]")
        plt.ylabel("Vitesse [km/s]")
        plt.grid()
        plt.legend()
        # Affichage accélération :
        plt.subplot(3, 1, 3)
        plt.plot(t, a/g, "-k", label="Accélération du train") # Accélération normalisée en g.
        plt.xlabel("Temps [s]")
        plt.ylabel("Accélération [g]")
        plt.grid()
        plt.legend()


        plt.figure("Puissance, tension et courant dans le train")
        # Affichage de la puissance :
        plt.subplot(3, 1, 1)
        plt.plot(t, P_train/1000000, "-k", label="Puissance consommée par le train") # P_train normalisée en MW.
        plt.title("Puissance, tension et courant dans le train")
        plt.xlabel("Temps [s]")
        plt.ylabel("Puissance [MW]")
        plt.legend()
        plt.grid()

        # Affichage de la tension:
        plt.subplot(3, 1, 2)
        plt.plot(t, V_train, "-k", label="Tensions aux bornes de la locomotive")
        plt.legend()
        plt.xlabel("Temps [s]")
        plt.ylabel("Tension [V]")
        plt.grid()

        # Affichage du courant:
        plt.subplot(3, 1, 3)
        plt.plot(t, I_train, "-k", label="Courant traversant le train")
        plt.xlabel("Temps [s]")
        plt.ylabel("Courant [A]")
        plt.legend()
        plt.grid()


        plt.figure("Courants dans les branches du circuit")
        # Affichage de I_1, I_2 et I_train:
        plt.plot(t, I_1, "-b", label="I_1")
        plt.plot(t, I_2, "-g", label="I_2")
        plt.plot(t, I_train, "-k", label="I_train")
        plt.xlabel("Temps [s]")
        plt.ylabel("Courant [A]")
        plt.title("Courants dans les branches du circuit")
        plt.grid()
        plt.legend()


        plt.figure("Puissances mises en jeu")
        # Affichage Puissance des stations et du train.
        plt.subplot(3, 1, 1)
        plt.plot(t, P_SST1/1000000, "-b", label="Sous-station 1") # normalisé en MW
        plt.plot(t, P_SST2/1000000, "-g", label="Sous-station 2") # normalisé en MW
        plt.xlabel("Temps [s]")
        plt.ylabel("Puissance [MW]")
        plt.title("Puissances mises en jeu")
        plt.legend()
        plt.grid()

        # Affichage de la puissance des stations avec pertes.
        plt.subplot(3, 1, 2)
        plt.plot(t, P_SST1_loss/1000000, "-m", label="Perte 1") # normalisé en MW
        plt.plot(t, P_SST2_loss/1000000, "-c", label="Perte 2") # normalisé en MW
        plt.grid()
        plt.legend()
        plt.xlabel("Temps [s]")
        plt.ylabel("Puissance [MW]")

        # Affichage de la puissance des stations avec pertes.
        plt.subplot(3, 1, 3)
        plt.plot(t, (P_SST1+P_SST2-P_SST1_loss-P_SST2_loss)/1000000, "-r", label="Puissance des stations, avec pertes") # normalisé en MW
        plt.plot(t, P_train/1000000, "-k", label="Puissance consommée par le train") # P_train normalisé en MW
        plt.grid()
        plt.legend()
        plt.xlabel("Temps [s]")
        plt.ylabel("Puissance [MW]")


        plt.figure("Gestion de la Batterie")
        # Affichage de l'énergie de la batterie, de la consommation du rhéostat et de l'énergie du train ainsi que de sa tension.
        plt.subplot(3, 1, 1)
        plt.plot(t, Bat_E/1000000, "-b", label="Énergie dans la batterie")
        plt.plot(t, Rheo_P/1000000, "-g", label="Énergie perdue dans le rhéostat")
        plt.title("Gestion de la Batterie")
        plt.xlabel("Temps [s]")
        plt.ylabel("Énergie [MJ]")
        plt.grid()
        plt.legend()

        plt.subplot(3, 1, 2)
        plt.plot(t, P_train/1000000, "-k", label="Puissance consommée par le train")
        plt.xlabel("Temps [s]")
        plt.ylabel("Puissance [MW]")
        plt.grid()
        plt.legend()

        plt.subplot(3, 1, 3)
        plt.plot(t, V_train, "-k", label="Tension aux bornes du train")
        plt.xlabel("Temps [s]")
        plt.ylabel("Tension [V]")
        plt.grid()
        plt.legend()


        plt.figure("Indication de la Qualité du Système")
        # Comparaison, dans un histogramme, de Ind_qual et IndCrit.
        plt.bar(["Indicateur Critique\n(à ne pas dépasser)", "Indicateur Actuelle"], [IndCrit, Ind_qual], color="grey", edgecolor="black")
        plt.title("Indication de la Qualité du Système")
        plt.ylabel("Différence de Tension [V]")
        plt.grid()
        plt.legend()



    return Ind_qual # Indicateur de qualité = Chute de tension maximale.



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

def choisir_meilleur_point(pareto_front, poids_capacite, poids_chute):
    scores = []
    for point in pareto_front:
        capacite, _ , chute = point
        score = poids_capacite * capacite + poids_chute * chute
        scores.append(score)

    meilleur_indice = np.argmin(scores)
    return pareto_front[meilleur_indice]



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
    print_label = False
    for pop in process[:-1]:
        for individual in pop:
            plt.plot(individual[0]/1000, individual[1]/1000, "+k")
    for individual in process[-1]:
        if not print_label:
            plt.plot(individual[0]/1000, individual[1]/1000, "+r", label="Dernière génération")
            print_label = True
        else:
            plt.plot(individual[0]/1000, individual[1]/1000, "+r")
    plt.title("Espace des Solutions / Espace des Objectifs")
    plt.xlabel("Capacités de la Batterie [kWh]")
    plt.ylabel("Seuils [kW]")
    plt.legend()
    plt.grid()
    plt.subplot(2, 1, 2)
    print_label = False
    for pop in process[:-1]:
        for individual in pop:
            plt.plot(individual[0]/1000, individual[2], "+k")
    for individual in process[-1]:
        if not print_label:
            plt.plot(individual[0]/1000, individual[2], "+r", label="Dernière génération")
            print_label = True
        else:
            plt.plot(individual[0]/1000, individual[2], "+r")
    plt.xlabel("Capacités de la Batterie [kWh]")
    plt.ylabel("Chutes de Tension [V]")
    plt.legend()
    plt.grid()





"""
    Teste de SIMULATION
    ===================
"""

Simulation(10000, 300000, True)



"""
    MONTE-CARLO
    ===========
"""

solutions = MonteCarlo()
sol = non_dominant_sort(solutions)
non_dominés = sol[0]

capacites = [solu[0] for solu in solutions]
seuils = [solu[1] for solu in solutions]
chutes = [solu[2] for solu in solutions]

pareto_capacites = [solu[0] for solu in non_dominés]
pareto_chutes = [solu[2] for solu in non_dominés]
pareto_seuils = [solu[1] for solu in non_dominés]

best_capacite,best_seuil,best_chute = choisir_meilleur_point(non_dominés, 1, 17)

fig, axs = plt.subplots(1, 2, figsize=(12, 5))  # Create a figure with 2 subplots side by side
plt.get_current_fig_manager().set_window_title('Monte-Carlo') # Nom de la fenêtre



# First subplot: Capacité vs Chute de tension maximale
axs[0].scatter(capacites, chutes, alpha=0.3, label="Solutions Monte Carlo")
axs[0].scatter(pareto_capacites, pareto_chutes, color="red", label="Front de Pareto")
axs[0].scatter(best_capacite, best_chute, color="green", label="Meilleures solutions")
axs[0].set_xlabel("Capacité de la batterie (Wh)")
axs[0].set_ylabel("Chute de tension maximale (V)")
axs[0].set_title("Espace des objectifs")
axs[0].legend()
axs[0].grid()

# Second subplot: Capacité vs Puissance Seuil
axs[1].scatter(capacites, seuils, alpha=0.3, label="Solutions Monte Carlo")
axs[1].scatter(pareto_capacites, pareto_seuils, color="red", label="Front de Pareto")
axs[1].scatter(best_capacite, best_seuil, color="green", label="Meilleures solutions")
axs[1].set_xlabel("Capacité de la batterie (Wh)")
axs[1].set_ylabel("Puissance Seuil (W)")
axs[1].set_title("Espace des solutions")
# Uncomment if legend is desired
# axs[1].legend()
axs[1].grid()



"""
    NSGA-II
    =======
"""

# NSGA2([1000, 15000], 1000, [0, 350000], 10000, 10, 25)


plt.show()
