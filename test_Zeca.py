import numpy as np
import matplotlib.pyplot as plt
import random


# IMPORTANT: Les unités sont toutes en SI, sauf indication contraire explicite.
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
P_LAC=(Vsst**2)/(4*Req)
P_train = np.copy(P_meca)

for i in range(len(P_meca)):
    if P_meca[i]>0:
        P_train[i]=(P_meca[i]/rendement) + conso
    else:
        P_train[i]=(P_meca[i]*rendement) + conso

P_train_copy=P_train.copy()
V_train=0.5*(Vsst+np.sqrt(Vsst**2-4*P_train*Req))
I_train=(Vsst-V_train)/Req

def Simulation_2(Bat_cap,P_seuil,P_train, temps,vitesse):
    # Bat_cap = 100000 # Capacité de la batterie (Wh).
    Bat_E = np.zeros(len(P_train)) # Énergie contenue dans la batterie (à un instant t).
    Rheo_P = np.zeros(len(P_train)) # Puissance dissipée par le rhéostat.
    P_train_copy = np.copy(P_train)

    for k in range(1, len(P_train)):
        if P_train_copy[k] < 0 and Bat_E[k-1] < Bat_cap: # Le train freine.
            gain_energie = abs(P_train_copy[k])/3600 #dt = 1
            Bat_E[k]= min(Bat_E[k-1]+gain_energie,Bat_cap)
            #P_train[k] = 0
            if Bat_E[k] == Bat_cap: # Dépasse-t-on les limites de la batterie?
                Rheo_P[k] = P_train_copy[k]# On dissipe dans le rhéostat.
            
        elif P_train_copy[k] >= P_seuil and Bat_E[k-1] < Bat_cap: # Si le train demande de l'énergie.
            energie_dispo = Bat_E[k]*3600
            if P_train_copy[k] > energie_dispo:
                P_train_copy[k] -= energie_dispo
                Bat_E[k]=0
            else:
                Bat_E[k]= Bat_E[k-1] - P_train_copy[k]/3600
                P_train_copy[k] = 0
        
        elif vitesse[k] == 0 and Bat_E[k-1] < Bat_cap:
            recharge = P_LAC[k]/3600
            Bat_E[k]= min(Bat_E[k-1]+recharge,Bat_cap)
            if Bat_E[k] > Bat_cap:
                Bat_E[k] = Bat_cap
                Rheo_P[k] = Rheo_P[k-1] + P_train_copy[k] + Bat_E*3600
    
        else:
            Bat_E[k] = Bat_E[k-1] # La batterie conserve son énergie.
            
    V_train = 0.5 * (Vsst + np.sqrt(Vsst**2 - 4 * P_train_copy * Req))
    Ind_qual = Vsst - np.min(V_train)
    
    # plt.figure("Puissance, tension et courant dans le train")
    # # Affichage de la puissance :
    # plt.subplot(3, 1, 1)
    # plt.plot(temps, P_train/1000000, "-k", label="Puissance consommée par le train") # P_train normalisée en MW.
    # plt.title("Puissance, tension et courant dans le train")
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Puissance [MW]")
    # plt.legend()
    # plt.grid()
    
    # # Affichage de la tension:
    # plt.subplot(3, 1, 2)
    # plt.plot(temps, V_train, "-k", label="Tensions aux bornes de la locomotive")
    # plt.legend()
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Tension [V]")
    # plt.grid()
    
    # # Affichage du courant:
    # plt.subplot(3, 1, 3)
    # plt.plot(temps, I_train, "-k", label="Courant traversant le train")
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Courant [A]")
    # plt.legend()
    # plt.grid()
    # plt.show()
    
    return Ind_qual#, P_train_copy, Rheo_P, Bat_E

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
    
    # plt.subplot(3, 1, 2)
    # plt.plot(t, P_train/1000000, "-k", label="Puissance consommée par le train")
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Puissance [MW]")
    # plt.grid()
    # plt.legend()
    
    # plt.subplot(3, 1, 3)
    # plt.plot(t, V_train, "-k", label="Tension aux bornes du train")
    # plt.xlabel("Temps [s]")
    # plt.ylabel("Tension [V]")
    # plt.grid()
    # plt.legend()
    # plt.show()
    
    return Ind_qual

#Simulation_2(0,0,P_train,temps,vitesse)

def Paires_P_seuil_CapaBat():
    Chute_tension = []
    Capa = np.linspace(0,14,1000) #kWh
    P_seuil = np.linspace(0,1000,1000) #kW
    solutions = [(C,P) for C in Capa for P in P_seuil]
    echantillonnage = random.sample(solutions,1000) #Distribution normale de couples
    C = [sol[0] for sol in echantillonnage]
    P = [sol[1] for sol in echantillonnage]
    for sol in echantillonnage:
        Chute_tension.append(Simulation_2(sol[0]*10e3,sol[1]*10e3,P_train,temps,vitesse))
        
    # Visualisation
    
    # plt.figure(figsize=(6, 4))
    # plt.scatter(C, P)  # Scatter plot de l'espace des solutions
    # plt.title("Espace des solutions", fontsize=14)
    # plt.xlabel("Capacité de la batterie (kWh)", fontsize=12)
    # plt.ylabel("Puissance seuil (kW)", fontsize=12)
    # plt.grid(True)
    # plt.show()
    
    plt.figure(figsize=(8, 4))
    
    plt.subplot(221)
    plt.scatter(C, P)  # Scatter plot de l'espace des solutions
    plt.title("Espace des objectifs", fontsize=14)
    plt.xlabel("Capacité de la batterie (kWh)", fontsize=12)
    plt.ylabel("Puissance seuil (kW)", fontsize=12)
    plt.grid(True)
    
    plt.subplot(222)
    plt.scatter(C, Chute_tension)  # Scatter plot de l'espace des solutions
    plt.title("Espace des objectifs", fontsize=14)
    plt.xlabel("Capacité de la batterie (kWh)", fontsize=12)
    plt.ylabel("Chute de tension (V)", fontsize=12)
    plt.grid(True)
    plt.show()

Paires_P_seuil_CapaBat()


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
        Simulation d'un cas particulier, pour une capacité de batterie Bat_cap connue et un seuil de demande d'énergie à la batterie Train_Seuil donnée également.
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
