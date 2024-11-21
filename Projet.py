import numpy as np
import matplotlib.pyplot as plt


# IMPORTANT: Les unités sont toutes en SI.


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
V_SST = 790 # Tension délivrée par la sous-station.
R_SST = 33*10**(-3) # Résistance interne de la sous-station.
rho_LAC = 131*10**(-6) # Résistance linéique de la LAC.
rho_rail = 18*10**(-6) # Résistance linéique des rails.
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
    TRAIN & RÉSEAU
    ==============
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


"""
    BATTERIE
    ========
"""


# Si appelé, simule le fonctionnement de la batterie.
Bat_cap = 15000 # Capacité de la batterie (Wh).
Bat_E = np.zeros(len(P_train)) # Énergie contenue dans la batterie.
Bat_Charge = np.zeros(len(P_train)) # Charge de la batterie, au cours du temps.
Bat_Charge[0] = 0 # Charge de départ de la batterie.
if Bat_Charge[0] > Bat_cap: # On vérifie qu'on ne dépasse pas le capacité de la batterie.
    raise ValueError("On dépasse les capacités de la batterie!")
Train_Seuil = 250000 # Seuil de puissance (négatif) auquel le train demande de l'énergie.
Rheo_P = np.zeros(len(P_train)) # Puissance dissipée par le rhéostat.
Rheo_P[0] = 0 # Le rhéostat ne dissipe rien au départ.

for k in range(1, len(t)):
    if P_train[k-1] < 0 and Bat_E[k-1] < Bat_cap: # Le train freine et la batterie n'est pas pleine.
        Bat_E[k] = Bat_E[k-1]-P_train[k-1]
        P_train[k] = 0
        if Bat_E[k] > Bat_cap: # On dépasse la capacité de la batterie?
            Diff = Bat_E[k] - Bat_cap
            Bat_E[k] = Bat_cap
            Rheo_P[k] = Diff # On dissipe dans le rhéostat.
    if P_train[k] >= Train_Seuil and Bat_E[k-1] > 0: # Si le train demande de l'énergie et que la batterie n'est pas vide.
        Bat_E[k] = Bat_E[k-1] - (P_train[k]-Train_Seuil)
        P_train[k] = Train_Seuil
        if Bat_E[k] < 0: # Dépasse-t-on les limites de la batterie?
            Diff = -Bat_E[k]
            Bat_E[k] = 0
            P_train[k] += Diff

# On recalcule ce qui a été modifié.
P_LAC = np.zeros(len(P_train))
for k in range(len(P_train)):
    if P_train[k] > (V_SST**2)/(4*R_eq[k]):
        P_LAC[k] = (V_SST**2)/(4*R_eq[k])-SysBor
    else:
        P_LAC[k] = P_train[k]


"""
    AFFICHAGE
    =========
"""


plt.figure("X,V,T du train")
# Affichage position :
plt.subplot(3, 1, 1)
plt.plot(t, x/1000, "-k", label="Position du train") # position normalisé en km
plt.title("Position, vitesse, accélération du train en fonction du temps")
plt.xlabel("Temps [s]")
plt.ylabel("Longueur [km]")
plt.grid()
plt.legend()

# Affichage vitesse :
plt.subplot(3, 1, 2)
plt.plot(t, v/1000, "-k", label="Vitesse du train") # vitesse normalisé en km/s
plt.xlabel("Temps [s]")
plt.ylabel("Vitesse [km/s]")
plt.grid()
plt.legend()
# Affichage accélération :
plt.subplot(3, 1, 3)
plt.plot(t, a/g, "-k", label="Accélération du train") # accélération normalisé en g
plt.xlabel("Temps [s]")
plt.ylabel("Accélération [g]")
plt.grid()
plt.legend()


plt.figure("P,V,I du train")
# Affichage de la puissance :
plt.subplot(3, 1, 1)
plt.plot(t, P_train/1000000, "-k", label="Puissance consommée") # P_train normalisé en MW
plt.title("Puissance, tension et courant dans le train")
plt.xlabel("Temps [s]")
plt.ylabel("Puissance [MW]")
plt.legend()
plt.grid()

# Affichage de la tension :
plt.subplot(3, 1, 2)
plt.plot(t, V_train, "-k", label="Tensions aux bornes de la locomotive")
plt.legend()
plt.xlabel("Temps [s]")
plt.ylabel("Tension [V]")
plt.grid()

# Affichage du courant :
plt.subplot(3, 1, 3)
plt.plot(t, I_train, "-k", label="Courant traversant le train")
plt.xlabel("Temps [s]")
plt.ylabel("Courant [A]")
plt.legend()
plt.grid()


plt.figure("I_1, I_2, I_train")
# Affichage de I_1, I_2 et I_train:
plt.plot(t, I_1, "-b", label="I_1")
plt.plot(t, I_2, "-g", label="I_2")
plt.plot(t, I_train, "-k", label="I_train")
plt.xlabel("Temps [s]")
plt.ylabel("Courant [A]")
plt.title("Courants dans les branches du circuit")
plt.grid()
plt.legend()


plt.figure("P_SS1, P_SS2, P_train")
# Affichage Puissance des stations et du train
plt.subplot(3, 1, 1)
plt.plot(t, P_SST1/1000000, "-b", label="Sous-station 1") # normalisé en MW
plt.plot(t, P_SST2/1000000, "-g", label="Sous-station 2") # normalisé en MW
plt.xlabel("Temps [s]")
plt.ylabel("Puissance [MW]")
plt.title("Puissances mises en jeu")
plt.legend()
plt.grid()

# Affichage de la puissance des stations avec pertes
plt.subplot(3, 1, 2)
plt.plot(t, P_SST1_loss/1000000, "-m", label="Perte 1") # normalisé en MW
plt.plot(t, P_SST2_loss/1000000, "-c", label="Perte 2") # normalisé en MW
plt.grid()
plt.legend()
plt.xlabel("Temps [s]")
plt.ylabel("Puissance [MW]")

# Affichage de la puissance des stations avec pertes
plt.subplot(3, 1, 3)
plt.plot(t, (P_SST1+P_SST2-P_SST1_loss-P_SST2_loss)/1000000, "-r", label="Puissance des stations, avec pertes") # normalisé en MW
plt.plot(t, P_train/1000000, "-k", label="Puissance consommée") # P_train normalisé en MW
plt.grid()
plt.legend()
plt.xlabel("Temps [s]")
plt.ylabel("Puissance [MW]")

# Affichage de l'énergie de la batterie, de la consommation du rhéostat et de l'énergie du train ainsi que de sa tension.
plt.figure()
plt.subplot(3, 1, 1)
plt.plot(t, Bat_E/1000000, "-b", label="Énergie dans la batterie")
plt.plot(t, Rheo_P/1000000, "-g", label="Somme des énergies perdues dans le rhéostat")
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


plt.show()
