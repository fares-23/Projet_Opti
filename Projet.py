import numpy as np
import matplotlib.pyplot as plt


# IMPORTANT: Les unités sont toutes en SI.


"""
    VARIABLES
    =========
"""

M = 70000000 # Masse du train.
A_0 = 780
A_1 = 6.4*10**(-6)
B_0 = 0
B_1 = 0.14*3600*10**(-9)
C_0 = 0.3634*3600**2*10**(-6)
C_1 = 0
alpha = 0 # Pente que le train doit gravir (compté positivement).
g = 9.81 # Accélération de la pesanteur.
V_SST = 790 # Tension délivrée par la sous-station.
R_SST = 33*10**(-3) # Résistance interne de la sous-station.
rho_LAC = 131*10**(-6) # Résistance linéique de la LAC.
rho_rail = 19*10**(-6) # Résistance linéique des rails.
x_Tot = 5000 # Distance totale que le train doit parcourir durant son trajet.
SysBor = 35000 # Consommation du système de bord.


"""
    GRANDEURS CONNUES
    =================
"""

marche = np.loadtxt("marche.txt")
t = marche[:, 0] # Temps.
x = marche[:, 1] # Distance parcourue par le train.


"""
    GRANDEURS CALCULÉES
    ===================
"""

# Vitesse v.
v = np.zeros(len(x))
for k in range(len(x)-1):
    v[k] = x[k+1]-x[k]
v[-1] = v[-2]

# Accélération a.
a = np.zeros(len(v))
for k in range(len(v)-1):
    a[k] = v[k+1]-v[k]
a[-1] = a[-2]

# F_resistive.
F_resistive = (A_0+A_1*M)+(B_0+B_1*M)*v+(C_0+C_1*M)*v**2

# F_motrice.
F_motrice = M*a+M*g*np.sin(alpha)+F_resistive

# P_mecanique.
P_mecanique = F_motrice*v

# Puissance totale consommée par le train: P_Tot.
P_Tot = P_mecanique + SysBor


"""
    AFFICHAGE
    =========
"""

plt.figure(1)
plt.plot(t, P_Tot*10**(-6), "-b", label="Puissance mécanique consommé par le train")
plt.legend()
plt.xlabel("Temps [s]")
plt.ylabel("Puissance [MW]")
plt.grid()
plt.show()
