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
rho_rail = 18*10**(-6) # Résistance linéique des rails.
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
v = np.gradient(x, t[1]-t[0])

# Accélération a.
a = np.gradient(v, t[1]-t[0])

#F_resistive.
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

plt.figure("données initiales")

# Graphique de la position
plt.subplot(3, 1, 1)
plt.plot(t, x, label='Position (x)')
plt.xlabel('Temps (s)')
plt.ylabel('Position (m)')
plt.title('Position en fonction du temps')
plt.legend()
plt.grid()

# Graphique de la vitesse
plt.subplot(3, 1, 2)
plt.plot(t, v, label='Vitesse (v)', color='orange')
plt.xlabel('Temps (s)')
plt.ylabel('Vitesse (m/s)')
plt.title('Vitesse en fonction du temps')
plt.legend()
plt.grid()

# Graphique de l'accélération
plt.subplot(3, 1, 3)
plt.plot(t, a, label='Accélération (a)', color='green')
plt.xlabel('Temps (s)')
plt.ylabel('Accélération (m/s²)')
plt.title('Accélération en fonction du temps')
plt.legend()
plt.grid()


plt.figure("Puissance mécanique consommé par le train")
plt.plot(t, P_Tot*10**(-6), "-b", label="Puissance mécanique consommé par le train")
plt.legend()
plt.xlabel("Temps [s]")
plt.ylabel("Puissance [MW]")
plt.grid()


plt.show()