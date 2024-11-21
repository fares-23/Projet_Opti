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
rend = 0.8 # Rendement du moteur+convertisseurs.


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

# Résistances du circtui.
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
# V_train = 0.5*(V_SST+np.sqrt(V_SST**2-4*R_eq*P_LAC))


"""
    AFFICHAGE
    =========
"""

# plt.figure()
# plt.subplot(4, 1, 1)
# plt.plot(t, x/1000, "-k", label="Position du train")
# plt.title("Position du train en fonction du temps")
# plt.xlabel("Temps [s]")
# plt.ylabel("Longueur [km]")
# plt.grid()
# plt.legend()
#
# plt.subplot(4, 1, 2)
# plt.plot(t, v/1000, "-k", label="Vitesse du train")
# plt.title("Vutesse du train en fonction du temps")
# plt.xlabel("Temps [s]")
# plt.ylabel("Vitesse [km/s]")
# plt.grid()
# plt.legend()
#
# plt.subplot(4, 1, 3)
# plt.plot(t, a/g, "-k", label="Accélération du train")
# plt.title("Accélération du train en fonction du temps")
# plt.xlabel("Temps [s]")
# plt.ylabel("Accélération [g]")
# plt.grid()
# plt.legend()
#
# plt.subplot(4, 1, 4)
# plt.plot(t, P_train/1000000, "-b", label="Puissance consommée")
# plt.legend()
# plt.xlabel("Temps [s]")
# plt.ylabel("Puissance [MW]")
# plt.title("Puissance consommée par le train en fonction du temps")
# plt.grid()
# plt.show()


plt.figure()
plt.plot(t, (V_SST**2-4*R_eq*P_LAC), "-b", label="Truc")
plt.legend()
plt.xlabel("Temps [s]")
plt.ylabel("SI")
plt.grid()


v_train = 0.5*(V_SST+np.sqrt(V_SST**2-4*R_eq*P_LAC))
plt.figure()
plt.plot(t, v_train, "-b", label="Truc")
plt.legend()
plt.xlabel("Temps [s]")
plt.ylabel("SI")
plt.grid()
plt.show()

# plt.figure()
# plt.plot(t, V_train, "-b", label="Tensions aux bornes de la locomotive")
# plt.legend()
# plt.xlabel("Temps [s]")
# plt.ylabel("Tension [V]")
# plt.grid()
# plt.show()
