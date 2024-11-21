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

P_LAC_max = (V_SST**2)/(4*min(R_eq))

# Tensions aux bornes du train.
V_train = 0.5*(V_SST+np.sqrt(V_SST**2-4*R_eq*P_LAC))

# Courant tranversant le train.
I_train = P_train/V_train

# Paramètres de la batterie

E_batterie_max = 10000 # Capacité maximale de la batterie (10 kWh)
P_batterie = np.zeros_like(t)  # Puissance échangée avec la batterie
P_rheostat = np.zeros_like(t)  # Puissance dissipée dans le rhéostat
E_batterie = np.zeros_like(t)  # Énergie stockée dans la batterie
dt = t[1]-t[0]
seuil = 0.85 * np.max(P_train)  # Seuil de puissance pour décharger la batterie
# Gestion de la batterie
for i in range(1,len(t)):
    if P_train[i] < 0:  # Train en freinage
        if E_batterie[i - 1] < E_batterie_max: #Energie batterie inférieur à la capacité max
            P_batterie[i] = P_train[i]
        else:
            P_batterie[i] = 0
            P_rheostat[i] = -P_train[i] - P_batterie[i]
    if P_train[i] > 0:  # Train en consommation
        if E_batterie[i - 1] > 0:  # Batterie dispo
            P_batterie[i] = min(P_train[i], E_batterie[i - 1]/dt)
            if P_train[i] > P_batterie[i]: #Puissance insufisante
                P_train[i] -= P_batterie[i] #Train prend la puissance dispo
                E_batterie[i] = 0 #vide la batterie
            else:
                P_batterie[i] = P_batterie[i-1] - P_train[i] #On vide la batterie en fonction de la puissance demandée par le train
        else:
            P_batterie[i] = 0
            P_rheostat[i] = 0
    # Mise à jour de l'énergie de la batterie
    E_batterie[i] = E_batterie[i - 1] - P_batterie[i]*3600*dt

plt.figure(figsize=(12, 6))
plt.plot(t, E_batterie, label='Charge de la batterie')
plt.xlabel('Temps (s)')
plt.ylabel('Charge de la batterie (kWh)')
plt.title('Évolution de la charge de la batterie en fonction du temps')
plt.legend()
plt.grid()
plt.show()
