import numpy as np
import matplotlib.pyplot as plt

# Lire les données depuis le fichier marche.txt
data = np.loadtxt("marche.txt")
temps = data[:, 0]  # Première colonne : temps
position = data[:, 1]  # Deuxième colonne : position


# Calculer la vitesse en utilisant la différence finie
vitesse = np.gradient(position, temps[1]-temps[0])

# Calculer l'accélération en utilisant la différence finie
acceleration = np.gradient(vitesse, temps[1]-temps[0])

#Valeurs fixes données dans l'énoncé
Vsst=790
Rsst=33e-3
rho_lac=131e-6
rho_rail=18e-6
M=70000
A0=780
A1=6.4*10**(-3)
B0=0
B1=0.14*10**(-3)/3.6
C0=0.3634/(3.6**2)
C1=0
g=9.81
rendement=0.8
conso=35000

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

P_train=np.zeros(len(P_meca))

for i in range(len(P_meca)):
    if P_meca[i]>0:
        P_train[i]=(P_meca[i]/rendement) + conso
    else:
        P_train[i]=(P_meca[i]*rendement) + conso
        

#On condidère delta > 0
V_train=0.5*(Vsst+np.sqrt(Vsst**2-4*P_train*Req))
I_train=(Vsst-V_train)/Req


# Afficher les résultats
plt.figure(figsize=(12, 8))

# Graphique de la position
plt.subplot(3, 1, 1)
plt.plot(temps, position, label='Position (x)')
plt.xlabel('Temps (s)')
plt.ylabel('Position (m)')
plt.title('Position en fonction du temps')
plt.legend()
plt.grid()

# Graphique de la vitesse
plt.subplot(3, 1, 2)
plt.plot(temps, vitesse, label='Vitesse (v)', color='orange')
plt.xlabel('Temps (s)')
plt.ylabel('Vitesse (m/s)')
plt.title('Vitesse en fonction du temps')
plt.legend()
plt.grid()

# Graphique de l'accélération
plt.subplot(3, 1, 3)
plt.plot(temps, acceleration, label='Accélération (a)', color='green')
plt.xlabel('Temps (s)')
plt.ylabel('Accélération (m/s²)')
plt.title('Accélération en fonction du temps')
plt.legend()
plt.grid()

#Graphique de la puissance
plt.figure()
plt.plot(temps, P_train, label='Puissance (p)', color='blue')
# plt.plot(temps, P_meca, label='Puissance (p)', color='red')
plt.xlabel('Temps (s)')
plt.ylabel('Puissance (W)')
plt.title('Puissance du train en fonction du temps')
plt.legend()
plt.grid()

# Tension du train au cours du temps
plt.figure()
plt.subplot(2,1,1)
plt.plot(temps, V_train, label='Tension du train')
plt.xlabel('Temps (s)')
plt.ylabel('Tension (V)')
plt.title('Tension du train en fonction du temps')
plt.legend()
plt.grid()

# Courant du train au cours du temps
plt.subplot(2,1,2)
plt.plot(temps, I_train, label='Courant du train')
plt.xlabel('Temps (s)')
plt.ylabel('Courant (A)')
plt.title('Courant du train en fonction du temps')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()
