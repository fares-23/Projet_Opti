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
