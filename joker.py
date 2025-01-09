# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 17:33:03 2024

@author: rayan
"""
import numpy as np
import matplotlib.pyplot as plt
def simutrain(Capabatterie,Pseuil):
    
    
    M = 70000 #kg
    A0 = 780 #N
    A1 = 6.4*(10**-3) #N/kg
    B1 = (0.14*(10**-3))/3.6
    C0 = 0.3634/((3.600)**2) #N/(km/s)^2
    g = 9.81


    fichier = "marche.txt"  # Nom du fichier de données
        
    data = np.loadtxt(fichier)  # Charge les données sous forme de tableau numpy
    temps = data[:, 0]  # Première colonne
    position = data[:, 1]  # Deuxième colonne

    vitesse = np.gradient(position,temps[1]-temps[0]) # Vitesse = dPosition/dt
    acceleration= np.gradient(vitesse,temps[1]-temps[0]) # Accélération = dVitesse/dt
    Fresist = A0+A1*M+B1*M*vitesse+(C0)*(vitesse**2) 
    Fmotrice = M*acceleration+Fresist
    Pmeca = Fmotrice*vitesse
    Ptrain = Pmeca

    for i in range (0,len(Pmeca)): # Modélisation de la puissance mécanique du train
        if Pmeca[i] >= 0: 
            Ptrain[i] = (Pmeca[i]/0.8)+35000 # train demande de la puissance
        else: 
            Ptrain[i] = (Pmeca[i]*0.8)+35000 # train donne de l'énergie
        
    # Vérification par affichage 
    # plt.figure()
    # plt.plot(temps[:140],Ptrain[:140])
    # plt.figure()
    # plt.plot(temps[:150],vitesse[:150])



    #### Puissance et tension du train SANS batterie ####

    Vsst = 790 # V
    Rsst = 33*(10**-3) # Ohm
    pLac = 95*(10**-6) # Ohm/m
    pRail = 10*(10**-6) # Ohm/m
    Rlac1 = pLac*position # Ohm
    Rlac2 = pLac*(position[-1]-position) # Ohm
    Rrail1 = pRail*position # Ohm
    Rrail2 = pRail*(position[-1]-position) # Ohm
    Req1 = Rsst + Rlac1 + Rrail1 # Ohm
    Req2 = Rsst + Rlac2 + Rrail2 # Ohm
    rbatterie = 0.9
    Req = (Req1*Req2)/(Req1 + Req2) # Ohm
    I1 = (Vsst)/Req1 # A
    I2 = (Vsst)/Req2 # A
    Itrain = I1+I2 # A

    Psst1 = Vsst*I1 # W
    Psst2 = Vsst*I2 # W

    Vtrain = np.zeros(len(temps)) # Création d'un vecteur de temps

    Vtrain = np.zeros(len(temps)) # Création d'un vecteur de temps
    #### Puissance et tension du train AVEC batterie ####
    Capabatterie = Capabatterie*1e3*3600 # en Wh

    Pbatterie = np.zeros(len(temps)) # Création du vecteur de la puissance de la batterie en fonction du temps

    Plac = np.zeros(len(temps)) # Création du vecteur
    Pseuil = Pseuil*10**3
    Ebatterie = np.zeros(len(temps))
    Ebatterie[0] = Capabatterie
    # Définition du régime du train (accélération ou freinage)
    for i in range(1,len(temps)):
        
        if Ptrain[i] > Pseuil:
            Plac[i] = Pseuil
            Pbatterie[i] = (Ptrain[i]-Plac[i])/rbatterie
            Ebatterie[i] = Ebatterie[i-1]-(Pbatterie[i]+Pbatterie[i-1])*0.5
            if Ebatterie[i]<=0:
                Ebatterie[i] = 0 
                Plac[i] = Ptrain[i]
                Pbatterie[i] = 0
               
        if Ptrain[i] < Pseuil and Ptrain[i] > 0:
            Plac[i] = Ptrain[i]
            Pbatterie[i] = 0
            Ebatterie[i] = Ebatterie[i-1]
            
        if Ptrain[i] <= 0:
            Plac[i] = 0
            if Ebatterie[i-1]>=Capabatterie:
                Pbatterie[i] = 0
                Ebatterie[i] = Capabatterie
            else:
                Pbatterie[i] = Ptrain[i]*rbatterie
                Ebatterie[i] = Ebatterie[i-1]-(Pbatterie[i]+Pbatterie[i-1])*0.5
                if Ebatterie[i]>Capabatterie:
                    Ebatterie[i] = Capabatterie
                
    #Calcul Vtrain           
    for i in range(len(temps)):
        if (Vsst**2)-(4*Req[i]*Plac[i])>0:
            Vtrain[i] = (1/2)*(Vsst+np.sqrt((Vsst**2)-(4*Req[i]*Plac[i])))
                
        else:
            Vtrain[i] = 0   
    
    return max(Vsst-Vtrain)


N = 1000  # Nombre de points
Capabatterie = np.random.uniform(0, 14, N)  # Capacité entre 0.5 et 14 kWh
Pseuil = np.random.uniform(0, 1000, N)  # Puissance seuil entre 0 et 1000 kW
deltaV =  []
for i in range(0,len(Capabatterie)):
    
    deltaV.append(simutrain(Capabatterie[i],Pseuil[i]))
        
deltaV = np.array(deltaV)




def non_dominated_sort(capacite, puissance):
    pareto_points = []
    for i in range(len(capacite)):
        is_dominated = False
        for j in range(len(capacite)):
            if (capacite[j] < capacite[i] and 
                puissance[j] < puissance[i] and 
                (capacite[j] < capacite[i] or puissance[j] > puissance[i])):
                is_dominated = True
                break
        if not is_dominated:
            pareto_points.append(i)
    return np.array(pareto_points)

pareto_indices = non_dominated_sort(Capabatterie, deltaV)

#4. Visualisation des résultats

plt.figure()
# Espace des solutions

plt.scatter(Capabatterie, Pseuil/1000)
plt.scatter(Capabatterie[pareto_indices], Pseuil[pareto_indices]/1000, color='orange')
plt.xlabel("Capacité (kWh)")
plt.ylabel("Puissance seuil (MW)")
plt.title("Espace des solutions")
plt.grid(visible=True)



# Espace des solutions
plt.figure()
plt.scatter(Capabatterie, deltaV)
plt.scatter(Capabatterie[pareto_indices], deltaV[pareto_indices], color='orange')
plt.xlabel("Capacité [kWh]")
plt.ylabel("Différence de tension [V]")
plt.title("Espace des objectifs")
plt.grid(visible=True)

plt.show()

print(simutrain(100000, 350000))