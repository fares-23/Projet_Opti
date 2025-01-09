import numpy as np
import matplotlib.pyplot as plt

# --- Étape 1 : Fonction Monte Carlo pour générer des données ---
def generate_solutions(n_samples, capacity_min, capacity_max, voltage_min, voltage_max):
    """
    Génère des solutions simulées pour un problème bi-objectifs.
    
    :param n_samples: Nombre de solutions à générer
    :param capacity_min: Capacité minimale de la batterie (kWh)
    :param capacity_max: Capacité maximale de la batterie (kWh)
    :param voltage_min: Chute de tension minimale (V)
    :param voltage_max: Chute de tension maximale (V)
    :return: Tableau numpy des solutions (capacités, chutes de tension)
    """
    capacities = np.random.uniform(capacity_min, capacity_max, n_samples)
    voltage_drops = np.random.uniform(voltage_min, voltage_max, n_samples)
    return np.column_stack((capacities, voltage_drops))

# --- Étape 2 : Fonction pour trouver les solutions du front de Pareto ---
def non_dominant_sort(solutions):
    """
    Identifie les solutions non dominées (front de Pareto) parmi un ensemble de solutions.
    
    :param solutions: Tableau numpy (2 colonnes : objectifs 1 et 2)
    :return: Tableau numpy des solutions du front de Pareto
    """
    pareto_front = []
    for i, sol_i in enumerate(solutions):
        dominated = False
        for j, sol_j in enumerate(solutions):
            if all(sol_j <= sol_i) and any(sol_j < sol_i):
                dominated = True
                break
        if not dominated:
            pareto_front.append(sol_i)
    return np.array(pareto_front)

# --- Étape 3 : Générer les solutions et identifier le front de Pareto ---
# Paramètres pour la simulation
n_samples = 1000
capacity_min, capacity_max = 100, 15000  # Plages de capacité batterie (kWh)
voltage_min, voltage_max = 50, 300       # Plages de chute de tension maximale (V)

# Génération des solutions
solutions = generate_solutions(n_samples, capacity_min, capacity_max, voltage_min, voltage_max)

# Identification du front de Pareto
pareto_front = non_dominant_sort(solutions)

# Identification de la meilleure solution (exemple : minimisation de la somme des deux objectifs)
best_solution = pareto_front[np.argmin(np.sum(pareto_front, axis=1))]

# --- Étape 4 : Visualisation ---
plt.figure(figsize=(12, 6))

# Points Monte Carlo
plt.scatter(solutions[:, 0], solutions[:, 1], alpha=0.3, label="Solutions Monte Carlo", c="blue")

# Points du front de Pareto
plt.scatter(pareto_front[:, 0], pareto_front[:, 1], color="red", label="Front de Pareto")

# Meilleure solution
plt.scatter(best_solution[0], best_solution[1], color="green", label="Meilleure solution", s=100)

# Personnalisation du graphique
plt.xlabel("Capacité de la batterie (kWh)")
plt.ylabel("Chute de tension maximale (V)")
plt.title("Optimisation bi-objectifs : Capacité vs Chute de tension maximale")
plt.legend()
plt.grid()

# Affichage
plt.show()
