import numpy as np
from scipy.optimize import minimize

# Constants
E_MeV = 16.86  # Beam energy
p = np.sqrt(E_MeV**2 - 0.511**2)  # Momentum in MeV/c

# Initial Twiss parameters at entrance
beta_x0 = 6.713
alpha_x0 = -9.771

# Beamline geometry
L_quad = 0.08  # Quadrupole length
L_drift1 = 0.08  # Drift between Q1 and Q2
L_drift2 = 0.08  # Drift between Q2 and Q3
L_drift3 = 2.4 - 2.24059  # Drift from end of Q3 to chicane entrance

# Matrix definitions
def drift(L):
    return np.array([[1, L], [0, 1]])

def quad(L, k):
    if k == 0:
        return drift(L)
    sqrt_k = np.sqrt(abs(k))
    if k > 0:
        return np.array([[np.cos(sqrt_k * L), np.sin(sqrt_k * L) / sqrt_k],
                         [-sqrt_k * np.sin(sqrt_k * L), np.cos(sqrt_k * L)]])
    else:
        return np.array([[np.cosh(sqrt_k * L), np.sinh(sqrt_k * L) / sqrt_k],
                         [sqrt_k * np.sinh(sqrt_k * L), np.cosh(sqrt_k * L)]])

def propagate_twiss(beta, alpha, M):
    gamma = (1 + alpha**2) / beta
    Twiss = np.array([[beta, -alpha], [-alpha, gamma]])
    Twiss_new = M @ Twiss @ M.T
    beta_new = Twiss_new[0, 0]
    alpha_new = -Twiss_new[0, 1]
    return beta_new, alpha_new

# Objective function with symmetry constraint: k1 = k3 = -0.5 * k2
def objective_constrained(k2_array):
    k2 = k2_array[0]
    k1 = -0.5 * k2
    k3 = -0.5 * k2

    M_total = (
        drift(L_drift3) @ quad(L_quad, k3) @ drift(L_drift2) @
        quad(L_quad, k2) @ drift(L_drift1) @ quad(L_quad, k1)
    )

    beta_f, alpha_f = propagate_twiss(beta_x0, alpha_x0, M_total)
    return (beta_f - 0.1)**2 + (alpha_f)**2

# Optimization
initial_guess = [0.0]
bounds = [(-0.1, 0.1)]  # Bounds for k2 only

result = minimize(objective_constrained, initial_guess, bounds=bounds)

# Extract and report results
k2 = result.x[0]
k1 = -0.5 * k2
k3 = -0.5 * k2

M_total = (
    drift(L_drift3) @ quad(L_quad, k3) @ drift(L_drift2) @
    quad(L_quad, k2) @ drift(L_drift1) @ quad(L_quad, k1)
)

beta_final, alpha_final = propagate_twiss(beta_x0, alpha_x0, M_total)

print("Optimal quadrupole strengths (m^-2):")
print(f"  k1 = {k1:.3f}")
print(f"  k2 = {k2:.3f}")
print(f"  k3 = {k3:.3f}")
print("\nMatched Twiss parameters at chicane entrance:")
print(f"  beta_x = {beta_final:.3f} m")
print(f"  alpha_x = {alpha_final:.3f}")