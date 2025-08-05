
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi, c, e, physical_constants
from mpl_toolkits.mplot3d import Axes3D

r_e = 2.8179403227e-15

# Beam parameters - TEX @ 2.4m
I_0 = 22.66  # Peak current (A) calculated for a 30 pC bunch charge and a 158.5e-6 m bunch length (528fs)
gamma_lorentz = 42.6  # For a 21.77 MeV beam 
sigma_delta = 1.595e-3  # Relative energy spread calculated for a 22 MeV beam with 40 keV spread
sigma_s = 158.5e-6  # Bunch length

# Wiggler parameters
K_vals = np.linspace(10.0, 100.0, 50)
lambda_w_vals = np.linspace(0.1, 0.5, 25)
N_w_vals = np.array([4, 5, 6, 7, 8, 9, 10])  # Number of wiggler periods

# Create 3D meshgrid
K_grid, lambda_w_grid, N_w_grid = np.meshgrid(K_vals, lambda_w_vals, N_w_vals, indexing='ij')

# Calculate k_w for each lambda_w
k_w_grid = 2 * pi / lambda_w_grid

# Calculate A_hat
A_hat_numerator = pi * r_e * K_grid * N_w_grid * I_0
#A_hat_denominator = np.sqrt(2) * gamma_lorentz * sigma_delta * (sigma_s * K_grid * k_w_grid * gamma_lorentz**2) ** (1/3)
A_hat_denominator = c * e * np.sqrt(2) * gamma_lorentz * sigma_delta * (sigma_s * K_grid * k_w_grid * gamma_lorentz**2) ** (1/3)
A_hat = A_hat_numerator / A_hat_denominator

# Flatten for scatter plot
K_flat = K_grid.flatten()
lambda_w_flat = lambda_w_grid.flatten()
N_w_flat = N_w_grid.flatten()
A_hat_flat = A_hat.flatten()

# Plot
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(K_flat, lambda_w_flat, N_w_flat, c=A_hat_flat, cmap='viridis')
ax.set_xlabel('K')
ax.set_ylabel('$λ_w$ (m)')
ax.set_zlabel('$N_w$')
ax.set_title('3D Scatter of Â(K, $λ_w$, $N_w$)')
fig.colorbar(sc, label='Â')
plt.savefig("A_vs_wiggler_3D.png")
plt.show()

# Target Â value and tolerance
A_target = 2  
tolerance = 0.1 

# Maximum wiggler length
L_max = 1.5 # meters

# Maximum dipole strength
B_max = 2.5  # Tesla

# Find indices where A_hat is close to target and total length is within maximum wiggler length
indices = np.where(np.abs(A_hat_flat - A_target) <= tolerance)[0]
indices = indices[np.where(lambda_w_flat[indices] * N_w_flat[indices] <= L_max)[0]]
indices = indices[np.where(K_flat[indices] / (93.4 * lambda_w_flat[indices]) <= B_max)[0]]

# Print matching parameter sets
print(f"Parameter sets where Â ≈ {A_target} (±{tolerance}):")
for idx in indices:
    print(f"K = {K_flat[idx]:.3f}, λ_w = {lambda_w_flat[idx]:.3f} m, N_w = {int(N_w_flat[idx])}, Â = {A_hat_flat[idx]:.4f}, Total_Length = {lambda_w_flat[idx] * N_w_flat[idx]:.3f} m, B of dipoles = {K_flat[idx]/(93.4 * lambda_w_flat[idx]):.3f} T")
