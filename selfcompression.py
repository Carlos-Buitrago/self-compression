import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import pi, c, e, physical_constants
from scipy.special import gamma, hyp1f1
from scipy.integrate import simpson
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

# Physical constants
r_e = 2.8179403227e-15

# Wiggler parameters
K = 60  # K-value of wriggler 
N_w = 6 # Number of wiggler periods 
lambda_w = 0.25  # Wiggler period (m)
k_w = 2*pi/lambda_w  # Wiggler wave number (1/m)

# Beam parameters - TEX @ 2.4m
I_0 = 22.66  # Peak current (A) calculated for a 30 pC bunch charge and a 158.5e-6 m bunch length (528fs)
I_0 = (I_0/(c * e))
gamma_lorentz = 42.6  # For a 21.77 MeV beam 
sigma_delta = 1.595e-3  # Relative energy spread calculated for a 22 MeV beam with 40 keV spread
sigma_s = 158.5e-6  # Bunch length


def compute_fwhm(x, y):

    y = np.array(y)
    x = np.array(x)

    half_max = np.max(y) / 2.0
    # Find where y crosses half maximum
    indices = np.where(y >= half_max)[0]
    
    if len(indices) < 2:
        raise ValueError("Data does not appear to have a clear peak.")
    
    left_idx = indices[0]
    right_idx = indices[-1]
    
    # Linear interpolation for better precision
    def interp_half(x1, y1, x2, y2):
        return x1 + (half_max - y1) * (x2 - x1) / (y2 - y1)
    
    # Left half-max point
    if left_idx == 0:
        x_left = x[left_idx]
    else:
        x_left = interp_half(x[left_idx-1], y[left_idx-1], x[left_idx], y[left_idx])
    
    # Right half-max point
    if right_idx == len(x) - 1:
        x_right = x[right_idx]
    else:
        x_right = interp_half(x[right_idx], y[right_idx], x[right_idx+1], y[right_idx+1])
    
    return x_right - x_left

# Define the form function
def F(x):
    term1 = 2**(5/6) * gamma(4/3) * hyp1f1(7/6, 3/2, -x**2 / 2) * x
    term2 = 2**(4/3) * gamma(5/6) * hyp1f1(2/3, 1/2, -x**2 / 2)
    return term1 - term2

# Grid definition
s = np.linspace(-5 * sigma_s, 5 * sigma_s, 500)         # Longitudinal position
delta = np.linspace(-200 * sigma_delta, 200 * sigma_delta, 2000)  # Relative energy deviation
S, D = np.meshgrid(s, delta, indexing='ij')

# Define R56_hat and corresponding R56
R56_hat_values = np.linspace(-0.4, 0.4, 1000)
R56_values = R56_hat_values * sigma_s / sigma_delta

# Normalized longitudinal position s_hat
s_hat = s / sigma_s 

# Plot the form function
plt.figure(figsize=(10, 6))
plt.plot(s_hat, F(s_hat), label=r'$F(\hat{s})$')
plt.text(0.25, -3.27, r'$\hat{s}_1 \approx -0.38$', ha='center', va='center', fontsize=10)
plt.text(2.7, 0.5, r'$\hat{s}_2 \approx 2.10$', ha='center', va='center', fontsize=10)
plt.axvline(x=-0.38, color='red', linestyle='--')
plt.axvline(x=2.10, color='red', linestyle='--')
plt.xlabel(r'$\hat{s} = s / \sigma_s$')
plt.ylabel(r'$F(s / \sigma_s)$')    
plt.grid(True)
plt.show()

# Define initial beam distribution f0(s, δ) as function of s and δ
def f0(s, delta):
    #coeff = I_0 / (np.sqrt(2 * pi) * c * e * sigma_delta)
    coeff = I_0 / (np.sqrt(2 * pi) * sigma_delta)
    return coeff * np.exp(-s**2 / (2 * sigma_s**2) - delta**2 / (2 * sigma_delta**2))

# Evaluate f0 over the grid
f0_values = f0(S, D)

# Integrate f0(s, δ) over δ to get the current I(s) (same as lambda(s))
I_s = simpson(f0_values, delta, axis=1)

# Normalize I(s) by I0
#I_s_normalized = I_s / np.max(I_s)  # Change to I_0?
I_s_normalized = I_s / I_0

# Plot the normalized initial current profile
plt.figure(figsize=(10, 6))
plt.plot(s_hat, I_s_normalized, 'r--', linewidth=2, label=r'$I(s)/I_0$')
plt.xlabel(r'$\hat{s} = s/\sigma_s$')
plt.ylabel(r'$I/I_0$')
plt.title('Normalized Current Profile from $f_0$')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Calculate A_hat (normalized amplitude for CSR-induced energy change)
A_hat_numerator = pi * r_e * K * N_w * I_0
#A_hat_denominator = c * e * np.sqrt(2) * gamma_lorentz * sigma_delta * (sigma_s * K * k_w * gamma_lorentz**2)**(1/3)
A_hat_denominator = np.sqrt(2) * gamma_lorentz * sigma_delta * (sigma_s * K * k_w * gamma_lorentz**2)**(1/3)
A_hat = A_hat_numerator / A_hat_denominator
print(f"A_hat: {A_hat:.2f}")

# DEBUGGING - check that maximum energy shift is in delta range
F_vals = F(s / sigma_s)
delta_shift = A_hat * sigma_delta * F_vals
print("Max delta shift: ", np.max(np.abs(delta_shift)))
print("Delta grid covers: ", delta[0], "to", delta[-1])

# Define CSR-induced energy change along the beam (normalized to gamma_lorentz)
def normalized_delta_gamma(s):
    return A_hat * sigma_delta * F(s / sigma_s) 

# Modified beam distribution with CSR-induced energy chirp f1(s, δ)
def f1(s, delta):
    #coeff = I_0 / (np.sqrt(2 * pi) * c * e * sigma_delta)
    coeff = I_0 / (np.sqrt(2 * pi) * sigma_delta)
    return coeff * np.exp(-s**2 / (2 * sigma_s**2) - (delta - normalized_delta_gamma(s))**2 / (2 * sigma_delta**2))

# Initialize matrix to store I(s)/I0 profiles
I_matrix = []
fwhm_values = []
I_max_values = []

# Loop over R56 values
valid_R56_hat_values = []
valid_R56_values = []
for i, R56 in enumerate(R56_values):
    f1_vals = f1(S - R56 * D, D)
    lambda_2 = simpson(f1_vals, delta, axis=1)
    I_over_I0 = lambda_2 / I_0
    try:
        fwhm = compute_fwhm(s, I_over_I0)
        fwhm_values.append(fwhm)
        I_max_values.append(np.max(I_over_I0))  # You can append this here too
        I_matrix.append(I_over_I0)
        valid_R56_hat_values.append(R56_hat_values[i])
        valid_R56_values.append(R56_values[i])
    except ValueError:
        # Optionally print something, or just silently skip
        print(f"  Skipping R56_hat = {R56_hat_values[i]:.4f} due to unclear peak")
        continue


    # Show progress of loop
    if i % 50 == 0:
        print(f"  → Processed {i}/{len(R56_values)} R56 values...")

I_matrix = np.array(I_matrix)  # Shape: (len(R56_hat_values), len(s))

# Create heatmap
plt.figure(figsize=(8, 5))
extent = [s_hat[0], s_hat[-1], valid_R56_hat_values[0], valid_R56_hat_values[-1]]
plt.imshow(I_matrix, aspect='auto', extent=extent, origin='lower', cmap='viridis')
plt.colorbar(label=r'$I(s)/I_0$')
plt.xlabel(r'$\hat{s} = s / \sigma_s$')
plt.ylabel(r'$\hat{R}_{56} = R_{56} \sigma_\delta / \sigma_s$')
plt.tight_layout()
plt.savefig("Heatmap_IoverI0_vs_R56_hat_vs_s_hat.png")
plt.show()

# Plot I_max vs R56_hat
plt.figure(figsize=(10, 6))
plt.plot(valid_R56_hat_values, I_max_values, label="$I_m/I_0$")
plt.xlabel(r"$\hat{R}_{56}$ [m]")
plt.ylabel("I_max$ [A]")
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig("Imax_vs_R56_hat.png")
plt.show()

# Create FWHM plot
fwhm_values_normalized = np.array(fwhm_values) / sigma_s 
plt.figure(figsize=(10, 6))
plt.plot(valid_R56_hat_values, fwhm_values_normalized)
plt.xlabel(r"$\hat{R}_{56}$ [m]")
plt.ylabel("FWHM/$\sigma_s$ [m]")
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig("FWHM_vs_R56_hat.png")
plt.show()


# Create smoothed FWHM plot
fwhm_smoothed = gaussian_filter1d(fwhm_values_normalized, sigma=2)
plt.figure(figsize=(10, 6))
plt.plot(valid_R56_hat_values, fwhm_smoothed)
plt.xlabel(r"$\hat{R}_{56}$ [m]")
plt.ylabel("FWHM/$\sigma_s$ [m]")
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.savefig("FWHMSmooth_vs_R56_hat.png")
plt.show()


# Convert to NumPy arrays
valid_R56_hat_values = np.array(valid_R56_hat_values)
fwhm_smoothed = np.array(fwhm_smoothed)
valid_R56_values = np.array(valid_R56_values)

# Define search windows around visual compression regions
left_mask = (valid_R56_hat_values > -0.2) & (valid_R56_hat_values < -0.1)
right_mask = (valid_R56_hat_values > 0.2) & (valid_R56_hat_values < 0.3)

# Find index of minimum within each window
left_idx = np.argmin(fwhm_smoothed[left_mask])
right_idx = np.argmin(fwhm_smoothed[right_mask])

# Map to global indices
left_true_idx = np.where(left_mask)[0][left_idx]
right_true_idx = np.where(right_mask)[0][right_idx]
top_two = [left_true_idx, right_true_idx]

# Print obtained R56 values
for idx in top_two:
    R56_hat = valid_R56_hat_values[idx]
    R56 = valid_R56_values[idx]
    print(f"  R56_hat for minimum FWHM = {R56_hat:.4f}")
    print(f"  R56 for minimum FWHM = {R56:.4e}")

# 2 current profiles for the 2 chosen R56 values
R56_1 = valid_R56_values[top_two[0]]
f1_vals_1 = f1(S - R56_1 * D, D)
lambda_2_1 = simpson(f1_vals_1, delta, axis=1)
#I_over_I0_1 = lambda_2_1 / np.max(I_s)
I_over_I0_1 = lambda_2_1 / I_0

R56_2 = valid_R56_values[top_two[1]]
f1_vals_2 = f1(S - R56_2 * D, D)
lambda_2_2 = simpson(f1_vals_2, delta, axis=1)
#I_over_I0_2 = lambda_2_2 / np.max(I_s)
I_over_I0_2 = lambda_2_2 / I_0

# Current profile for R56 = 0
f1_vals_3 = f1(S , D)
lambda_2_3 = simpson(f1_vals_3, delta, axis=1)
#I_over_I0_2 = lambda_2_2 / np.max(I_s)
I_over_I0_3 = lambda_2_3 / I_0

# Plot the three current profiles
plt.figure(figsize=(8, 5))
plt.plot(s_hat, I_over_I0_1, label=fr"${{R}}_{{56}} = {R56_values[top_two[0]]:.3f}$")
plt.plot(s_hat, I_over_I0_2, label=fr"${{R}}_{{56}} = {R56_values[top_two[1]]:.3f}$")
plt.plot(s_hat, I_over_I0_3, label=fr"${{R}}_{{56}} = 0$")
plt.xlabel(r"$\hat{s} = s / \sigma_s$")
plt.ylabel(r"$I(s) / I_0$")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("Current_Profiles.png")
plt.show()

# Store results for both R56 minima
peak_currents_R56_1 = []
peak_currents_R56_2 = []

# Current profiles different values of R56
plt.figure(figsize=(10, 6))
for idx in range(0, len(valid_R56_values), 50):
    plt.plot(s_hat, I_matrix[idx])

plt.plot(s_hat, I_over_I0_3, 'k--', label=r'After chirp, no compression ($R_{56} = 0$)', linewidth=2)

plt.xlabel(r"$\hat{s} = s / \sigma_s$")
plt.ylabel(r"$I(s) / I_0$")
plt.title("Current Profiles for Various $R_{56}$ Values")
plt.legend(fontsize=8)
plt.grid(True)
plt.tight_layout()
plt.savefig("Profiles_vsR56.png")
plt.show()
