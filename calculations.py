import numpy as np

# Constants
G = 6.674*10**(-11) # Graviational constant in m^3 /(kg s^2)
Msol = 1.988*10**(30) # Solar mass in kg

# Conversion factors
kpc_to_m = 3.086e19  # kpc is 3.086e+19 m
pc_to_m = 3.086e16 # pc is 3.086e+16 m
Mpc_to_m = 1e3*kpc_to_m

# IC parameters
Omega_m = 0.3077
H0 = 67.81  # in km/s/kpc
L = 20.0  # in Mpc
N_particles = 64**3  # 64^3 particles

# Convert H0 to  m /s / kpc
H0_km_kpc = H0 *1000

# Convert H0 to s^-1
H0_s = H0_km_kpc/kpc_to_m  

# Calculate critical density in SI units
rho_c = (3 * H0_s**2) / (8 * np.pi * G)  

# Calculate total mass of the simulation box in SI units
V = (L*Mpc_to_m)**3  # Volume in m^3

# Mass in kg
M_total = Omega_m * rho_c * V

# Calculate particle mass
M_particle = M_total / N_particles  # in kg

M_particle_Msol = M_particle/Msol

print(f"Particle mass: {M_particle_Msol:.4e} Msol")

# Number of particles needed to resolve a halo
N_particles_per_halo = 32

# Minimum halo mass
M_min_halo = N_particles_per_halo * M_particle # in kg

M_min_halo_Msol = M_min_halo/Msol

print(f"Minimum halo mass that the simulation can resolve: {M_min_halo_Msol:.4e} Msol")



