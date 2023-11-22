import numpy as np

# Constants
H0_code = np.sqrt(8 * np.pi / 3)  # Hubble constant in code units
L_code = 1  # Box size in code units
L_phys = 20000  # Box size in kpc

# Conversion factor for length from code to physical units
dKpcUnit = L_phys / L_code

# Gravitational constant in physical units km2 Mpc MSun-1 s-2
G_phys = 4.301e-9

# G_phys from https://lweb.cfa.harvard.edu/~dfabricant/huchra/ay145/constants.html

# Hubble constant in physical units
H0_phys = 67.81 # in km/s/kpc
#For H0 these are the units in the assignment, however I think they maybe are supposed to be in km/s/Mpc, as this is the closest I get to the correct answer.

# Critical density
rho_crit_phys = (3 * H0_phys**2 )/ (8 * np.pi * G_phys)

# Mass unit conversion from code to solar masses, converting back to Mpc
dMsolUnit = ((rho_crit_phys) * (dKpcUnit/1000)**3)

print(f'{dMsolUnit:.6g}')
"""
Output:
1.02091e+15
"""
