#!/usr/bin/env python3
"""
Generate redshift_time.txt using cosmological parameters.
Called during build process via Makefile.
"""

import sys
import os
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import numpy as np

def generate_redshift_table(omega_m, omega_l, h0, output_file):
    """
    Generate redshift, age, and lookback time table.
    
    Parameters:
    -----------
    omega_m : float
        Matter density parameter
    omega_l : float
        Dark energy density parameter
    h0 : float
        Hubble constant in km/s/Mpc
    output_file : str
        Path to output file
    """
    
    # Create cosmology
    cosmo = FlatLambdaCDM(H0=h0*u.km/u.s/u.Mpc, Om0=omega_m)
    
    z_array = np.linspace(0, 30, 10000)  
    
    # Compute age of universe at each redshift
    age_array = cosmo.age(z_array).to(u.Gyr).value
    
    # Lookback time = age(z=0) - age(z)
    lookback_array = cosmo.lookback_time(z_array).to(u.Gyr).value
    
    # Write to file
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, 'w') as f:
        for z, age, lookback in zip(z_array, age_array, lookback_array):
            f.write(f"{z:.5f} {age:.5f} {lookback:.5f}\n")
    
    print(f"Generated {output_file} with {len(z_array)} entries")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python generate_redshift_table.py <omega_m> <omega_l> <h0> <output_file>")
        sys.exit(1)
    
    omega_m = float(sys.argv[1])
    omega_l = float(sys.argv[2])
    h0 = float(sys.argv[3])
    output_file = sys.argv[4]
    
    generate_redshift_table(omega_m, omega_l, h0, output_file)
