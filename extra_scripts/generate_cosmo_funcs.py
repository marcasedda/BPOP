#!/usr/bin/env python3
"""
Generate redshift_time.txt using cosmological parameters.
Called during build process via Makefile.
"""

import sys
import os
from astropy.cosmology import FlatLambdaCDM
import dd_SFRD
import astropy.units as u
import numpy as np
import jax.numpy as jnp
from tqdm import tqdm

def generate_cosmo_funcs(omega_m, omega_l, h0, output_file1, output_file2, output_file3):
    """
    Generate redshift, age, lookback time, SFRD(z) and dp/dlogZ(Z|z) 
    and save to tables
    
    Parameters:
    -----------
    omega_m : float
        Matter density parameter
    omega_l : float
        Dark energy density parameter
    h0 : float
        Hubble constant in km/s/Mpc
    output_file1 : str
        Path to output file for redshift_time.txt
    output_file2 : str
        Path to output file for sfrd function

    """
    
    # Create cosmology
    cosmo = FlatLambdaCDM(H0=h0*u.km/u.s/u.Mpc, Om0=omega_m)
    
    z_array = np.linspace(1e-6, 30, 10000)
    
    # Compute age of universe at each redshift
    age_array = cosmo.age(z_array).to(u.Gyr).value
    
    # Lookback time
    lookback_array = cosmo.lookback_time(z_array).to(u.Gyr).value
    
    cosmo_correction, cosmo_time_fid, cosmo_time_gen, dL, dV_dz_dtdet, dz_dL = dd_SFRD.cosmo_func(h0, omega_m, jnp.array(z_array))
    zetalog, drhosfr_dlogZ = dd_SFRD.dSFRD_dlogZ(jnp.array(z_array), cosmo_correction, cosmo_time_gen, kind="AM13", sigma_FMR=0.15)
    
    sfrd_array = np.trapz(drhosfr_dlogZ, zetalog, axis = 1)
    sfrd_array[0] = sfrd_array[1]
    dp_dlogZ = drhosfr_dlogZ / sfrd_array[:, None]

    # Write to file
    os.makedirs(os.path.dirname(output_file1), exist_ok=True)
    os.makedirs(os.path.dirname(output_file2), exist_ok=True)
    os.makedirs(os.path.dirname(output_file3), exist_ok=True)

    with open(output_file1, 'w') as f:
        for z, age, lookback in zip(z_array, age_array, lookback_array):
            f.write(f"{z:.5f} {age:.5f} {lookback:.5f}\n")
    
    print(f"Generated {output_file1} with {len(z_array)} entries")
    
    with open(output_file2, 'w') as f:
        for z, sfrd in zip(z_array, sfrd_array):
            f.write(f"{z:.5f} {sfrd:.5f} \n")

    print(f"Generated {output_file2} with {len(z_array)} entries")

    combined = np.column_stack([z_array.reshape(-1, 1), dp_dlogZ])
    np.savetxt(output_file3, combined, fmt='%.5f', delimiter=' ')

    print(f"Generated {output_file3} with {len(z_array)} entries")


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python generate_redshift_table.py <omega_m> <omega_l> <h0> <output_file1> <output_file2> <output_file3>")
        sys.exit(1)
    
    omega_m = float(sys.argv[1])
    omega_l = float(sys.argv[2])
    h0 = float(sys.argv[3])
    output_file1 = sys.argv[4]
    output_file2 = sys.argv[5]
    output_file3 = sys.argv[6]

    generate_cosmo_funcs(omega_m, omega_l, h0, output_file1, output_file2, output_file3)
