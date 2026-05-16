import jax.numpy as jnp
from jax import jit, vmap
from functools import partial

c = 299792.458 #km/s
s_Gyr = 977.792 #second to Gyr conversion
nzeta=100
zetalog= jnp.linspace(-4., -1., nzeta)

#LambdaCDM functions with jit
@jit
def cosmo_func(h, Om, z):
    #We approximately not consider Omega rad.
    H = h * jnp.sqrt(Om * (1+z)**3 + 1-Om)
    H_fid = 70. * jnp.sqrt(0.30 * (1+z)**3 + 0.7)

    base = jnp.linspace(0, 1, 1000)  
    z_mat = z[:, None] * base  
    fid_H_mat = 70. * jnp.sqrt(0.30 * (1+z_mat)**3 + 1 - 0.30)
    gen_H_mat = h * jnp.sqrt(Om * (1+z_mat)**3 + 1 - Om)
    gen_integrand = c / gen_H_mat
    fid_integrand = c / fid_H_mat

    fid_comoving_distance = jnp.trapezoid(fid_integrand, z_mat, axis = 1)
    gen_comoving_distance = jnp.trapezoid(gen_integrand, z_mat, axis = 1)

    cosmo_corr = (fid_comoving_distance/gen_comoving_distance)**2 * (H/H_fid)

    #analytical expression for the age of the universe in the flatLambdaCDM taken from astropy
    gen_prefactor = (2./3) * (s_Gyr/h) / jnp.sqrt(1-Om) 
    gen_arg = jnp.arcsinh( jnp.sqrt( (1 / Om - 1 + 0j) / (z + 1.0) ** 3 )) 
    t_gen = (gen_prefactor * gen_arg).real#Gyr

    fid_prefactor = (2./3) * (s_Gyr/70) / jnp.sqrt(1-0.3)
    fid_arg = jnp.arcsinh( jnp.sqrt( (1 / 0.3 - 1 + 0j) / (z + 1.0) ** 3 ))
    t_fid = (fid_prefactor * fid_arg).real #Gyr

    dL = (1+z) * gen_comoving_distance
    dV_dz_dtdet = ( (4*jnp.pi*c)/H ) * (dL**2 / (1 + z)**3)*1.e-9 #Gpc^3
    dz_dL = 1/( dL / (1 + z) + c * (1 + z) / H )

    #colossus_cosmo da sostituire con un calcolatore della funzione inversa di t(z), cioè z(t)
    return cosmo_corr, t_fid, t_gen, dL, dV_dz_dtdet, dz_dL

@jit
def db_Sch(log_M, log_M_tp, alph_1, alph_2, phi_1, phi_2):
        
    return jnp.log(10)*jnp.exp(-10**( log_M - log_M_tp ))*( phi_1*((10**( log_M - log_M_tp ))**(alph_1+1)) + phi_2*((10**( log_M - log_M_tp ))**(alph_2+1)) )

@jit
def SMF(mstarlog, z):  
    log_M_tp_vals = jnp.expand_dims(jnp.array([10.67, 10.83, 10.82, 10.87, 10.93, 10.86, 10.88, 10.75, 10.41, 10.32]), 1)
    alph_1_vals   = jnp.expand_dims(jnp.array([-1.41, -1.40, -1.35, -1.35, -1.49, -1.47, -1.47, -1.47, -1.47, -1.47]), 1)
    alph_2_vals   = jnp.expand_dims(jnp.array([0.14, -0.36, 0.06, 0.79, -0.56, -0.03, 0., 0., 0., 0.]), 1)
    phi_1_vals    = jnp.expand_dims(jnp.array([0.87, 0.70, 0.84, 0.74, 0.25, 0.24, 0.25, 0.22, 0.21, 0.13])*1.e-3, 1)
    phi_2_vals    = jnp.expand_dims(jnp.array([0.63, 0.32, 0.39, 0.13, 0.31, 0.13, 0., 0., 0., 0.])*1.e-3, 1)
        
    SMF_tmp = jnp.log10(db_Sch(mstarlog, log_M_tp_vals, alph_1_vals, alph_2_vals, phi_1_vals, phi_2_vals))

    z_edges_SMF   = jnp.array([0.2, 0.5, 0.8, 1.1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.5])
    z_centers_SMF = (z_edges_SMF[1:] + z_edges_SMF[:-1])/2. 

    def interpolate_fn(z_centers, SMF_tmp, j):
        return jnp.interp(z, z_centers, SMF_tmp[:, j], left="extrapolate", right='extrapolate')
        
    return vmap(interpolate_fn, in_axes = (None, None, 0))(z_centers_SMF, SMF_tmp, jnp.arange(len(mstarlog)))

@jit
def N_SB(mstarlog, z):

    zref=jnp.array([0.75,1.5,2.5,4.4])
    a1ref=jnp.array([-0.3,-0.3,-0.3,-0.3])
    b1ref=jnp.array([2.52,2.73,2.82,3.045])
    a2ref=jnp.array([-0.0067,-0.05,-0.075,-0.131])
    b2ref=jnp.array([0.1,0.555,0.817,1.408])
    m0logref=jnp.array([8.25,8.7,8.9,9.7])
    mfixlogref=jnp.array([7.6,7.99,8.24,9.25])
 
    a1= jnp.interp(z, zref, a1ref, left="extrapolate", right = "extrapolate")    
    b1= jnp.interp(z, zref, b1ref, left="extrapolate", right = "extrapolate")
    a2= jnp.interp(z, zref, a2ref, left="extrapolate", right = "extrapolate") 
    b2= jnp.interp(z, zref, b2ref, left="extrapolate", right = "extrapolate") 
    m0log= jnp.interp(z, zref, m0logref, left="extrapolate", right = "extrapolate") 
    mfixlog= jnp.interp(z, zref, mfixlogref, left="extrapolate", right = "extrapolate") 
    
    res=0.03*mstarlog/mstarlog

    res = jnp.where(mstarlog<=mfixlog, a1*mfixlog+b1, res)

    j1= jnp.logical_and(mstarlog>=mfixlog, mstarlog<=m0log)
    masked_mat = ((j1*mstarlog)*a1 + b1)*j1
    res_tmp = jnp.copy(res)
    res = jnp.where(~j1, jnp.nan, res) * 0     
    res = masked_mat + res
    res = jnp.where(jnp.isnan(res), 0, res)
    res = ~j1*res_tmp + res

    j2= jnp.logical_and(mstarlog>=m0log, mstarlog<=10.5)
    masked_mat = ((j2*mstarlog)*a2 + b2)*j2
    res_tmp = jnp.copy(res)
    res = jnp.where(~j2, jnp.nan, res) * 0 
    res = masked_mat + res
    res = jnp.where(jnp.isnan(res), 0, res)
    res = ~j2*res_tmp + res
    
    res = jnp.where(res<0.03, 0.03, res)  
    return res

@jit
def log_SFR_MS(log_M_st, cosmo_time):   
    a0_SFR = 0.20
    a1_SFR = -0.034
    b0_SFR = -26.134
    b1_SFR = 4.722
    b2_SFR = -0.1925
    return jnp.array((a1_SFR*cosmo_time+b1_SFR)*log_M_st + b2_SFR*log_M_st**2 + b0_SFR + a0_SFR*cosmo_time)

@jit
def resc_gauss(x, N, mu, sig):
    return N/jnp.sqrt(2*jnp.pi*sig**2)*jnp.exp(-((x-mu)**2/(2*sig**2)))
    
@jit
def SFR_dist(log_SFR, mstarlog, z, cosmo_time):
    sigma_MS = 0.2
    sigma_SB = 0.25
    offset_SB = 0.6
    
    frac_SB = N_SB(mstarlog, z)
    frac_MS = 1 - frac_SB
    log_psi_MS = log_SFR_MS(mstarlog, cosmo_time)
    log_psi_SB = log_psi_MS + offset_SB
    return resc_gauss(log_SFR, frac_MS, log_psi_MS, sigma_MS) + resc_gauss(log_SFR, frac_SB, log_psi_SB, sigma_SB)

@partial(jit, static_argnums=3)
def FMR(mstarlog, sfrlog, z , kind):

    offset_z1 = jnp.array([0.,0.,1.,2.,3.,4.255,5.775,7.665])
    offset_z2 = jnp.array([0.,1.,2.,3.,4.,4.98,6.86,8.47])
    offset_fmr1 = jnp.array([0.,0.,0.,0.,0.,-2.284e-1,-3.430e-1,-5.147e-1])
    offset_fmr2 = jnp.array([0.,0.,0.,0.,0.,-0.096,-0.086,-0.34])

    def Cu20_func(mstarlog, sfrlog, z):
        return 8.779-0.31/2.1*jnp.log10(1.+(10.**mstarlog/(10.**10.11*10.**(0.56*sfrlog)))**(-2.1))-8.76+jnp.log10(0.0153)
        
    def Cu23_func(mstarlog, sfrlog, z):
        return Cu20_func(mstarlog, sfrlog, z) + jnp.interp(z, offset_z1, offset_fmr1, left = "extrapolate", right = "extrapolate")

    def AM13_func(mstarlog, sfrlog, z):
        return 0.43*(mstarlog-0.66*sfrlog)+4.585-8.76+jnp.log10(0.0153) + jnp.interp(z, offset_z2, offset_fmr2, left = "extrapolate", right = "extrapolate")
        
    FMR_dict = {}

    FMR_dict["Cu23"] = Cu23_func
    FMR_dict["Cu20"] = Cu20_func
    FMR_dict["AM13"] = AM13_func
    
    return FMR_dict[kind](mstarlog, sfrlog, z)

@partial(jit, static_argnums=(4,5))
def Z_dist(log_Z, log_M_st, log_psi, z, kind, sigma_FMR):
    #sigma_FMR: 0.15 or 0.35
    #kind: "Cu23", "Cu20", "AM13"
    return resc_gauss(log_Z, 1., FMR(log_M_st, log_psi, z, kind), sigma_FMR) 

@partial(jit, static_argnums=(3,4))
def dSFRD_dlogZ(z, cosmo_correction, cosmo_time, kind, sigma_FMR):
    nmstar = 100
    mstarlog = jnp.linspace(6., 12., nmstar)

    nsfr = 100
    sfrlog = jnp.linspace(-4., 4., nsfr)
    
    def loop(sfrlog, mstarlog, zetalog, z, cosmo_correction, cosmo_time, k):
        
        q1_tmp = SFR_dist(sfrlog[:, jnp.newaxis], mstarlog[:, jnp.newaxis, jnp.newaxis], z[k], cosmo_time[k])*(10**sfrlog[:,jnp.newaxis]) 
        q2_tmp = Z_dist(zetalog, mstarlog[:, jnp.newaxis, jnp.newaxis], sfrlog[:, jnp.newaxis], z[k], kind, sigma_FMR)                      
        Q1_tmp = jnp.trapezoid(jnp.einsum('ijk, ijl ->lij', q1_tmp, q2_tmp), sfrlog)
        Q2_tmp = 10**jnp.array(SMF(mstarlog, z[k]))*cosmo_correction[k]
        
        return jnp.trapezoid(Q1_tmp * Q2_tmp, mstarlog, axis = 1)
        
    drhosfr_dlogZ = vmap(loop, in_axes = (None, None, None, None, None, None, 0))(sfrlog, mstarlog, zetalog, z, cosmo_correction, cosmo_time, jnp.arange(len(z)))
    
    return zetalog, drhosfr_dlogZ

@jit
def dSFRD_dlogM(z, cosmo_correction, cosmo_time):
    nmstar = 100
    mstarlog = jnp.linspace(6., 12., nmstar)

    nsfr = 100
    sfrlog = jnp.linspace(-4., 4., nsfr)

    def loop(sfrlog, mstarlog, z, cosmo_correction, cosmo_time, k):
        
        q1_tmp = SFR_dist(sfrlog, mstarlog[:, jnp.newaxis], z[k], cosmo_time[k])*(10**sfrlog)
        Q1_tmp = jnp.trapezoid(q1_tmp, sfrlog, axis = 1)
        Q2_tmp = 10**jnp.array(SMF(mstarlog, z[k]))*cosmo_correction[k]
        return Q1_tmp * Q2_tmp

    drhosfr_dlogM = vmap(loop, in_axes = (None, None, None, None, None, 0))(sfrlog, mstarlog, z, cosmo_correction, cosmo_time, jnp.arange(len(z)))
           
    return mstarlog, drhosfr_dlogM 

@jit
def full_SFRD(z, cosmo_correction, cosmo_time):
    
    mstarlog, drhosfr_dlogM = dSFRD_dlogM(z, cosmo_correction, cosmo_time)
    rhosfr = jnp.trapezoid(drhosfr_dlogM, mstarlog, axis = 1)
    return rhosfr

@partial(jit, static_argnums=(3,4))
def ddSFRD_dlogM_dlogZ(z, cosmo_correction, cosmo_time, kind, sigma_FMR):
    nmstar = 100
    mstarlog = jnp.linspace(6., 12., nmstar)

    nsfr = 100
    sfrlog = jnp.linspace(-4., 4., nsfr)

    def loop(sfrlog, mstarlog, zetalog , z, cosmo_correction, cosmo_time, k):
        
        q1_tmp = SFR_dist(sfrlog[:, jnp.newaxis], mstarlog[:, jnp.newaxis, jnp.newaxis], z[k], cosmo_time[k])*(10**sfrlog[:,jnp.newaxis]) 
        q2_tmp = Z_dist(zetalog, mstarlog[:, jnp.newaxis, jnp.newaxis], sfrlog[:, jnp.newaxis], z[k], kind, sigma_FMR)                      
        Q1_tmp = jnp.trapezoid(jnp.einsum('ijk, ijl ->lij', q1_tmp, q2_tmp), sfrlog)
        Q2_tmp = 10**jnp.array(SMF(mstarlog, z[k]))*cosmo_correction[k]
        return Q1_tmp * Q2_tmp
        
    ddrhosfr_dlogM_dlogZ = vmap(loop, in_axes = (None, None, None, None, None, None, 0))(sfrlog, mstarlog, zetalog, z, cosmo_correction, cosmo_time, jnp.arange(len(z)))

    return mstarlog, zetalog, ddrhosfr_dlogM_dlogZ
