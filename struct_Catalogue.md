**Structure of file Catalogue.txt**

```#ID``` : unique identification number of the merger

```Metal```: metallicity

```Nrec```: number of previous mergers

```EnvType```: type of environment ("none" for isolated binaries, "YC", "GC", "NC" for young, globular, and nuclear clusters, respectively)

```lab```: unused label

```m1[Msun]```: mass of the primary in Msun

```m2[Msun]```: mass of the secondary in Msun

```a1```: dimensionless spin of the primary

```a2```: dimensionless spin of the secondary

```Mfin[Msun]``` afin xeff vGW[km/s]: remnant mass (in Msun), dimensionless spin, effective spin parameter, and GW recoil (in km/s), respectively

```tfor[yr]```: formation time of the environment in which BH progenitors formed in yr

```tlast_mer[yr]```: delay time (time elapsed from formation to merger. Note: for repeated mergers, this time includes all previous mergers' timeline)

```Mclu_t0[Ms] Rclu_t0[pc] Vesc[km/s]```: initial mass (in msun), half-mass radius (in pc), and escape velocity (in km/s) of the cluster hosting the merger

```BinaryStatus```: status of the binary ("ejected", "GWrecoiled", "in-cluster", "highlander")

```aeje[AU] aGW[AU]```: semi-major axis (in AU) below which the binary is ejected or GW-driven evolution dominates, respectively

```nBHs```: number of BHs initially available in the host cluster

```Mcore_th[Ms] rcore_th[pc]```: initial mass (in Msun) and radius (in pc) of the cluster core

```redshift_merger redshift_formation```: redshift at formation and last merger

```tSMBH[yr]``` redshiftSMBH: time at which the BH mass exceeds a threshold, and corresponding redshift

```mprog[Ms]```: mass of the progenitor BH (in case of repeated mergers this differs from the primary)

```eccentricity semimajoraxis[AU] acrit[AU]```: eccentricity and semimajor axis (in AU) of the merging BBH

```tmerger[yr]```: time elapsed from decoupling of the BBH from the dynamics of the cluster and the merger (in yr)

```cos(angle_s1s2) cos(angle_s1L) cos(angle_s2L)```: angles between the two spin vectors, the spin of the primary and the orbital angular momentum, the spin of the secondary and the orbital angular
momentum, respectively
