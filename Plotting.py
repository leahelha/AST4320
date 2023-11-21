import numpy as np
import pynbody
import matplotlib.pyplot as plt
from pynbody . analysis import hmf
import glob

'''DISCLAIMER: Although the filenames of the simulations are L20_64_dark.000xxx, the simulations ran for nSteps = 128, not 64'''
"""
L20_64_dark.000011, a = 0.16856184932537036, z = 4.932540512590884  X
L20_64_dark.000014, a = 0.19782413880283825, z = 4.05499483557288   X
L20_64_dark.000020, a = 0.25102850656084447, z = 2.9836113184923057 X
L20_64_dark.000030, a = 0.33012887286629417, z = 2.029120086703265  X
L20_64_dark.000054, a = 0.4974857324852427, z = 1.0101078979780866  X
L20_64_dark.000080, a = 0.6678047077993748, z = 0.49744377109186977 X
L20_64_dark.000128, z = 0
"""

ahf_folders = ['AHF_z0', 'AHF_z05', 'AHF_z1', 'AHF_z2', 'AHF_z3', 'AHF_z4', 'AHF_z5'] 
files = ['L20_64_dark.000128', 'L20_64_dark.000080', 'L20_64_dark.000054', 'L20_64_dark.000030', 'L20_64_dark.000020', 'L20_64_dark.000014', 'L20_64_dark.000011']

redshift = ['0', '0.5', '1', '2', '3', '4', '5']
col = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#17becf', '#e377c2', '#7f7f7f', '#bcbd22']

halo_masses = []

fontsize=12
for i in range(len(ahf_folders)):
    snap = pynbody.load(f"{ahf_folders[i]}/{files[i]}")
    #print('a =', snap.properties['a'])

    halos = snap.halos()
    # put snapshot into physical units ,
    # and define which units to use
    snap.physical_units(distance='kpc a h^-1', mass='Msol h^-1')
    
    'HMF'
    dmbin = 0.1 #0.05 #0.1

    # theoretical Press - Schechter ( PS ) halo mass function
    M_theory, _, HMF_theory = hmf.halo_mass_function(snap, kern="PS", log_M_min=10, log_M_max=15, delta_log_M=dmbin)
    # halo mass function from data
    #M_data, HMF_data, error = hmf.simulation_halo_mass_function(snap)

    M_data, HMF_data, error = hmf.simulation_halo_mass_function(snap, log_M_min=10, log_M_max=15, delta_log_M=dmbin)
    #print(M_data)
    
    
  
    
    plt.figure()
    plt.loglog(M_theory, HMF_theory, color='k')
    #plt.errorbar(M_data, HMF_data, yerr=error, fmt='o', capthick='0.1', elinewidth='1', color='red', markersize=3)
    plt.errorbar(M_data, HMF_data, yerr=error, fmt='bo', markersize=3, label='Our results')
    

    if i == 90:
        '''Going through other student data + Sijing's data'''
        s = ['15502', '15506', '15507', '15509', '15510', 'Shen']
        
        for k in s:
            file = glob.glob(f'Comparisons/{k}/*000128') #Lists any file ending in .000128
            #print(f"{file}")
            print(k)
            snap2 = pynbody.load(file[0])
            # put snapshot into physical units ,
            # and define which units to use
            halos = snap2.halos()
            mass = snap2.properties
            print(mass)
            snap2.physical_units(distance='kpc a h^-1', mass='Msol h^-1')
            
            
            '''HMF'''
            #dmbin = 0.1

            # theoretical Press - Schechter ( PS ) halo mass function
            M_theory, _, HMF_theory = hmf.halo_mass_function(snap2, kern="PS", log_M_min=10, log_M_max=15, delta_log_M=dmbin)
            # halo mass function from data
            #M_data, HMF_data, error = hmf.simulation_halo_mass_function(snap)

            M_data, HMF_data, error = hmf.simulation_halo_mass_function(snap2, log_M_min=10, log_M_max=15, delta_log_M=dmbin)
            
            #print(M_data)
            
            #help(hmf.halo_mass_function)
            
            #plt.loglog(M_theory, HMF_theory)
            #plt.errorbar(M_data, HMF_data, yerr=error, fmt='o', capthick='0.1', elinewidth='1', color='red', markersize=3)
            plt.errorbar(M_data, HMF_data, yerr=error, fmt='o', markersize=3, label=f'{k}') 
   
    plt.legend()
    plt.xlabel('$log_{10}(M)$  $[M_{sol} h^{-1}]$', fontsize=fontsize)
    plt.ylabel('$dN/(dlog_{10}M)$  $[Mpc^{-3} h^3]$', fontsize=fontsize)
    plt.title(f'z = {redshift[i]}', fontsize=fontsize)
    plt.savefig(f'./plots/Testing_HMF/{ahf_folders[i]}_halo_mass_function_none')

    exit()
  
    #plt.show()
    

    'DM density profiles'

    from pynbody . analysis . theoretical_profiles import NFWprofile
        # load in halos
    halos = snap.halos()
    #print(f'Nr of halos {len(halos)}')
    #halos are numbered from 1 , not 0
    
    'Plotting 10 most massive halos'
    plt.figure()
    for j in range(10):
        print(halos[j+1].properties['Mhalo'],f'\n')
        mass = halos[j+1].properties['Mhalo']
        halo_masses.append(mass)

        halo = halos[j+1]
        halo_properties = halo.properties

        shift_halo = pynbody.analysis.angmom.halo.center(halos[j+1], 'pot', vel=False, wrap=True) #centers whole simulation to the halo center

        # calculate a theoretical profile from halo properties
        NFW_theory = NFWprofile(halo_radius = halo_properties["Rhalo"], concentration = halo_properties["cNFW"], halo_mass = halo_properties["Mhalo"])

        theoretical_profile = NFW_theory.profile_functional
        
        profiles = pynbody.analysis.profile.Profile(snap.d, rmin='.01 kpc', rmax=str(halos[j+1].properties['Rhalo']) + ' kpc', ndim=3, nbins=20) #*** Changed nbins from 100


        # results, can be used instead of r-theory 
        #color = col[i]
        r = profiles['rbins'].in_units('kpc a h**-1')   # in_units('kpc') #tried: in_units('kpc a h**-1')  
        #print(profiles._properties)
        rho = profiles['density'].in_units('h**2 Msol a**-3 kpc**-3')#.in_units('Msol h**-1') #in_units('Msol kpc**-3')
        #print(profiles['density'].units)
        #print(r)
        
        rho_theory = theoretical_profile(r)

        
        #_error = profiles['density']*profiles['n']**(-1/2)

        color = col[j]
        plt.semilogy(r, rho, linestyle=':', color=color, alpha = 0.4)
        #plt.errorbar(r, rho, yerr=_error, fmt='o', markersize=1, color=color)
        
        if j == 0:
            plt.semilogy(r, rho_theory, linestyle='-', color='black', label = f'Most massive theoretical halo')
            plt.semilogy(r, rho, linestyle=':', color='black', label = f'Most massive halo')

    plt.legend()
    plt.title(f'z = {redshift[i]}',fontsize=fontsize)
    plt.xscale('log')
    plt.xlabel('Radius [kpc a $h^{-1}$]',fontsize=fontsize)
    plt.ylabel('Density [$h^2 M_{sol} a^{-3} kpc^{-3}$]',fontsize=fontsize)
    plt.savefig(f'./plots/Testing_DM4/{ahf_folders[i]}_10_most_massive_halos_nerror')
    #plt.show()

    '''Column density plot'''
    # visualize the dark matter in your snapshot
    pynbody.plot.image(snap.d, width = snap.properties ['boxsize'],
    av_z = True , units =('h**2 Msol kpc**-3 a**-3'))
    plt.title(f'z = {redshift[i]}', fontsize=fontsize)
    plt.savefig(f'./plots/Column_density/{ahf_folders[i]}_Column denisty map centered')
    #plt.show()

    
    



exit()
'Plotting everything together'
plt.figure()
for i in range(len(ahf_folders)):
    snap = pynbody.load(f"{ahf_folders[i]}/{files[i]}")
    print('a =', snap.properties['a'])

    # put snapshot into physical units ,
    # and define which units to use
    snap.physical_units(distance='kpc a h^-1', mass='Msol h^-1')
    
    # load in halos
    halos = snap.halos()
    #halos are numbered from 1 , not 0
    #print(halos[1].properties)
    dmbin = 0.1

    # theoretical Press - Schechter ( PS ) halo mass function
    M_theory, _, HMF_theory = hmf.halo_mass_function(snap, kern="PS", log_M_min=10, log_M_max=15, delta_log_M=dmbin)
    # halo mass function from data
    #M_data, HMF_data, error = hmf.simulation_halo_mass_function(snap)

    M_data, HMF_data, error = hmf.simulation_halo_mass_function(snap, log_M_min=10, log_M_max=15, delta_log_M=dmbin)
    #print(M_data)


    #help(hmf.halo_mass_function)
    
#     plt.loglog(M_theory, HMF_theory, color='k')
#     #plt.errorbar(M_data, HMF_data, yerr=error, fmt='o', capthick='0.1', elinewidth='1', color='red', markersize=3)
#     plt.errorbar(M_data, HMF_data, yerr=error, fmt='o', markersize=3, label=f'z{i}')

# plt.xlabel('log(M) [Msol]')
# plt.ylabel('$dN/(d\log_{10}M)$  [Mpc$^{-3} h^3$]')
# plt.legend()

# plt.savefig(f'./plots/All_halo_mass_function')



'Attempting to plot top 10 massive dark matter density profiles'

"""

_r = np.zeros((len(ahf_folders), 10))
_rho = np.zeros((len(ahf_folders), 10))
_rho_theory = np.zeros((len(ahf_folders), 10))
_error = np.zeros((len(ahf_folders), 10))
for j in range(10):
    for i in range(len(ahf_folders)):
        snap = pynbody.load(f"{ahf_folders[i]}/{files[i]}")
        #print('a =', snap.properties['a'])

        # put snapshot into physical units ,
        # and define which units to use
        snap.physical_units(distance='kpc a h^-1', mass='Msol h^-1')
        
        # load in halos
        halos = snap.halos()
        #halos are numbered from 1 , not 0
        #print(halos[1].properties)

        from pynbody . analysis . theoretical_profiles import NFWprofile

        halo = halos[1]
        halo_properties = halo.properties

        shift_halo = pynbody.analysis.angmom.halo.center(halos[1], 'pot', vel=False, wrap=True) #centers whole simulation to the halo center

        # calculate a theoretical profile from halo properties
        NFW_theory = NFWprofile(halo_radius = halo_properties["Rhalo"], concentration = halo_properties["cNFW"], halo_mass = halo_properties["Mhalo"])

        theoretical_profile = NFW_theory.profile_functional

        profiles = pynbody.analysis.profile.Profile(snap.d, rmin='.01 kpc', rmax=str(halos[1].properties['Rhalo']) + ' kpc', ndim=3, nbins=100)

        # results, can be used instead of r-theory 

        r = profiles['rbins'].in_units('kpc')
        rho = profiles['density'].in_units('Msol kpc**-3')
        rho_theory = theoretical_profile(r)
        error = profiles ['density']*profiles['n']**(-1/2)

        # _r.append(r)
        # _rho.append(rho)
        # _rho_theory.append(rho_theory)
        # _error.append(error)
        _r[i, j] = r[j]
        _rho[i, j] = rho[j]
        _rho_theory[i, j] = rho_theory[j]
        _error[i, j] = error[j]


        #plt.semilogy(r, rho)

plt.errorbar(_r[0,:], _rho[0, :], yerr=_error[0, :], fmt='o', markersize=2)
plt.semilogy(_r[0, :], _rho_theory[0, :])
    
plt.legend()
plt.xlabel('kpc')
plt.ylabel('Msol kpc^-3')
plt.show()

print(np.shape(_r))
# plt.savefig(f'./plots/All_NFW_profile_halomass_top_10')
   """ 

"""
Note: the fitting step will fail if there are zeros in the density array rho.
To avoid this, you should add a step that checks removes places where ρ = 0.
Alternatively, you can try to make the bins larger when you calculate the density.
"""



"""
Group session questions:
- Column density map for L20_64_dark.000128 was supposed to be at a redshift 0, snap.properties['a'] returns value of a=0.9999999999999978

- Error messages:
    - UserWarning: Halo finder masses not provided. Calculating them (might take a while...)
    - warnings.warn("Your bin range does not encompass the full range of halo masses"
    - RuntimeWarning: divide by zero encountered in power
    - error = profiles["density"]* profiles["n"]**(-1/2)
    - 
   
    
"""

"""
What i got for running ahf-v1.0-110/bin/AHF-v1.0-110 AHF.000128.input:

These are the input parameters you provided, please check them carefully again:
===============================================================================
ic_filename       = /uio/hume/student-u49/leaheh/AST4320_project2/L20_64_dark.000128
ic_filetype       = 90
outfile_prefix    = L20_64_dark.000128
time    = 1.000000
nbodies = 262144
ndim    = 3
nsph    = 0
ndark   = 262144
nstar   = 0
"""

"""
Notes for report:
- Column density plot. Discuss qualitatively what you observe in these snapshots in terms of number of halos, their sizes, the cosmic filaments etc.
- *** check units on NWF plot,  'kpc a h**-1' is in the properties

"""

"""
redshifted files I chose

L20_64_dark.000011, a = 0.16856184932537036, z = 4.932540512590884  
L20_64_dark.000014, a = 0.19782413880283825, z = 4.05499483557288   
L20_64_dark.000020, a = 0.25102850656084447, z = 2.9836113184923057 
L20_64_dark.000030, a = 0.33012887286629417, z = 2.029120086703265  
L20_64_dark.000054, a = 0.4974857324852427, z = 1.0101078979780866  
L20_64_dark.000080, a = 0.6678047077993748, z = 0.49744377109186977 

#Code for finding this
files = []
for i in range(9):
    s = f'00{i+1}'
    files.append(s)

for i in range(10,100):
    s = f'0{i}'
    files.append(s)

for i in range(100,128):
    s = f'{i}'
    files.append(s)

z_ = [5, 4, 3, 2, 1, 0.5]
with open('redshifts.txt', 'w') as outfile:
    for i in files:
        snap = pynbody.load(f"AST4320_project2/L20_64_dark.000{i}")
        a = snap.properties['a']
        z = 1/a - 1

        outfile.write(f'L20_64_dark.000{i}, a = {str(a)}, z = {str(z)}\n')


"""

#/mn/stornext/d17/extragalactic/personal/shens/ahf-v1.0-110/

### comamands for AMH finder

# emacs -nw  AHF.000128.input #open and edit can replace emacs with vim
#    # run module intel first
# ahf-v1.0-110/bin/AHF-v1.0-110 AHF.000128.input # run AHF, change the 000128 bit after editing file

#/mn/stornext/d17/extragalactic/personal/shens/L40_128_dark