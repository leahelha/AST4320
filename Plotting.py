import numpy as np
import pynbody
import matplotlib.pyplot as plt
from pynbody . analysis import hmf
import glob

'''DISCLAIMER: Although the filenames of the simulations are L20_64_dark.000xxx, the simulations ran for nSteps = 128, not 64'''


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
    plt.savefig(f'./Plots/Testing_HMF/{ahf_folders[i]}_halo_mass_function_none')

    
  
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
    plt.savefig(f'./Plots/Testing_DM4/{ahf_folders[i]}_10_most_massive_halos_nerror')
    #plt.show()

    '''Column density plot'''
    # visualize the dark matter in your snapshot
    pynbody.plot.image(snap.d, width = snap.properties ['boxsize'],
    av_z = True , units =('h**2 Msol kpc**-3 a**-3'))
    plt.title(f'z = {redshift[i]}', fontsize=fontsize)
    plt.savefig(f'./Plots/Column_density/{ahf_folders[i]}_Column denisty map centered')
    #plt.show()

    



"""
redshifted files I chose

L20_64_dark.000011, a = 0.16856184932537036, z = 4.932540512590884  
L20_64_dark.000014, a = 0.19782413880283825, z = 4.05499483557288   
L20_64_dark.000020, a = 0.25102850656084447, z = 2.9836113184923057 
L20_64_dark.000030, a = 0.33012887286629417, z = 2.029120086703265  
L20_64_dark.000054, a = 0.4974857324852427, z = 1.0101078979780866  
L20_64_dark.000080, a = 0.6678047077993748, z = 0.49744377109186977 

#Code for writing redshifts.txt and finding the redshifts
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