#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Derived from work by RJ Graham on AEOLUS

import numpy as np
import scipy.integrate as spint
import os, sys, time
from io import StringIO

# Function with necessary information for each species to calculate Rayleigh scattering
def species_info(species):
    info_dict = {}
    if species.lower() == 'co2':
        info_dict['Delta'] = 0.0805 #Vardavas and Carver 1984 for CO2
        info_dict['A'] = 43.9e-5 #Allen's Astrophysical Quantities (2002) Table 5.2
        info_dict['B'] = 6.4e-3 #Allen's
        info_dict['mu'] = 44e-3 #kg/mol
        
    elif species.lower() == 'n2':
        info_dict['Delta'] = 0.0305 #V&C 84 for N2
        info_dict['A'] = 29.06e-5 #Allen's Table 5.2
        info_dict['B'] = 7.7e-3 #Allen's  Table 5.2
        info_dict['mu'] = 28e-3 #kg/mol
    elif species.lower() == 'h2o':
        cross_section_hydrogen_1um=2.49e-6 #m^2/kg; mass-weighted cross section of hydrogen at 1 micron from Pierrehumbert 2010 Table 5.2
        info_dict['cross_section'] = 0.3743 * cross_section_hydrogen_1um #m^2/kg; mass weighted rayleigh cross section of water, from same table as h2
        info_dict['normalization_wavelength'] = 1e-6 #m; 1 micron
    else:
        raise Exception("Invalid species '%s'" % species)
    return info_dict

# Function for calculating raleigh scattering coefficient at a given wavelength 
# for a species w/ given index of refraction info (info which is obtained from species_info function)
def cross_section(wavelength, info_dict):#wavelength in m 
    if 'Delta' in info_dict:
        Delta = info_dict['Delta']
        A = info_dict['A']
        B = info_dict['B']
        mu = info_dict['mu']
        microns = wavelength * 10**6
        delta = (6+3*Delta)/(6-7*Delta) #V&C 1984
        coefficient = 4.577e-21 * delta/microns ** 4 * (A * (1 + B / microns**2)) ** 2 * (10**-4 * 6.022e23) / mu #V&C 1984; m^2/kg
    else:
        cross_section_1um = info_dict['cross_section']
        normalization_wavelength = info_dict['normalization_wavelength']
        coefficient = cross_section_1um * normalization_wavelength ** 4. / wavelength ** 4
    return coefficient


# This function integrates the Rayleigh cross-section for each species over a given band to allow for 
# proper averaging and then adds the integrated values together, weighted by molar mixing ratio
def band_integrator(species_list, molar_mixing_ratio_list, wavelength1_list, wavelength2_list):
    total_coefficient = 0
    for n, species in enumerate(species_list):
        info_tup = species_info(species)
        
        #Coefficient w/o mixing ratio weight
        partial_coefficient = np.vectorize(spint.quad,excluded=['args'])(cross_section,wavelength1_list,wavelength2_list, args = info_tup)[0]
        
        #Coefficient weighted  by mixing ratio
        total_coefficient +=  molar_mixing_ratio_list[n] * partial_coefficient
        
    return total_coefficient

# This is the function that adds  Rayleigh coefficients to spectral files
def rayleigh_coeff_adder(species_list = ['co2'], mixing_ratio_list = [1.], spectral_file_path= './spectral_files/sp_b318_HITRAN_a16_RS/sp_b318_HITRAN_a16'):
    spectral_file = open(spectral_file_path,'r')
    
    wavelength_band_string = ""
    
    #Make a loop that runs through the lines in the spectral file until it finds BLOCK 1,
    #at which point it begins to "pay attention" and read the relevant lines w/ the wavelength data in that file
    #into a new file that stores the wavelength bands'''
    pay_attention = False
    
    beginning_block_text_list = []
    
    for line in spectral_file:
        
        if pay_attention:
            #If the line says END and pay_attention is True, then we've reached the
            #end of block 1, so the loop can terminate because we don't need to read
            # the rest of the spectral file
            if 'END' in line:
                break
            # The first couple of lines in the block  have unnecessary text;
            # We ignore those by specifying that the loop only writes a line out
            # that starts with a space, which makes the loop ignore the unnecessary
            # bits and only grab the wavelength bands
            if line[0]==' ':
                wavelength_band_string += line
                if line[-1] != "\n":
                    wavelength_band_string += "\n"
            else:
                beginning_block_text_list.append(line)
        if '*BLOCK: TYPE =    1' in line:
            #The lines following this line will describe the wavelength data
            pay_attention = True
    
    spectral_file.close()

    #Now we generate an array from the wavelength band file
    wavelength_bands = np.genfromtxt(StringIO(wavelength_band_string),usecols=np.arange(0,3))

    # This function calculates the rayleigh scattering cross section for a species at a given wavelength
    
    #calculate the change in wavelength over each band
    dlambda = wavelength_bands[:,2]-wavelength_bands[:,1]
    
    #Integrate the cross section function over the band
    cross_section_integral = band_integrator(species_list,mixing_ratio_list,wavelength_bands[:,1],wavelength_bands[:,2])
    
    #Divide each integral by dlambda to get the average cross section over each band
    cross_section_avgs = cross_section_integral/dlambda
    
    #Put the numerical labels and the cross section averages together in one array
    cross_section_list = np.vstack((wavelength_bands[:,0],cross_section_avgs)).T
    
    #A list of the first couple of lines in BLOCK 3, which contains the Rayleigh scattering coefficients
    beginning_block_3 = ['*BLOCK: TYPE =    3: SUBTYPE =    0: VERSION =    0\n',\
                         'Rayleigh mass scatering coefficients at STP: unit m**2/kg\n',\
                         'Band        Rayleigh coefficient\n',\
                         '                 m2/kg\n']
    # This list will hold all the lines in Block 3
    block3_list = []
    # This loop adds the first few lines to block 3's list
    for n in range(len(beginning_block_3)):
        block3_list.append(beginning_block_3[n])
    
    # This loop puts the cross sections in the proper format and then adds them to the list line by line
    for n in range(cross_section_list.shape[0]):
        #print('{0:5d}'.format(np.array(cross_section_list[n,0],dtype=int))+'\t\t %.9E\n'%(cross_section_list[n,1]))
        block3_list.append('{0:5d}'.format(np.array(cross_section_list[n,0],dtype=int))+'        %.9E\n'%(cross_section_list[n,1]))
    
    # Adding the end statement of block 3
    block3_list.append('*END\n')
    # Opening the spectral file and reading all of its lines out
    spectral_file = open(spectral_file_path,'r')
    spectral_file_lines = spectral_file.readlines()
    spectral_file.close()
    
    # We're going to loop through the spectral file's lines and put them into lists
    # depending on whether they occur before or after block 3 -- 
    block3_flag = False
    block4_and_up_flag = False
    
    # Wwe make two more lists, the first will contain each line preceding block 3 and 
    # the second will contain each line after block 3
    first_blocks_list = []
    last_blocks_list = []
    for line in spectral_file_lines:
        
        #If we encounter this string in a line, we've entered block 3
        if '*BLOCK: TYPE =    3' in line:
            block3_flag = True
            
        # If we encounter this string in a line, we've entered block 4
        elif '*BLOCK: TYPE =    4' in line:
            block4_and_up_flag = True
            
            # block3_Flag is false now because we've left block 3
            block3_flag = False
            
        # Before we reach block 3 or block 4, lines go into first_blocks_list
        if not block3_flag and not block4_and_up_flag:
            first_blocks_list.append(line)
            
        # After we've made it through block 3 and entered block 4, lines go into
        # last_blocks_list
        elif block4_and_up_flag:
            last_blocks_list.append(line)
    
    #now we create a dummy file to read all of the blocks into
    temp_spectral_file_path =  spectral_file_path+'_temp_spectral_file'
    temp_spectral_file = open(temp_spectral_file_path,'w')
    
    #First, blocks 0 to 2 are read in
    for n in range(len(first_blocks_list)):
        temp_spectral_file.write(first_blocks_list[n])
    
    #Then, block 3 is read in
    for n in range(len(block3_list)):
        temp_spectral_file.write(block3_list[n])
    
    # Finally, the remaining blocks are read in  
    for n in range(len(last_blocks_list)):
        temp_spectral_file.write(last_blocks_list[n])
        
    temp_spectral_file.close()
    time.sleep(20/1000.0) # for filesystem to catch up, wait 20 ms
    
    # The last step is to delete the original spectral file and replace it with the new one
    os.remove(spectral_file_path)
    os.rename(temp_spectral_file_path,spectral_file_path)

# For executing as standalone script.
# The arguments are:
# 1. spectral file
# 2. CO2 mixing raito
# 3. N2 mixing ratio
# 4. H2O mixing ratio
if __name__ == "__main__":

    print("Python: Inserting Rayleigh scattering")

    args = sys.argv

    if len(args) > 5:
        raise Exception("Invalid number of cmd line arguments")

    # Arg 1 
    spectral_file_path = str(args[1])

    # Args 2,3,4
    species_list = []
    mixing_ratio_list = []

    if float(args[2]) > 0:
        species_list.append("co2")
        mixing_ratio_list.append(float(args[2]))

    if float(args[3]) > 0:
        species_list.append("n2")
        mixing_ratio_list.append(float(args[3]))

    if float(args[4]) > 0:
        species_list.append("h2o")
        mixing_ratio_list.append(float(args[4]))

    # Normalise mixing ratios
    mr_tot = np.sum(mixing_ratio_list)
    mixing_ratio_list = list( np.array(mixing_ratio_list)/mr_tot )

    rayleigh_coeff_adder(species_list=species_list, 
                         mixing_ratio_list=mixing_ratio_list, 
                         spectral_file_path=spectral_file_path
                         )
