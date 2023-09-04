#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Derived from work by Harrison Nicholls on AEOLUS

import numpy as np
import shutil, os, sys, subprocess
from scipy.stats import binned_statistic

def LoadStellarSpectrum(spectrum):
    """Load a stellar spectrum from a file.

    Loads a stellar spectrum from a CSV file.

    Parameters
    ----------
        spectrum : str
            File containing stellar spectrum

    Returns
    ----------
        spec_wl : list
            Wavelength [nm]
        spec_fl : list
            Flux at 1 AU [erg s-1 cm-2 nm-1]
    """

    if os.path.isfile(spectrum):
        spec_data = np.loadtxt(spectrum, skiprows=2,delimiter='\t').T
        spec_wl = spec_data[0]
        spec_fl = spec_data[1]
    else:
        raise Exception("Cannot find stellar spectrum!")

    return spec_wl, spec_fl


def PrepareStellarSpectrum(wl, fl, star_file, nbins_max=95000):
    """Write a stellar spectrum.

    This function supplements InsertStellarSpectrum by writing a stellar 
    spectrum to a file in the format that SOCRATES expects. The flux needs to 
    be provided at 1 AU. Bins the spectrum down to at most nbins_max bins.

    Parameters
    ----------
        wl : list
            Wavelength [nm]
        fl : list
            Flux at 1 AU [erg s-1 cm-2 nm-1]
        star_file : str
            Path to output file, which will contain the stellar spectrum.
        nbins_max : int
            Number of spectral bins to down-sample to
            
    """

    socrates_nbins_max = int(1e5 - 1)
    nbins_max = min(nbins_max, len(wl)-1)

    # Validate
    if (len(wl) != len(fl)):
        raise Exception("Stellar wavelength and flux arrays have different lengths")
    if (len(wl) < 500):
        print("WARNING: Loaded stellar spectrum is very short!")
    if (nbins_max > socrates_nbins_max):
        raise Exception("Too many bins requested for stellar spectrum (maximum is %d)" % socrates_nbins_max)
    
    # Down-sample spectrum when necessary or requested
    if (len(wl) > socrates_nbins_max):
        print("Rebinning stellar spectrum")

        # Parameters
        max_retry           = 10    # Number of times to down-sample before giving up
        nbins_floor         = 250   # Minimum number of bins

        # Store old wl,fl arrays
        wl_orig = wl
        fl_orig = fl

        if (nbins_max < 500):
            print("WARNING: Requested number of bins is small (%d bins)" % nbins_max)

        nbins_max = min( int(socrates_nbins_max), nbins_max) # Must be fewer than 100k

        # Perform rounds of down-sampling until nbins < nbins_max
        # Rarely performs more than 2 rounds
        bins_retry = True
        downsample_factor = len(wl)/nbins_max
        count_retry = 0  
        while bins_retry:

            # Prevent loop from getting stuck
            if count_retry > max_retry:
                raise Exception("Giving up downsampling spectrum after %d rounds" % count_retry)
            
            # Ensure that the bins can be subdivided by truncating the array
            downsample_factor = int(downsample_factor)
            crop_amount = len(wl_orig) - int(len(wl_orig) % downsample_factor)
            wl = wl_orig[:crop_amount]

            # Downsample
            wl = wl.reshape(-1, downsample_factor).mean(axis=1)  # https://saturncloud.io/blog/the-optimal-approach-to-downsampling-a-numpy-array/
            nbins = len(wl)

            if (nbins < nbins_floor):
                raise Exception("Too few bins in stellar spectrum (%d bins)" % nbins)
            else:
                values, edges, _ = binned_statistic(wl_orig, fl_orig, statistic='mean', bins=wl)
                wl = 0.5 * (edges[:-1] + edges[1:])
                fl = values

            if np.isnan(fl).any() or np.any(fl <= 0):
                # Try again with fewer bins (occurs when there are empty bins)
                downsample_factor *= 1.5
            else:
                # Try again if still too many samples (unlikely case, but does happen sometimes)
                bins_retry = bool( nbins > nbins_max )
                downsample_factor *= 1.2

            count_retry += 1
            
    # Convert units
    wl = np.array(wl) * 1.0e-9  # [nm] -> [m]
    fl = np.array(fl) * 1.0e6   # [erg s-1 cm-2 nm-1] -> [W m-3]

    # Store header
    content = ""
    content += "Star spectrum at 1 AU. Created using PrepareStellarSpectrum() \n"
    content += "      WAVELENGTH        IRRADIANCE\n"
    content += "          (m)               (W/m3)\n"
    content += "*BEGIN_DATA\n"

    # Store body of data
    for i in range(len(wl)):
        content += str("      %1.7e      %1.7e\n" % (wl[i],fl[i]))

    # Store footer
    content += "*END\n"
    content += " "

    # Write content to file
    with open(star_file,'w') as handle:
        handle.write(content)


def InsertStellarSpectrum(orig_file:str, star_file:str, outp_file:str):
    """Insert a stellar spectrum.

    It's nice to be able to switch out the stellar spectrum for a different one. 
    This function takes in an original spectra file, with opacity data, and 
    inserts a stellar spectrum into a copy of it.

    Parameters
    ----------
        orig_file : str
            Path to original spectral file WITHOUT stellar spectrum.
        star_file : str
            Path to file containing stellar spectrum in the SOCRATES format.
        outp_file : str
            Path to output file, containing both opacity and stellar data.
            
    """

    # k files
    orig_filek = orig_file+"_k"
    outp_filek = outp_file+"_k"

    # Delete "new" files if they already exist
    if os.path.exists(outp_file):
        os.remove(outp_file)
    if os.path.exists(outp_filek):
        os.remove(outp_filek)

    # Copy original files to new location (retain old files)
    shutil.copyfile(orig_file,  outp_file)
    shutil.copyfile(orig_filek, outp_filek)

    # Run prep_spec from SOCRATES
    binary = "socrates/bin/prep_spec"
    inputs = [outp_file,'a','6','n','T','100 4000','100','2','n',star_file,'y','-1','EOF']
    p = subprocess.run([binary], stdout=subprocess.PIPE, input='\n'.join(inputs), encoding='ascii')
    if (p.returncode != 0):
        print("WARNING: prep_spec returned with code %d" % p.returncode)
    
# For standalone execution
# The arguments are:
# 1. stellar spectrum tsv file
# 2. socrates spectral file (source)
# 3. socrates spectral file (output)
if __name__=='__main__':

    print("Python: inserting stellar spectrum")

    args = sys.argv

    if len(args) != 4:
        raise Exception("Invalid number of arguments")
    
    star = str(args[1]) # path to stellar spectrum
    socr = str(args[2]) # path to socrates spectral file
    outp = str(args[3]) # path to output spectral file

    temp = ".star_temp"
    
    # Load data
    wl, fl = LoadStellarSpectrum(star)

    # Write data in SOCRATES format
    PrepareStellarSpectrum(wl, fl, temp)

    # Insert data into spectral file
    InsertStellarSpectrum(socr, temp, outp)

    # Remove temp file
    os.remove(temp)

