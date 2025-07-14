#!/usr/bin/env bash
# Download and unpack required and/or optional data
# All files can be found at https://zenodo.org/communities/proteus_framework

# Exit script if any of the commands fail
set -e

# Check that curl is installed
if ! [ -x "$(command -v curl)" ]; then
  echo "ERROR: curl is not installed" >&2
  echo "You must install curl in order to use this script" >&2
  exit 1
fi

# Check internet connectivity
header=$(curl -Is  https://zenodo.org | head -n 1)
if ! [[ $header == "HTTP/1.1 2"* || $header == "HTTP/1.1 3"* ]]; then
    # Return error if we don't get a positive HTTP response from Zenodo
    echo "ERROR: Failed to establish a connection to Zenodo"
    echo "Response: $header"
    exit 1
fi

# Root and resources folders
root=$(dirname $(realpath $0))
root=$(realpath "$root/..")
res="$root/res"
spfiles=$res/spectral_files
stellar=$res/stellar_spectra
surface=$res/surface_albedos
thermo=$res/thermodynamics
parfiles=$res/parfiles

# Make basic data folders
mkdir -p $res
mkdir -p $spfiles
mkdir -p $stellar
mkdir -p $surface
mkdir -p $thermo
mkdir -p $parfiles

# Help strings
help_basic="Get the basic data required to run the model"
help_highres="Get a spectral file with many high-resolution opacities"
help_steam="Get pure-steam spectral files"
help_anyspec="Get a particular spectral file by name, passing it as an argument"
help_stellar="Get a collection of stellar spectra"
help_surfaces="Get a collection of surface single-scattering albedos"
help_parfiles="Get a collection of gas linelist par files"
help_thermo="Get lookup data for thermodynamics (heat capacities, etc.)"
help="\
Download and unpack data used to run the model.

Call structure:
    $ get_data.sh [TARGET]

Where [TARGET] can be any of the following:
    basic
        $help_basic
    highres
        $help_highres
    steam
        $help_steam
    stellar
        $help_stellar
    anyspec
        $help_anyspec
    surfaces
        $help_surfaces
    parfiles
        $help_parfiles
    thermodynamics
        $help_thermo\
"

# Generic Zenodo downloader function
function zenodo {
    # $1 = Zenodo identifier for Record
    # $2 = target folder (on disk)
    # $3 = target filename (on disk and in Record)

    # target file path
    tgt="$2/$3"

    # exists?
    # if [[ -f "$tgt" ]]; then
    #     echo "    $1 > file already exists"
    #     return 0
    # fi

    # get data
    echo "    $1 > $tgt"
    mkdir -p $2
    curl -LsS "https://zenodo.org/records/$1/files/$3" > $tgt

    # check if command failed
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to download $1. Issue with curl command."
        exit 1
    fi

    # check file exists
    if [[ ! -f "$tgt" ]]; then
        echo "ERROR: Failed to download $1. File not found on disk."
        exit 1
    fi

    # check if file contains error message (replace NULL with blank)
    # header=$(head --bytes 100 $tgt)
    # if [[ $header == *"error"* || $header == *"Error"* ]]; then
    #     echo "ERROR: Failed to download from Zenodo Record $1"
    #     echo "Response: $header ..."
    #     exit 1
    # fi

    return 0
}

# Get zip file and extract it
function get_zip {
    # $1 = Zenodo record
    # $2 = target folder on disk
    # $3 = name of zip file in the Zenodo record

    zippath="$2/$3"

    zenodo $1 $2 $3
    unzip -oq $zippath -d $2
    rm $zippath
}

# Get a spectral file by name
function anyspec {
    # $1 = Codename (e.g. Honeyside)
    # $2 = Number of bands (e.g. 48)

    # This could be made much neater using associative arrays,
    #    but unfortunately they are not supported on MacOS

    # Key provided by user, used for match statement below
    name="$1$2"

    # Default filenames for .sf and .sf_k parts of spectral file
    sf_h="$1.sf"
    sf_k="$1.sf_k"

    # Get record on Zenodo containing the .sf and .sf_k files
    case $name in
        "Frostflow16" )
            rec="15799743"
            ;;
        "Frostflow48" )
            rec="15696415"
            ;;
        "Frostflow256")
            rec="15799754"
            ;;
        "Frostflow4096")
            rec="15799776"
            ;;

        "Dayspring16" )
            rec="15799318"
            ;;
        "Dayspring48" )
            rec="15721749"
            ;;
        "Dayspring256")
            rec="15799474"
            ;;
        "Dayspring4096")
            rec="15799495"
            ;;

        "Honeyside16" )
            rec="15799607"
            ;;
        "Honeyside48" )
            rec="15799652"
            ;;
        "Honeyside256")
            rec="15799731"
            ;;
        "Honeyside4096")
            rec="15696457"
            ;;

        "Oak318" )
            rec="15743843"
            ;;

        "Legacy318" )
            rec="15806343"
            sf_h="sp_b318_HITRAN_a16_no_spectrum"
            sf_k="sp_b318_HITRAN_a16_no_spectrum_k"
            ;;
        * )
            echo "ERROR: Unknown spectral file requested ($name)"
            exit 1
            ;;
    esac

    # Download the files
    zenodo $rec $spfiles/$1/$2 $sf_h
    zenodo $rec $spfiles/$1/$2 $sf_k
}

# Handle request for downloading a group of data
function handle_request {
    case $1 in
        "basic")
            echo $help_basic

            anyspec Oak 318
            anyspec Dayspring 48
            anyspec Honeyside 256

            zenodo 15721440 $stellar sun.txt

            get_zip 15805460 $thermo gases.zip
            ;;

        "highres")
            echo $help_highres
            anyspec Honeyside 4096
            ;;

        "steam")
            echo $help_steam
            anyspec Frostflow 16
            anyspec Frostflow 48
            anyspec Frostflow 256
            ;;

        "anyspec")
            echo $help_anyspec
            anyspec $2 $3
            ;;

        "stellar")
            echo $help_stellar
            rec="15721440"
            zenodo $rec $stellar trappist-1.txt
            zenodo $rec $stellar eps-eri.txt
            zenodo $rec $stellar hd97658.txt
            zenodo $rec $stellar gj1214.txt
            zenodo $rec $stellar sun.txt
            zenodo $rec $stellar l-98-59.txt
            ;;

        "surfaces")
            echo $help_surfaces
            rec="15805460"
            zenodo $rec $surface andesite.dat
            zenodo $rec $surface basalt_glass.dat
            zenodo $rec $surface basalt_tuff.dat
            zenodo $rec $surface diorite.dat
            zenodo $rec $surface gabbro.dat
            zenodo $rec $surface granite.dat
            zenodo $rec $surface harzburgite.dat
            zenodo $rec $surface hematite.dat
            zenodo $rec $surface lherzolite.dat
            zenodo $rec $surface lunar_anorthosite.dat
            zenodo $rec $surface lunar_marebasalt.dat
            zenodo $rec $surface mars_basalticshergottite.dat
            zenodo $rec $surface mars_breccia.dat
            zenodo $rec $surface norite.dat
            zenodo $rec $surface phonolite.dat
            zenodo $rec $surface pyrite.dat
            zenodo $rec $surface rhyolite.dat
            zenodo $rec $surface tephrite.dat
            zenodo $rec $surface tholeiitic_basalt.dat
            ;;

        "thermodynamics")
            echo $help_thermo
            get_zip 15805460 $thermo gases.zip
            ;;

        "parfiles")
            echo $help_parfiles
            rec="15806626"
            zenodo $rec $parfiles h2o-co2_4000-5000.par
            zenodo $rec $parfiles mixture_100-50000.par
            ;;

        *)
            echo "$help"
            return 0
            ;;
    esac
    return 0
}

handle_request $@

exit 0
