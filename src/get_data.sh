#!/usr/bin/env bash
# Download and unpack required and/or optional data
# All files can be found at https://zenodo.org/communities/proteus_framework

# Exit script if any of the commands fail
# set -e

# Check that wget is installed
if ! [ -x "$(command -v wget)" ]; then
  echo "ERROR: wget is not installed" >&2
  echo "You must install wget in order to use this script" >&2
  exit 1
fi

# Check that unzip is installed
if ! [ -x "$(command -v unzip)" ]; then
  echo "ERROR: unzip is not installed" >&2
  echo "You must install unzip in order to use this script" >&2
  exit 1
fi

# Check internet connectivity
header=$(wget --spider -S "https://zenodo.org" 2>&1 | grep "HTTP/")
# echo $header
if ! [[ $header == *"HTTP/1.1 2"* || $header == *"HTTP/1.1 3"* ]]; then
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
help_dryrun="Test the get_data script"
help_basic="Get the basic data required to run the model"
help_highres="Get a spectral file with many high-resolution opacities"
help_steam="Get pure-steam spectral files"
help_anyspec="Get a particular spectral file by name, passing it as an argument"
help_stellar="Get a collection of stellar spectra"
help_surf_standard="Get a basic collection of surface reflectance data"
help_surf_extended="Get an extended collection of surface reflectance data"
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
        $help_surf_standard
    surfaces_extended
        $help_surf_extended
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

    # target url
    url="https://zenodo.org/records/$1/files/$3"

    # get data
    echo "    zenodo/$1 > $tgt"
    mkdir -p $2
    wget -qO $tgt $url

    # check if command failed or if file does not exist
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to download $1. Issue with wget command"
    elif [[ ! -f "$tgt" ]]; then
        echo "ERROR: Failed to download $1. File not found on disk."
    else
        return 0
    fi

    # try again at downloading the file?
    echo "Trying again to download the file"
    sleep 1
    wget -qO $tgt $url

    # check if command failed or if file does not exist
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to download $1. Issue with wget command"
        exit 1
    elif [[ ! -f "$tgt" ]]; then
        echo "ERROR: Failed to download $1. File not found on disk."
        exit 1
    fi

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

        "dryrun")
            echo $help_dryrun
            echo "Sleeping for 3 seconds..."
            sleep 3
            ;;

        "basic")
            echo $help_basic

            anyspec Oak 318
            # anyspec Dayspring 16
            anyspec Dayspring 48
            # anyspec Honeyside 256

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
            anyspec Frostflow 4096
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
            echo $help_surf_standard
            rec="15880455"
            zenodo $rec $surface andesite.dat
            zenodo $rec $surface basaltglass.dat
            zenodo $rec $surface basalttuff.dat
            zenodo $rec $surface diorite.dat
            zenodo $rec $surface gabbro.dat
            zenodo $rec $surface granite.dat
            zenodo $rec $surface harzburgite.dat
            zenodo $rec $surface hematite.dat
            zenodo $rec $surface lherzolite.dat
            zenodo $rec $surface lunaranorthosite.dat
            zenodo $rec $surface lunarmarebasalt.dat
            zenodo $rec $surface marsbasalticshergottite.dat
            zenodo $rec $surface marsbreccia.dat
            zenodo $rec $surface norite.dat
            zenodo $rec $surface phonolite.dat
            zenodo $rec $surface pyrite.dat
            zenodo $rec $surface rhyolite.dat
            zenodo $rec $surface tephrite.dat
            zenodo $rec $surface tholeiiticbasalt.dat
            ;;

        "surfaces_extended")
            echo $help_surf_extended
            get_zip 15881238 $surface ecostress.zip
            get_zip 15881496 $surface lavaworld.zip
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
            ;;
    esac
    return 0
}

handle_request $@

exit 0
