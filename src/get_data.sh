#!/usr/bin/env bash
# Download and unpack required and/or optional data
# All files can be found on OSF at https://osf.io/8dumn/

# Exit script if any of the commands fail
set -e

# Check that curl is installed
if ! [ -x "$(command -v curl)" ]; then
  echo 'ERROR: curl is not installed.' >&2
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

# Generic OSF downloader function
function osf {
    # $1 = OSF identifier
    # $2 = target folder
    # $3 = target filename

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
    curl -LsS "https://osf.io/download/$1/" > $tgt

    # check file exists
    if [[ ! -f "$tgt" ]]; then
        echo "ERROR: Failed to download $1"
        exit 1
    fi

    # check if file contains error message
    header=$(head --bytes 80 $tgt)
    if [[ $header == *"message_short"* ]]; then
        echo "ERROR: Failed to download $1"
        echo "Response: $header ..."
        exit 1
    fi

    return 0
}

# Get zip file
function get_zip {
    # $1 = OSF identifier
    # $2 = target folder

    zipfile=".temp.zip"
    osf $1 $2 $zipfile
    unzip -oq $2/$zipfile -d $2
    rm $2/$zipfile
}

# Get a spectral file by name (associative arrays not supported on MacOS)
function anyspec {
    # $1 = Codename (e.g. Honeyside)
    # $2 = Number of bands (e.g. 48)

    name="$1$2"

    # Get identifiers
    case $name in
        "Frostflow16" )
            id_a="6rvfe"; id_b="kxve8"
            ;;
        "Frostflow48" )
            id_a="9n6mw"; id_b="xfap8"
            ;;
        "Frostflow256")
            id_a="mnvyq"; id_b="tzsgc"
            ;;
        "Frostflow4096")
            id_a="eyw6b"; id_b="ry6qz"
            ;;

        "Dayspring16" )
            id_a="uwfja"; id_b="j7f2w"
            ;;
        "Dayspring48" )
            id_a="heuza"; id_b="c5jv3"
            ;;
        "Dayspring256")
            id_a="b5gsh"; id_b="dn6wh"
            ;;
        "Dayspring4096")
            id_a="g5nh8"; id_b="htn3q"
            ;;

        "Honeyside16" )
            id_a="6cqnp"; id_b="f5snk"
            ;;
        "Honeyside48" )
            id_a="rxj5a"; id_b="xfc5h"
            ;;
        "Honeyside256")
            id_a="97436"; id_b="xny6v"
            ;;
        "Honeyside4096")
            id_a="p672d"; id_b="ujb4z"
            ;;

        "Oak318" )
            id_a="qmp4e"; id_b="5fxr7"
            ;;

        "Legacy318" )
            id_a="b7dvn"; id_b="m8zf5"
            ;;
        * )
            echo "ERROR: Unknown spectral file requested"
            exit 1
            ;;
    esac

    # Download the file
    osf $id_a $spfiles/$1/$2 $1.sf
    osf $id_b $spfiles/$1/$2 $1.sf_k
}

# Handle request for downloading a group of data
function handle_request {
    case $1 in
        "basic")
            echo $help_basic

            anyspec Oak 318
            anyspec Dayspring 48
            anyspec Honeyside 256
            osf 6k8ba $spfiles reference.pdf

            osf 2qdu8 $stellar sun.txt

            get_zip 4m5x8 $thermo
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

            osf mabp2 $stellar trappist-1.txt
            osf rk7mj $stellar eps-eri.txt
            osf agsrq $stellar hd97658.txt
            osf ehfsy $stellar gj1214.txt
            osf 2qdu8 $stellar sun.txt
            osf 45cjx $stellar l-98-59.txt
            ;;

        "surfaces")
            echo $help_surfaces
            osf 3af8u $surface andesite.dat
            osf mgae7 $surface basalt_glass.dat
            osf d92wk $surface basalt_tuff.dat
            osf pz54b $surface diorite.dat
            osf bz5jy $surface gabbro.dat
            osf hyqv5 $surface granite.dat
            osf hwkby $surface harzburgite.dat
            osf 2xpmb $surface hematite.dat
            osf 6ga5d $surface lherzolite.dat
            osf bdcte $surface lunar_anorthosite.dat
            osf 5x4fy $surface lunar_marebasalt.dat
            osf 8r3x2 $surface mars_basalticshergottite.dat
            osf 97zsh $surface mars_breccia.dat
            osf y6knb $surface norite.dat
            osf wqz48 $surface phonolite.dat
            osf qsntg $surface pyrite.dat
            osf q6ujb $surface rhyolite.dat
            osf usj7w $surface tephrite.dat
            osf aj6us $surface tholeiitic_basalt.dat
            ;;

        "thermodynamics")
            echo $help_thermo
            get_zip 4m5x8 $thermo
            ;;

        "parfiles")
            echo $help_parfiles
            osf me3uc $parfiles h2o-co2_4000-5000.par
            osf mucfd $parfiles mixture_100-50000.par
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
