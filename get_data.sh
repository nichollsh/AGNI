#!/bin/bash
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
res="$root/res/"
spfiles=$res/spectral_files
stellar=$res/stellar_spectra
surface=$res/surface_albedos
realgas=$res/realgas

# Make basic data folders
mkdir -p $res
mkdir -p $spfiles
mkdir -p $stellar
mkdir -p $surface
mkdir -p $realgas

# Help strings
help_basic="Get the basic data required to run the model"
help_highres="Get a spectral file with many high-resolution opacities"
help_steam="Get pure-steam spectral files"
help_stellar="Get a collection of stellar spectra"
help_surfaces="Get a collection of surface single-scattering albedos"
help_realgas="Get a real-gas EOS coefficients and lookup tables"
help="\
Helper script used to download and unpack data used to run the model.

Call structure:
    ./get_data [TARGET]

Where [TARGET] can be any of the following:
    basic
        $help_basic
    highres
        $help_highres
    steam
        $help_steam
    stellar
        $help_stellar
    surfaces
        $help_surfaces
    realgas
        $help_realgas
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

    return 0
}

# Handle request for downloading a group of data
function handle_request {
    case $1 in
        "basic")
            echo $help_basic

            osf qmp4e $spfiles/Oak/318/ Oak.sf
            osf 5fxr7 $spfiles/Oak/318/ Oak.sf_k

            osf heuza $spfiles/Dayspring/48/ Dayspring.sf
            osf c5jv3 $spfiles/Dayspring/48/ Dayspring.sf_k

            osf b5gsh $spfiles/Dayspring/256/ Dayspring.sf
            osf dn6wh $spfiles/Dayspring/256/ Dayspring.sf_k

            osf 2qdu8 $res/stellar_spectra sun.txt
            ;;

        "highres")
            echo $help_highres

            osf p672d $spfiles/Honeyside/4096/ Honeyside.sf
            osf ujb4z $spfiles/Honeyside/4096/ Honeyside.sf_k
            ;;

        "steam")
            echo $help_steam

            osf 6rvfe $spfiles/Frostflow/16/ Frostflow.sf
            osf kxve8 $spfiles/Frostflow/16/ Frostflow.sf_k

            osf 9n6mw $spfiles/Frostflow/48/ Frostflow.sf
            osf xfap8 $spfiles/Frostflow/48/ Frostflow.sf_k

            osf mnvyq $spfiles/Frostflow/256/ Frostflow.sf
            osf tzsgc $spfiles/Frostflow/256/ Frostflow.sf_k
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

        "realgas")
            echo $help_realgas

            # AQUA PT lookup
            zipfile=".aqua_temp.zip"
            osf uqrdx $realgas $zipfile
            unzip -oq $realgas/$zipfile -d $realgas
            rm $realgas/$zipfile

            ;;

        *)
            echo "$help"
            return 0
            ;;
    esac
    return 0
}

handle_request $1

echo "Done!"

exit 0
