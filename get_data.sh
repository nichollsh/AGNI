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

# Make basic data folders
mkdir -p $res
mkdir -p "$res/spectral_files"
mkdir -p "$res/stellar_spectra"

# Help strings
help_basic="Get the basic data required to run the model"
help_highres="Get a spectral file with many high-resolution opacities"
help_steam="Get pure-steam spectral files"
help_stellar="Get a collection of stellar spectra"
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
"

# Generic OSF downloader function
function osf {
    # $1 = OSF identifier
    # $2 = target folder
    # $3 = target filename

    # target file path
    tgt="$2/$3"

    # exists?
    if [[ -f "$tgt" ]]; then
        echo "    $1 > file already exists"
        return 0
    fi

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

            osf qmp4e $res/spectral_files/Oak/318/ Oak.sf
            osf 5fxr7 $res/spectral_files/Oak/318/ Oak.sf_k

            osf heuza $res/spectral_files/Dayspring/48/ Dayspring.sf
            osf c5jv3 $res/spectral_files/Dayspring/48/ Dayspring.sf_k

            osf b5gsh $res/spectral_files/Dayspring/256/ Dayspring.sf
            osf dn6wh $res/spectral_files/Dayspring/256/ Dayspring.sf_k

            osf 2qdu8 $res/stellar_spectra sun.txt
            ;;

        "highres")
            echo $help_highres

            osf p672d $res/spectral_files/Honeyside/4096/ Honeyside.sf
            osf ujb4z $res/spectral_files/Honeyside/4096/ Honeyside.sf_k
            ;;

        "steam")
            echo $help_steam

            osf 6rvfe $res/spectral_files/Frostflow/16/ Frostflow.sf
            osf kxve8 $res/spectral_files/Frostflow/16/ Frostflow.sf_k

            osf 9n6mw $res/spectral_files/Frostflow/48/ Frostflow.sf
            osf xfap8 $res/spectral_files/Frostflow/48/ Frostflow.sf_k

            osf mnvyq $res/spectral_files/Frostflow/256/ Frostflow.sf
            osf tzsgc $res/spectral_files/Frostflow/256/ Frostflow.sf_k
            ;;

        "stellar")
            echo $help_stellar
            osf mabp2 $res/stellar_spectra trappist-1.txt
            osf rk7mj $res/stellar_spectra eps-eri.txt
            osf agsrq $res/stellar_spectra hd97658.txt
            osf ehfsy $res/stellar_spectra gj1214.txt
            osf 2qdu8 $res/stellar_spectra sun.txt
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
