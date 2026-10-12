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

# User agent string to identify the script to Zenodo
os=$(uname -s)
arch=$(uname -m)
ua="AGNI/1.13 ($os $arch)"

# Data servers; a record that Zenodo fails to provide is read from its DataverseNL mirror
ZENODO_URL=${ZENODO_URL:-https://zenodo.org}
DATAVERSE_URL=${DATAVERSE_URL:-https://dataverse.nl}

# Check internet connectivity
header=$(wget --user-agent "'$ua'" -T 30 -t 2 --spider -S "$ZENODO_URL" 2>&1 | grep "HTTP")
# echo $header
use_mirror=""
if ! [[ $header == *"HTTP/1.1 2"* || $header == *"HTTP/1.1 3"* || $header == *"response 2"* || $header == *"response 3"* ]]; then
    # Use the mirrors if we don't get a positive HTTP response from Zenodo
    echo "WARNING: Failed to establish a connection to Zenodo, using DataverseNL"
    echo "Response: $header"
    use_mirror=1
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
scattering=$res/scattering
refractive=$res/refractive

# Make basic data folders
mkdir -p $res
mkdir -p $spfiles
mkdir -p $stellar
mkdir -p $surface
mkdir -p $thermo
mkdir -p $parfiles
mkdir -p $scattering
mkdir -p $refractive

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
help_scattering="Get lookup data and parameter files for aerosol and cloud scattering"
help_refractive="Get refractive indices of aerosol and cloud materials, for Mie calculations"
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
    scattering
        $help_scattering
    refractive
        $help_refractive
    thermodynamics
        $help_thermo\
"

# Get the DOI of the DataverseNL mirror of a Zenodo record
function mirror_doi {
    # $1 = Zenodo identifier for Record

    case $1 in
        15799743) echo 10.34894/DX7CDY ;;
        15696415) echo 10.34894/MFHZIN ;;
        15799754) echo 10.34894/V9KQKY ;;
        15799776) echo 10.34894/NATTTQ ;;
        15799318) echo 10.34894/VRDDRJ ;;
        15721749) echo 10.34894/ES4SKG ;;
        15799474) echo 10.34894/9FGULL ;;
        15799495) echo 10.34894/JXMS4R ;;
        15799607) echo 10.34894/WDE4CC ;;
        15799652) echo 10.34894/UBSSB2 ;;
        15799731) echo 10.34894/1KKSYT ;;
        15696457) echo 10.34894/2WO9EC ;;
        15743843) echo 10.34894/K3UKBX ;;
        15806343) echo 10.34894/LZUB3T ;;
        17981836) echo 10.34894/37SUKC ;;
        15721440) echo 10.34894/BC1DEH ;;
        15880455) echo 10.34894/8ARDN5 ;;
        19294180) echo 10.34894/6Z8Y0Q ;;
        23000222) echo 10.34894/PZFHP2 ;;
        15881238) echo 10.34894/QR4XR2 ;;
        15881496) echo 10.34894/FLOPVW ;;
        * )
            echo "ERROR: Failed to download $1 from Zenodo, and it has no mirror" >&2
            return 1
            ;;
    esac
}

# Wrapper around wget command
function fetch {
    # $1 = url
    # $2 = target file path, or - to print the content

    wget --user-agent "'$ua'" -T 60 -t 3 -qO "$2" "$1"
}

# Generic single file from Zenodo record, or from its mirror if Zenodo fails
function zenodo {
    # $1 = Zenodo identifier for Record
    # $2 = target folder (on disk)
    # $3 = target filename (on disk and in Record)

    # target file path
    tgt="$2/$3"

    # target url
    url="$ZENODO_URL/records/$1/files/$3"

    rm -rf $tgt
    mkdir -p $2

    # get data, trying twice
    if [ -z "$use_mirror" ]; then
        echo "    zenodo/$1 > $tgt"
        for attempt in 1 2; do
            if fetch $url $tgt; then
                return 0
            fi
            echo "WARNING: Failed to download $1 from Zenodo (attempt $attempt of 2)"
            sleep 1
        done
    fi

    # find the file in the mirror by its name, and get it
    doi=$(mirror_doi $1) || exit 1
    echo "    dataverse/$doi > $tgt"
    url="$DATAVERSE_URL/api/datasets/:persistentId/dirindex?persistentId=doi:$doi"
    id=$(fetch "$url" - | sed -n "s|.*datafile/\([0-9]*\)\">${3//./\\.}</a>.*|\1|p")
    if [ -z "$id" ] || ! fetch "$DATAVERSE_URL/api/access/datafile/$id" $tgt; then
        echo "ERROR: Failed to download $3 from DataverseNL (doi:$doi)"
        rm -f $tgt
        exit 1
    fi

    return 0
}

# Wrapper around unzip command
function unzip_wrap {
    # $1 = zip file path
    # $2 = target folder to unzip into

    # If a file to check for is in the zip, exclude it from unzip command
    exclude_flag=""
    for exclude_file in _readme.txt MANIFEST.TXT; do
        if unzip -l $1 | grep -q $exclude_file; then
            exclude_flag="${exclude_flag:--x} $exclude_file"
        fi
    done

    unzip -oq $1 -d $2 ${exclude_flag}
    retcode=$?
    rm $1

    return $retcode
}

# Get all files in a Zenodo record, except for the readme
function zenodo_all {
    # $1 = Zenodo identifier for Record
    # $2 = target folder (on disk) to download files into

    # get list of files in the record
    files=""
    if [ -z "$use_mirror" ]; then
        url="$ZENODO_URL/api/records/$1"
        files=$(fetch $url -)
        retcode=$?

        # check if command failed, and try again?
        if [ $retcode -ne 0 ]; then
            echo "WARNING: Failed to list files in $1. Issue with wget (code=$retcode)"
            echo "Trying again to list the files"
            sleep 1
            files=$(fetch $url -) || files=""
        fi
        files=$(echo "$files" | grep -o "/$1/files/[^/\"]*/content")
    fi

    # without a list of files, get the whole record as one zip from its mirror
    if [[ -z "$files" ]]; then
        tgt="$2/$1.zip"
        mkdir -p $2
        doi=$(mirror_doi $1) || exit 1
        echo "    dataverse/$doi > $tgt"
        url="$DATAVERSE_URL/api/access/dataset/:persistentId/?persistentId=doi:$doi"
        if ! fetch "$url" $tgt || ! unzip_wrap $tgt $2; then
            echo "ERROR: Failed to download $1 from DataverseNL (doi:$doi)"
            rm -f $tgt
            exit 1
        fi
        return 0
    fi
    echo "    downloading $(echo "$files" | wc -l) files"

    # download each file in turn
    for f in $files; do
        f=${f#/$1/files/}
        f=${f%/content}
        if [[ "$f" != "_readme.txt" ]]; then
            zenodo $1 $2 $f
        fi
    done

    return 0
}

# Get a zip file from within a Zenodo record, and extract it
function get_zip {
    # $1 = Zenodo record
    # $2 = target folder on disk
    # $3 = name of zip file in the Zenodo record

    zenodo $1 $2 $3
    unzip_wrap "$2/$3" $2
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
            anyspec Dayspring 48
            anyspec Dayspring 16

            zenodo 17981836 $stellar sun.txt

            handle_request "thermodynamics"
            handle_request "scattering"
            handle_request "refractive"

            zenodo 15806626 $parfiles h2o-co2_4000-5000.par
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

            # muscles spectra
            zenodo_all 15721440 $stellar

            # modern sun spectrum
            zenodo 17981836 $stellar sun.txt
            ;;

        "surfaces")
            echo $help_surf_standard
            zenodo_all 15880455 $surface    # Hammond+24 RELAB data
            ;;

        "surfaces_extended")
            echo $help_surf_extended
            get_zip 15881238 $surface ecostress.zip
            get_zip 15881496 $surface lavaworld.zip
            ;;

        "thermodynamics")
            echo $help_thermo
            get_zip 23278603 $thermo gases.zip
            ;;

        "parfiles")
            echo $help_parfiles
            rec="15806626"
            zenodo $rec $parfiles h2o-co2_4000-5000.par
            zenodo $rec $parfiles mixture_100-50000.par
            ;;

        "scattering")
            echo $help_scattering
            zenodo_all 19294180 $scattering
            ;;

        "refractive")
            echo $help_refractive
            zenodo_all 23000222 $refractive
            ;;

        *)
            echo "$help"
            ;;
    esac
    return 0
}

handle_request $@

exit 0
