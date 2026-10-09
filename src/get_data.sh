#!/usr/bin/env bash
# Download and unpack required and/or optional data
# All files can be found at https://zenodo.org/communities/proteus_framework

# When Zenodo fails, a record with a DataverseNL mirror is read from the mirror.
# Servers: ZENODO_URL (default https://zenodo.org), DATAVERSE_URL (https://dataverse.nl)

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
ua="AGNI/1.0 ($os $arch)"

ZENODO_URL=${ZENODO_URL:-https://zenodo.org}
DATAVERSE_URL=${DATAVERSE_URL:-https://dataverse.nl}

# Check internet connectivity
function reachable {
    header=$(wget --user-agent "'$ua'" --spider -S "$1" 2>&1 | grep "HTTP")
    [[ $header == *"HTTP/1.1 2"* || $header == *"HTTP/1.1 3"* || $header == *"response 2"* || $header == *"response 3"* ]]
}
if ! reachable "$ZENODO_URL"; then
    echo "WARNING: Failed to establish a connection to Zenodo, using the DataverseNL mirror"
    use_mirror=1
    if ! reachable "$DATAVERSE_URL"; then
        echo "ERROR: Failed to establish a connection to Zenodo and DataverseNL"
        exit 1
    fi
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

# DataverseNL mirror of a Zenodo record; add a pair here when a mirror is published
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
        *)
            echo "ERROR: Failed to download $1. It has no DataverseNL mirror." >&2
            return 1
            ;;
    esac
}

# Download a URL to a file, refusing an HTML page (an error or bot-check page, not data)
function fetch {
    # $1 = url
    # $2 = target file path
    hdr=$(wget --user-agent "'$ua'" -S -qO "$2" "$1" 2>&1)
    if [ $? -ne 0 ] || [[ ! -f "$2" ]]; then
        rm -f "$2"
        return 1
    fi
    if echo "$hdr" | grep -i "^ *content-type:" | tail -1 | grep -qi "text/html"; then
        echo "ERROR: $1 returned an HTML page, not data"
        rm -f "$2"
        return 1
    fi
    return 0
}

# Generic single file from Zenodo record
function zenodo {
    # $1 = Zenodo identifier for Record
    # $2 = target folder (on disk)
    # $3 = target filename (on disk and in Record)

    tgt="$2/$3"
    mkdir -p $2
    if [ -z "$use_mirror" ]; then
        echo "    zenodo/$1 > $tgt"
        fetch "$ZENODO_URL/records/$1/files/$3" $tgt && return 0
        echo "Trying again to download the file"
        sleep 1
        fetch "$ZENODO_URL/records/$1/files/$3" $tgt && return 0
        echo "ERROR: Failed to download $1 from Zenodo"
    fi

    doi=$(mirror_doi $1) || exit 1
    echo "    dataverse/$doi > $tgt"
    index="$DATAVERSE_URL/api/datasets/:persistentId/dirindex?persistentId=doi:$doi"
    id=$(wget --user-agent "'$ua'" -qO- "$index" | sed -n "s|.*datafile/\([0-9]*\)\">${3//./\\.}</a>.*|\1|p")
    if [ -z "$id" ] || ! fetch "$DATAVERSE_URL/api/access/datafile/$id" $tgt; then
        echo "ERROR: Failed to download $3 from DataverseNL doi:$doi"
        exit 1
    fi
    return 0
}

# Wrapper around unzip command, which first tests the archive
function unzip_wrap {
    # $1 = zip file path
    # $2 = target folder to unzip into

    if ! unzip -tq $1 > /dev/null; then
        echo "ERROR: $1 is not a valid zip archive"
        rm -f $1
        return 1
    fi

    # Exclude the readme and the DataverseNL manifest if the zip has them
    exclude=""
    for exclude_file in _readme.txt MANIFEST.TXT; do
        if unzip -l $1 | grep -q " $exclude_file$"; then
            exclude="$exclude $exclude_file"
        fi
    done

    unzip -oq $1 -d $2 ${exclude:+-x $exclude}
    rm $1

    return 0
}

# Get whole Zenodo record as Zip, and extract the files
function zenodo_all {
    # $1 = Zenodo identifier for Record
    # $2 = target folder (on disk) to extract files into

    tgt="$2/$1.zip"
    mkdir -p $2
    if [ -z "$use_mirror" ]; then
        echo "    zenodo/$1 > $tgt"
        fetch "$ZENODO_URL/api/records/$1/files-archive" $tgt && unzip_wrap $tgt $2 && return 0
        echo "Trying again to download the file"
        sleep 1
        fetch "$ZENODO_URL/api/records/$1/files-archive" $tgt && unzip_wrap $tgt $2 && return 0
        echo "ERROR: Failed to download $1 from Zenodo"
    fi

    doi=$(mirror_doi $1) || exit 1
    echo "    dataverse/$doi > $tgt"
    url="$DATAVERSE_URL/api/access/dataset/:persistentId/?persistentId=doi:$doi"
    if ! (fetch "$url" $tgt && unzip_wrap $tgt $2); then
        echo "ERROR: Failed to download $1 from DataverseNL doi:$doi"
        exit 1
    fi
    return 0
}

# Get a zip file from within a Zenodo record, and extract it
function get_zip {
    # $1 = Zenodo record
    # $2 = target folder on disk
    # $3 = name of zip file in the Zenodo record

    zenodo $1 $2 $3
    unzip_wrap "$2/$3" $2 || exit 1
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
            get_zip 21390786 $thermo gases.zip
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
