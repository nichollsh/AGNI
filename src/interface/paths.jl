# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0

"""
**Module for handling file paths and directories.**
"""
module paths

    # Load modules
    using LoggingExtras

    # AGNI root directory (constant)
    const ROOT_DIR::String = normpath(abspath(dirname(abspath(@__FILE__)), "..", ".."))
    export ROOT_DIR

    # Resources directory (constant)
    const RES_DIR::String = normpath(joinpath(ROOT_DIR, "res"))
    export RES_DIR

    # FWL_DATA folder (fall back to RES_DIR if not set)
    const FWL_DATA::String = normpath(joinpath(get(ENV, "FWL_DATA", RES_DIR)))
    export FWL_DATA

    # RAD_DIR (socrates root directory)
    const RAD_DIR::String = abspath(ENV["RAD_DIR"])
    export RAD_DIR

    """
    **Get path to other data dirs (can be overridden)**

    Arguments:
    - `name::String` name of the directory to get

    Returns:
    - `String` path to the requested directory, or `nothing` if the name is unknown.
    """
    function get_dir(name::String)::Union{String, Nothing}

        if name == "thermodynamics"
            return joinpath(RES_DIR, "thermodynamics")

        elseif name == "scattering"
            return joinpath(RES_DIR, "scattering")

        elseif name == "config"
            return joinpath(RES_DIR, "config")

        elseif name == "stellar_spectra"
            return joinpath(RES_DIR, "stellar_spectra")

        elseif name == "spectral_files"
            return joinpath(RES_DIR, "spectral_files")

        elseif name == "blobs"
            return joinpath(RES_DIR, "blobs")

        elseif name == "out"
            return joinpath(ROOT_DIR, "out")

        else
            @warn "Unknown directory name: $name"
            return nothing
        end
    end
    export get_dir

    """
    **Check if directory is 'safe' for removal**

    Arguments:
    - path::String                  the path to check

    Returns:
    - Bool                          true if the path is safe for removal
    """
    function is_safe_dir(path::String)::Bool
        # Do not allow empty paths
        isempty(path) && return false

        # Normalise path for other checks...
        path = normpath(abspath(path))

        # Contains git repo
        ispath(joinpath(path, ".git")) && return false

        # Is current working directory
        (path == pwd()) && return false

        # Is system root directory
        (path == normpath("/")) && return false

        # Is user home directory
        (path == homedir()) && return false

        # Is AGNI root directory
        (path == paths.ROOT_DIR) && return false

        # Is AGNI resources directory
        (path == paths.RES_DIR) && return false

        return true
    end
    export is_safe_dir

end
