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

    # Folders in res/ that get_dir resolves and AGNI_DIR_<name> can override
    const RES_NAMES = ("thermodynamics", "scattering", "refractive", "config",
                        "stellar_spectra", "spectral_files", "blobs")

    # RAD_DIR (socrates root directory)
    const RAD_DIR::String = abspath(ENV["RAD_DIR"])
    export RAD_DIR

    """
    **Resolve a folder of `res/` and the setting that placed it.**

    The folder is `AGNI_DIR_<name>` when that is set; otherwise it sits in `AGNI_DIR_res`
    when set, else in `res` (the config `[files] res_dir`), else in the `res/` of AGNI.
    Blank variables count as unset; values are made absolute (`~` is not expanded).

    Arguments:
    - `name::String` one of `RES_NAMES`
    - `res::Union{String,Nothing}` res root from the configuration, or `nothing`

    Returns:
    - `Tuple` the folder, and the setting that placed it (`nothing` for the default)
    """
    function resolve_dir(name::String; res::Union{String,Nothing}=nothing)
        # A set, non-blank variable as an absolute path without a trailing separator
        function env_dir(v::String)::Union{String, Nothing}
            dir = strip(get(ENV, v, ""))
            isempty(dir) && return nothing
            dir = abspath(dir)
            return (length(dir) > 1 && endswith(dir, "/")) ? dir[1:end-1] : dir
        end
        var = "AGNI_DIR_$name"
        dir = env_dir(var)
        isnothing(dir) || return (dir, var)
        root = env_dir("AGNI_DIR_res")
        isnothing(root) || return (joinpath(root, name), "AGNI_DIR_res")
        isnothing(res) || return (joinpath(res, name), "[files] res_dir")
        return (joinpath(RES_DIR, name), nothing)
    end

    """
    **Get path to other data dirs (can be overridden, see `resolve_dir`)**

    Arguments:
    - `name::String` name of the directory to get
    - `res::Union{String,Nothing}` res root from the configuration, or `nothing`

    Returns:
    - `String` path to the requested directory, or `nothing` if the name is unknown.
    """
    function get_dir(name::String;
                        res::Union{String,Nothing}=nothing)::Union{String, Nothing}
        name == "out" && return joinpath(ROOT_DIR, "out")
        name in RES_NAMES && return first(resolve_dir(name; res=res))
        @warn "Unknown directory name: $name"
        return nothing
    end
    export get_dir

    """
    **Check if directory is 'safe' for removal and can be written to.**

    Arguments:
    - path::String                  the path to check

    Returns:
    - Bool                          true if the path is safe
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

        # Have permissions to write to this path, if it exists
        (isdir(path) && !iswritable(path)) && return false

        return true
    end
    export is_safe_dir

    """
    **Get available disk space on file system mounted at a path**

    Arguments:
    - `path::String` the path to check

    Returns:
    - `Int64` available disk space in bytes
    """
    function get_avail_space(path::String)::Int64
        # Check if path exists - if not, default to the root directory of the file system
        if !ispath(path)
            path = "/"
        end
        return Int64(diskstat(path).available)
    end
    export get_avail_space
end
