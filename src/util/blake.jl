# This file is part of AGNI. License is GPL-3.0: https://www.gnu.org/licenses

"""
**Module for performing hashing with BLAKE2b.**

Blake2b is a cryptographic hash function, which can be used to verify the integrity
of files. More information about BLAKE2b here: https://www.blake2.net/

Not supported on all platforms.
"""
module blake

    # Import modules
    using LoggingExtras

    # Path to executable blob
    const exec_path::String = abspath(dirname(@__FILE__),
                                        "..", "..", "res", "blobs", "b2sum-amd64-linux")


    """
    **Calculate the BLAKE2b hash for a file.**

    Only works on AMD64 Linux.

    Arguments:
    - `fpath::String`       Path to the file to hash.
    - `quiet::Bool=false`   Suppresses warnings about missing files.

    Returns:
    - `String`: The BLAKE2b hash of the file.
    """
    function hash_file(fpath::String; quiet::Bool=false)::String
        if !Sys.islinux()
            return "ONLY_SUPPORTED_ON_LINUX"
        end
        if Sys.ARCH != :x86_64
            return "ONLY_SUPPORTED_ON_AMD64"
        end
        if !isfile(fpath)
            quiet || @warn "File not found '$fpath'"
            return "FILE_NOT_FOUND:$fpath"
        end
        content = strip(read(`$exec_path $fpath`, String))
        return split(content, " ")[1]
    end
    export hash_file

    """
    **Validate the integrity of a file using its BLAKE2b hash.**

    Assumes that there's a `.chk` file in the same folder.

    Arguments
    - `fpath::String`       Path to the file to validate
    - `quiet::Bool=false`   Suppresses warnings

    Returns
    - `Bool`                File is valid.
    """
    function valid_file(fpath::String; quiet::Bool=false)::Bool
        # Return true if unsupported
        if !Sys.islinux() || Sys.ARCH != :x86_64
            quiet || @debug "Skipping integrity check for '$fpath'"
            return true
        end

        # Get the actual hash
        hash_obs = hash_file(fpath)

        # Get the expected hash
        cpath::String = fpath*".chk"
        if !isfile(cpath)
            quiet || @warn "File not found '$cpath'"
            return false
        end
        hash_exp = strip(read(fpath*".chk", String))

        # Check if equal
        if Bool(hash_exp != hash_obs)
            quiet || @warn "File '$fpath' failed integrity check"
            quiet || @warn "    obs = '$hash_obs'"
            quiet || @warn "    exp = '$hash_exp'"
            return false
        end
        return true
    end
    export valid_file

end # end module

# If executed directly...
if abspath(PROGRAM_FILE) == @__FILE__
    import .blake
    if length(ARGS) != 1
        @error("Invalid arguments: $(ARGS)")
    elseif !Sys.islinux()
        @error("File hashing is only supported on Linux")
    elseif Sys.ARCH != :x86_64
        @error("File hashing is only supported on AMD64 architecture")
    else
        @info("Computing BLAKE2b hash...")
        hash_obs = blake.hash_file(ARGS[1])
        @info("$hash_obs")
    end
end
