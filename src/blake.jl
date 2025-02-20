# Contains module and functions for performing hashing with BLAKE2b

module blake

    # Import modules
    using LoggingExtras

    # Path to executable blob
    const exec_path::String = abspath(dirname(@__FILE__),
                                        "..", "res", "blobs", "b2sum-amd64-linux")


    # Hashing only works on AMD64 Linux
    function _is_linux()
        if !Sys.islinux()
            @warn "Cannot determine file hash - only supported on Linux"
            return false
        end
        return true
    end

    """
    **Calculate the BLAKE2b hash for a file.**

    Only works on AMD64 Linux.
    """
    function hash_file(fpath::String)::String
        if !_is_linux()
            return "ONLY_SUPPORTED_ON_LINUX"
        end
        if !isfile(fpath)
            @warn "File not found '$fpath'"
            return "FILE_NOT_FOUND:$fpath"
        end
        content = strip(read(`$exec_path $fpath`, String))
        return split(content, " ")[1]
    end

    """
    **Validate the integrity of a file using its BLAKE2b hash.**

    Assumes that there's a `.chk` file in the same folder.
    """
    function valid_file(fpath::String)::Bool
        # Return true if unsupported
        if !_is_linux()
            @debug "Skipping integrity check for '$fpath'"
            return true
        end

        # Get the actual hash
        hash_obs = hash_file(fpath)

        # Get the expected hash
        cpath::String = fpath*".chk"
        if !isfile(cpath)
            @warn "File not found '$cpath'"
            return false
        end
        hash_exp = strip(read(fpath*".chk", String))

        # Check if equal
        if Bool(hash_exp != hash_obs)
            @warn "File '$fpath' failed integrity check"
            @warn "    obs = '$hash_obs'"
            @warn "    exp = '$hash_exp'"
            return false
        end
        return true
    end

end # end module

# If executed directly...
if abspath(PROGRAM_FILE) == @__FILE__
    import .blake
    if length(ARGS) != 1
        @error("Invalid arguments: $(ARGS)")
    else
        @info("Computing BLAKE2b hash...")
        hash_obs = blake.hash_file(ARGS[1])
        @info("$hash_obs")
    end
end
