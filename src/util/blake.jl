# This file is part of AGNI. License is Apache-2.0: https://apache.org/licenses/LICENSE-2.0

"""
**Module for performing hashing with BLAKE2b.**

Blake2b is a cryptographic hash function, which can be used to verify the integrity
of files. More information about BLAKE2b here: https://www.blake2.net/

Not supported on all platforms.
"""
module blake

    # Import modules
    using LoggingExtras

    # Load the BLAKE2b implementation
    include("blake2b.jl")
    import .Blake2b: compute2b

    """
    **Calculate the BLAKE2b hash for a string.**

    Converts to/from bytes using UTF-8 encoding.

    Arguments:
    - `data::String`       The string to hash.

    Returns:
    - `hash::String`       The BLAKE2b hash of the string.
    """
    function hash_string(data::String)::String
        return join(string.(compute2b(Vector{UInt8}(data)), base=16, pad=2), "")
    end
    export hash_string


    """
    **Calculate the BLAKE2b hash for a file.**

    Arguments:
    - `fpath::String`       Path to the file to hash.

    Returns:
    - `String`: The BLAKE2b hash of the file.
    """
    function hash_file(fpath::String)::String
        # check file exists
        if !isfile(fpath)
            @warn "File not found '$fpath'"
            return "FILE_NOT_FOUND:$fpath"
        end

        # compute the hash
        return hash_string(read(fpath,String))
    end
    export hash_file

    """
    **Validate the integrity of a file using its BLAKE2b hash.**

    Assumes that there's a `.chk` file in the same folder.

    Arguments
    - `fpath::String`       Path to the file to validate

    Returns
    - `Bool`                File is valid.
    """
    function valid_file(fpath::String)::Bool
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
    export valid_file

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
