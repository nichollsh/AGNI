# Contains the globe module

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

"""
**Main globe module for handling multiple 1D-atmospheric columns**
"""
module globe

    # System libraries
    using Printf
    using Logging

    # Local files
    import ..atmosphere
    import ..phys


end
