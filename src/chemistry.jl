# Contains things relating to atmospheric chemistry

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module chemistry

    # System libraries
    using Printf
    using LinearAlgebra
    using Logging
    using LoopVectorization

    # Local files
    import ..atmosphere
    import ..phys

end # end module
