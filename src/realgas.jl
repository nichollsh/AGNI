# Contains routines for loading and working with real-gas equations of state

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end

module phys

    import ..phys

    using PCHIPInterpolation
    using LoggingExtras

    function read_aqua_pt(fpath::String)

    end

end
