# Contains functions for reading/writing data between disk and program

# Not for direct execution
if (abspath(PROGRAM_FILE) == @__FILE__)
    thisfile = @__FILE__
    error("The file '$thisfile' is not for direct execution")
end 

module fileio 

    using Printf

    function write_pt(atmos, fname)
        
        rm(fname, force=true)

        open(fname, "w") do f
            write(f, "# pl    , tmpl \n")
            write(f, "# [bar] , [K] \n")
            for i in 1:atmos.nlev_l
                @printf(f, "%1.5e , %1.5e \n", atmos.pl[i] * 1.0e-5, atmos.tmpl[i] * 1.0e-5)
            end
        end

    end

end 
