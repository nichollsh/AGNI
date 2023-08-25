# Not for direct execution

module plot_pt 

    # Import stuff
    using Plots

    function (atm, fname)
        """
        Plot the temperature-pressure profile for the 'atm' struct, saving it to
        a PDF file at 'fname'.

        """

        arr_T = atm.t_level[1, :]
        arr_P = atm.p_level[1, :]

        plot(arr_T, arr_P, yflip=true, xlabel="Temperature (K)", ylabel="pressure (Pa)")
        savefig("T.pdf")

    end
    

end

