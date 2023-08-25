"""
    ParseFortran

Very minimal Fortran parsing: parse Fortran Type fields,  module parameters.
"""
module ParseFortran

"""
    parse_type(fortran_file_name::AbstractString, type_name::AbstractString) 
        -> Vector[Dict(:name=>fname, :type=>ftype, :initializer=>finitializer, ...)), ...]

Parse `fortran_file_name` for fields (variable definitions) within Type `type_name`
"""
function parse_type(fortran_file_name::AbstractString, type_name::AbstractString; verbose=false)

    lines = scan_file(fortran_file_name)

    internal_lines = find_internal_lines(lines, "type", type_name)    

    return split_types(internal_lines)
end


"""
    parse_module_parameters(fortran_file_name::AbstractString; ftype="integer") 
        -> Vector[Dict(:name=>fname, :initializer=>finitializer), ...]

Parse `fortran_file_name` for parameters in module `module_name`
"""
function parse_module_parameters(fortran_file_name::AbstractString, module_name::AbstractString)
    lines = scan_file(fortran_file_name)

    internal_lines = find_internal_lines(lines, "module", module_name)

    internal_types = split_types(internal_lines)

    parameter_types = [
        l for l in internal_types
            if l[:parameter] && !contains(lowercase(l[:initializer]), "reshape") 
    ]

    return parameter_types
end

"scan a Fortran file, remove comments, merge continuation lines"
function scan_file(fortran_file_name::AbstractString)

    output = []
    lnbuf = ""
    open(fortran_file_name) do file
        for ln in eachline(file)
            # strip comment
            idx = findfirst('!', ln)
            if !isnothing(idx)
                ln = ln[1:idx-1]
            end
            # add to buffer
            lnbuf = lnbuf*rstrip(ln)
            if !isempty(lnbuf)
                # push to output unless has a continuation
                if lnbuf[end] == '&'
                    lnbuf = lnbuf[1:end-1]
                else
                    push!(output, lnbuf)
                    lnbuf = ""
                end
            end
        end
    end

    return output
end

"find Fortran lines within a type or module definition"
function find_internal_lines(input_lines, modtype, name)    
    startidx = findfirst(
        x -> split(lowercase(x)) == lowercase.([modtype, name]),
        input_lines
    )
    !isnothing(startidx) || 
        error("$modtype $name not found")

    endidx = findfirst(
        x -> split(lowercase(x)) == lowercase.(["end", modtype, name]),
        input_lines
    )
    !isnothing(endidx) || 
        error("end $modtype $name not found")

    return input_lines[startidx+1:endidx-1]
end    

"split Fortran lines defining variables into a Dict with name, type etc"
function split_types(input_lines)
    context = []

    for ln in input_lines
        if contains(ln, "::")
            lnstrip = strip(ln)
            ftype, frest = strip.(split(lnstrip, "::"))

            ftype = filter(x -> !isspace(x), lowercase(ftype))  # remove whitespace
            ftype_split = split(ftype, ",")
            ftype = ftype_split[1]
            allocatable=false
            parameter=false
            for type_extra in ftype_split[2:end]
                if type_extra == "allocatable"
                    allocatable = true
                elseif type_extra == "parameter"
                    parameter = true
                end
            end

            if contains(frest, "=")
                fname, finitializer = strip.(split(frest, "="))
            else
                fname = frest
                finitializer = nothing
            end
        
            if contains(fname, "(")
                fname, fdims = split(fname, "(")
                fdims = "("*fdims
            else
                fdims = ""
            end

            push!(
                context,
                Dict(
                    :name=>fname,
                    :type=>ftype,
                    :initializer=>finitializer,                    
                    :dims=>fdims,
                    :allocatable=>allocatable,
                    :parameter=>parameter,
                )
            )
        end
    end

    return context
end

end
