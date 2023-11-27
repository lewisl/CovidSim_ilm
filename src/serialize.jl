#####################################################################
#  Output simulation data in various kinds of text files
#####################################################################


"""
Utility function used by series_to_csv, popdat_to_csv, and modeldef_to_yaml
"""
function setpathstr(;pathstr="", idstr="", overwrite=false, usetimestamp=true, basedir=:current)
    datestr = usetimestamp ? string(round(Dates.now(), Dates.Minute(1))) : ""
    datestr = replace(datestr, ":" => "-") # no colons in filenames!

    basedirstr =    if basedir === :current
                        pwd()
                    elseif basedir === :home
                        homedir()
                    elseif basedir === :none
                        ""
                    else
                        throw(DomainError(basedir, "Argument must be :current, :home, or :none"))
                    end
    writepathstr = joinpath(basedirstr, pathstr)

    if !isdir(writepathstr)
        throw(ErrorException("FATAL: Path $writepathstr does not exist"))
    end

    return writepathstr, datestr
end

"""
Utility function to change the type of the keys in the input dict.
"""
function change_key_type(d; f=Symbol)
    Dict(f(k) =>
        if !(typeof(v) <: AbstractDict)
            v
        else
            change_key_type(v, f=f)
        end
         for (k, v) in d)
end


"""
     **series\\_to\\_csv(series; pathstr="", idstr="", locale=0, overwrite=false, usetimestamp=true, basedir=:current)**

Outputs simulation history series as csv. Each locale results in 
two csv files: one for the cum (cumulative) values and 
another csv file for the new (day-to-day net change) values.

Optional inputs as named parameters:
- pathstr: directory--relative to the basedir--where the file will be written. It is not necessary to end the path with a '/'.
- basedir: Use :current to save to the current directory or make path relative to the current directory. Use :home to use your home directory as the base directory. Use :none to use pathstr only.
- idstring: an optional string to describe the series that will become part of the filename. It could be a name for the scenario based on your inputs for vaccination, social distancing, etc.
- locale: use to select an individual locale from the data. Default is 0, a signal to output all locales.
- overwrite: if true, an existing file will be over-written without recourse. Default is false.
- usetimestamp: if true, a part of the filename will be date-time. You may use both an idstring and a timestamp. Default is true to enable writing output files with unique names.

"""
function series_to_csv(dat; pathstr="", idstr="", locale=0, overwrite=false, usetimestamp=true, basedir=:current)
    
    writepathstr, datestr = setpathstr(pathstr=pathstr, idstr=idstr, overwrite=overwrite, 
                              usetimestamp=usetimestamp, basedir=basedir)

    if locale == 0
        sourcelocales = keys(dat)
    else
        if in(locale, keys(dat))
            sourcelocales = [locale] # make scalar iterable
        else
            throw(DomainError(locale, "Locale $locale not in datafile"))
        end
    end

    for loc in sourcelocales
        for kind in (:cum, :new)

            fname = join(filter(!=(""), ("series", idstr, kind, "loc", loc, datestr, ".csv")), "_", "")
            filepathstr = joinpath(writepathstr, fname)    

            if (isfile(filepathstr)) & (!overwrite)
                throw(ErrorException("FATAL: Argument overwrite set to false: can't overwrite existing file"))
            end 

            CSV.write(filepathstr, getfield(dat[loc], kind))

        end
    end
    
end


"""
     **popdat\\_to\\_csv(dat; pathstr="", idstr="", locale=0, overwrite=false, usetimestamp=true, basedir=:current)**

Outputs simulation population data as csv. Note that this "popdat" refers to a moment in time during the simulation and includes limited history data.

Optional inputs as named parameters:
- pathstr: directory--relative to the basedir--where the file will be written. It is not necessary to end the path with a '/'.
- basedir: Use :current to save to the current directory or make path relative to the current directory. Use :home to use your home directory as the base directory. Use :none to use pathstr only.
- idstring: an optional string to describe the series that will become part of the filename. It could be a name for the scenario based on your inputs for vaccination, social distancing, etc.
- locale: use to select an individual locale from the data. Default is 0, a signal to output all locales.
- overwrite: if true, an existing file will be over-written without recourse. Default is false.
- usetimestamp: if true, a part of the filename will be date-time. You may use both an idstring and a timestamp. Default is true to enable writing output files with unique names.

"""
function popdat_to_csv(dat; pathstr="", idstr="", overwrite=false, usetimestamp=true, basedir=:current)

    writepathstr, datestr = setpathstr(pathstr=pathstr, idstr=idstr, overwrite=overwrite, 
                                       usetimestamp=usetimestamp, basedir=basedir)

    if !isdir(writepathstr)
        throw(ErrorException("FATAL: Path $writepathstr does not exist"))
    end
    
    for loc in keys(dat)

            fname = join(filter(!=(""), ("popdat", idstr, "loc", loc, datestr, ".csv")), "_", "")
            filepathstr = joinpath(writepathstr, fname)    

            if (isfile(filepathstr)) & (!overwrite)
                throw(ErrorException("FATAL: Argument overwrite set to false: can't overwrite existing file"))
            end 

            CSV.write(filepathstr, dat[loc])

    end
    
end

#=
A model is a named tuple of all of the data structures that are created and initialized before a
simulation is run. You'll never use these types explicitly. The data structures are created by 
function setup_model and functions it calls. Data for parameters are loaded from yaml files. The
content of the parameter files are converted to appropriate Juia data structures.

A model contains these elements of the type shown:
    :ndays              => Int64
    :day1               => Date
    :locales            => Vector{Int64} 
    :dat                => NamedTuple{(:popdat, :agegrp_idx)}
        :popdat         => Dict{Int64, LazyTable} # the int is a locale identifier
        :LazyTable      => contains columns:
                           :status      => Symbol
                           :agegrp      => Symbol
                           :cond        => Symbol
                           :duration    => Int64 
                           :variant     => Vector{Symbol}
                           :sickday     => Vector{Symbol}
                           :recovday    => Vector{Int64}
                           :deadday     => Vector{Int64} 
                           :ring        => Int64
                           :sdcase      => Int64 ???
                           :vaxstatus   => Symbol 
                           :vaxrcvd     => Vector{Symbol} 
                           :vaxday      => Vector{Int64}
                           :tested      => Bool
                           :testday     => Int64
                           :quar        => Bool
                           :quarday     => Int64
    :series, 
    :geo, 
    :progressionset, 
    :vaxset, 
    :vaxschedset, 
    :infectset, 
    :social, 
    :trvec, 
    :variantlist, 
    :vaxlist, 
    :indoor_seq, 
    :seriescolnames
=#

##############################################################################################
#  functions called by model_to_yaml
##############################################################################################

function geo_to_yaml!(io, geo, pad)
    colnames = Tables.columnnames(geo)
    num_names = length(colnames)
    # write header row
    for i in 1:(num_names-1)
        write(io, string(pad, colnames[i], ','))
    end
    write(io, string(pad, colnames[num_names], '\n'))
    # write rows
    for row in geo
        print(io, pad)
        for (i, v) in enumerate(values(row))
            print(io, v)
            if i < num_names
                write(io, ", ")
            else
                write(io, '\n')
            end
        end
    end
end


function vaccinefile_to_yaml!(io, vaxset, pad)
    for (k, v) in vaxset
        write(io, pad, k, ":\n")  # string(pad, k, ":\n")
        for name in fieldnames(Vaccineparams)
            if name === :infectfactor
                write(io, repeat(pad,2), name, ":\n")
                for (k, v) in getfield(v, name)
                    write(io, repeat(pad,3), string(k, ": ", v, '\n'))
                end
            elseif name === :effectiveness
                write(io, repeat(pad,2), name, ":\n")
                for (k, v) in getfield(v, name)
                    write(io, repeat(pad,3), string(k, ": \n"))
                    for (k2, v2) in v
                        write(io, string(repeat(pad,4), k2, ": ", v2, "\n"))
                    end
                end
            else
                write(io, string(repeat(pad,2), name, ": ", getfield(v, name), '\n'))
            end
        end
    end
end


function social_to_yaml!(io, social, pad)
    for name in fieldnames(SocialParams)
        write(io, pad)
        write(io, name, ": ")
        if (name === :contactfactors) | (name === :touchfactors)
            arr = getfield(social, name)
            mapfunc = name === :contactfactors ? mapcondition : maptouch
            write(io, '\n')
            for col in 1:size(arr, 2)  # columns = agegrps
                write(io, repeat(pad,2), string(mapagegrp(col), ": \n"))
                write(io, string(pad, "  ", "  ", '{'))
                for row in 1:size(arr, 1)   # rows = conditions
                    write(io, string(mapfunc(row), ": ", arr[row, col], ", "))
                end
                write(io, "}\n")
            end
        else
            print(io, getfield(social, name))
            write(io, '\n')
        end
    end
end


function variants_to_yaml!(io, infectset, progressionset, pad)
    pgset = progressionset
    infset = infectset
    for variant in keys(pgset)
        write(io, string("  ", variant, ":\n"))
        pgvar = pgset[variant]
        infvar = infset[variant]

        write(io, string(repeat(pad,2), "spread:\n"))
        write(io, string(repeat(pad,3), "sendrisk:  ["))
            for val in infvar.sendrisk
                write(io, string(val, ", "))
            end
            write(io, "]\n")
        #
        write(io, string(repeat(pad,3), "recvrisk:  ["))
            for val in infvar.recvrisk
                write(io, string(val, ", "))
            end
            write(io, "]\n")
        #
        write(io, string(repeat(pad,3), "basemultiplier: ", string(infvar.basemultiplier, "\n")))

        write(io, string(repeat(pad,2), "immunity:\n"))
            write(io, string(repeat(pad,3), "recovery_immunity:\n"))
            for (k, v) in infvar.recovery_immunity
                write(io, string(repeat(pad,4), k, ": ", v, "\n"))
            end
        #
        write(io, string(repeat(pad,3), "immunehalflife: ", string(infvar.immunehalflife, "\n")))

        write(io, string(repeat(pad,2), "progression_tree:\n"))
        for age in fieldnames(typeof(pgvar.tree))
            agetree = getproperty(pgvar.tree, age)
            write(io, repeat(pad,3), age, ":\n")
            for (duration, matrix) in sort(agetree) 
                write(io, repeat(pad, 4), string(duration, ":\n"))
                for i in eachindex(matrix[:,1])  # countup by row number--need row num to look up condition
                    write(io, repeat(pad, 5), mapcondition(i), ": [")  # start a flow vector
                    for val in matrix[i,:]  # write the column values of each row
                        write(io, string(val, ", "))
                    end
                    write(io, "]\n")  # close the flow vector
                end
            end
        end

        write(io, string(repeat(pad,2), "progression_factors:\n"))
            if !isnothing(pgvar.factors.riskadjust)  # no else: we can leave this key out if its value is empty
                write(io, repeat(pad,3), string("riskadjust: ["))
                for value in pgvar.factors.riskadjust
                    write(io, string(value), ", ")
                end
                write(io, "]\n") # close the flow vector
            end   
            write(io, repeat(pad,3), string("vaxhalflifeadjust:", "\n"))
            for (k,v) in pgvar.factors.vaxhalflifeadjust
                Base.write(io, repeat(pad, 4), string(k, ": ", v, "\n"))
            end

        write(io, "\n")
    end
end


function vaxsched_to_yaml!(io, allvaxscheds, pad)
    for (schedname,vaxsched) in allvaxscheds  # vaxsched is an instance of struct Vaxsched
        write(io, repeat(pad,2), string(schedname, ":\n"))
        for topfield in fieldnames(typeof(vaxsched))   # topfield is a Symbol of a field in vaxsched
            if topfield == :vaxesincluded   # this field is a Dict of struct
                write(io, repeat(pad, 3), string("vaxesincluded", ": \n"))
                vaxesincluded = getproperty(vaxsched, :vaxesincluded) # this is a dict of struct Vaxinclude
                for (vax, props) in vaxesincluded   # write all the props of vaxesincluded
                    write(io, repeat(pad, 4), string(vax, ": \n"))
                    for vp in fieldnames(typeof(props))
                        write(io, repeat(pad, 5), string(vp, ": ", getproperty(props, vp), "\n"))
                    end
                end
            elseif topfield == :dayrange  # write the range as a flow vector
                startday = string(getproperty(allvaxscheds[schedname], topfield)[1])  # start of range
                endday = string(getproperty(allvaxscheds[schedname], topfield)[end])  # end of range
                write(io, repeat(pad, 3), topfield, ": [", startday, ", ", endday, "]\n")
            elseif topfield == :filtervec  # must convert :val to "val" only to convert it back at load time!
                write(io, repeat(pad,3), topfield, ": ")
                # for val in allvaxscheds[schedname].filtervec
                write(io, string([string(x) for x in allvaxscheds[schedname].filtervec]), "\n")
            else  # all other fields
                write(io, repeat(pad, 3), string(topfield, ": ", getproperty(vaxsched, topfield), "\n"))
            end
        end
    end
end


function model_to_yaml(model; pathstr="", idstr="", overwrite=false, usetimestamp=true, basedir=:current )

    # setup output
        writepathstr, datestr = setpathstr(pathstr=pathstr, idstr=idstr, overwrite=overwrite,
            usetimestamp=usetimestamp, basedir=basedir)

        fname = join(filter(!=(""), ("modeldef", idstr, datestr, ".yml")), "_", "")
        filepathstr = joinpath(writepathstr, fname)

        if (isfile(filepathstr)) & (!overwrite)
            throw(ErrorException("FATAL: Argument overwrite set to false: can't overwrite existing file"))
        end

    pad = "  "  # 2 spaces

    io = IOBuffer()

    # write output
    scalars = Dict("day1" => string(model.day1), "ndays" => model.ndays, "locales" => model.locales,
                "dovax"=>model.dovax)

    YAML.write(io, scalars)
    write(io, "\n")
    
    write(io, string("geofile", ": |\n"))
        geo_to_yaml!(io, model.geo, pad)
        write(io, "\n")

    write(io, string("socialfile", ": |\n"))
        social_to_yaml!(io, model.social, pad)
        write(io, "\n")

    write(io, string("variantfile", ": |\n"))
        variants_to_yaml!(io, model.infectset, model.progressionset, pad)
        write(io, "\n")

    write(io, string("vaccinefile", ": |\n"))
        vaccinefile_to_yaml!(io, model.vaxset, pad)
        write(io, "\n")

    write(io, string("vaxscheds", ": |\n"))
        vaxsched_to_yaml!(io, model.vaxschedset, pad)
        write(io, "\n")

    write(io, "\n\n")
    flush(io)

    #write IOBuffer to the modeldef file
        seekstart(io)
        write(filepathstr, io)
        close(io)
    return filepathstr
end


function modelinputs_to_yaml(ndays::Int, locales::Vector{Int};  
            # for inputs
            day1,
            dovax,
            paramdir = "../sample_parameters",
            geofilename = "geo2data.csv", 
            socialfilename = "socialparams.yml",
            scheddir="vaccine_schedule",
            vaccinefilename = "vaccines.yml",
            variantfilename = "variants.yml",
            # for outputs
            pathstr="", idstr="", overwrite=false, usetimestamp=true, basedir=:current  
        )

    writepathstr, datestr = setpathstr(pathstr=pathstr, idstr=idstr, overwrite=overwrite, 
                              usetimestamp=usetimestamp, basedir=basedir)

    scalars = Dict("day1"=>string(day1), "ndays"=>ndays, "locales"=>locales, "dovax"=>dovax)

    parameterfiles = Dict("geofile"=>geofilename, "socialfile"=>socialfilename, "scheddir"=>scheddir,
                          "vaccinefile"=>vaccinefilename, "variantfile"=>variantfilename)

    io = IOBuffer()

    YAML.write(io, scalars)
    write(io, "\n")

    pad = "    "

    # write modeldef components to an IOBuffer
        for (key, pfname) in parameterfiles
            if key == "scheddir"
                write(io, string("vaxscheds", ": |\n\n"))
                fnames = readdir(joinpath(paramdir, scheddir), join=true)
                for filepath in fnames
                    schedname = basename(splitext(filepath)[1])
                    write(io, string(pad, schedname, ": |\n"))
                    for l in readlines(filepath)
                        write(io, string(pad, pad, l, "\n"))
                    end     
                    write(io, "\n\n")   
                end
                continue   # nothing left to do-->skip rest of loop body and get next (key, pfname)

            # elseif key == "geofile"
            #     filepath = pfname
            else
                filepath = joinpath(paramdir, pfname)
            end

            write(io, string(key, ": |\n"))
            for l in readlines(filepath)
                write(io, string(pad, l, "\n"))
            end
            write(io, "\n\n")
            
        end
        flush(io)
    

    # write IOBuffer to the modeldef file
        seekstart(io)
        fname = join(filter(!=(""), ("modeldef", idstr, datestr, ".yml")), "_", "")
        filepathstr = joinpath(writepathstr, fname)    

        if (isfile(filepathstr)) & (!overwrite)
            throw(ErrorException("FATAL: Argument overwrite set to false: can't overwrite existing file"))
        end     

        write(filepathstr, io)
        close(io)
    return filepathstr
end


"""
    yaml_to_model(fname::String; basedir=:home, pathstr="")

Deserialize the yaml of a model definition. Read in a previously saved YAML model definition to a dict that can be input to build a simulation model.
"""
function yaml_to_model(fname::String; basedir=:home, pathstr="")

    basedirstr =    if basedir === :current
                        pwd()
                    elseif basedir === :home
                        homedir()
                    elseif basedir === :none
                        ""
                    else
                        throw(DomainError(basedir, "Argument must be :current, :home, or :none"))
                    end


    readpathstr = joinpath(basedirstr, pathstr, fname)

    !isfile(readpathstr) && (throw(ErrorException("FATAL: File $writepathstr does not exist")))

    yaml_model = change_key_type(YAML.load_file(readpathstr), f=Symbol)  # return type is big ugly dict
end