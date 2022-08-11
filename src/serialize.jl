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


function modeldef_to_yaml(ndays::Int, locales::Vector{Int};  
            day1,
            dovax,
            paramdir = "../sample_parameters",
            geofilename = "../data/geo2data.csv", 
            socialfilename = "socialparams.yml",
            scheddir="vaccine_schedule",
            vaccinefilename = "vaccines.yml",
            variantfilename = "variants.yml",
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

            elseif key == "geofile"
                filepath = pfname
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

    yaml_model = YAML.load_file(readpathstr)

end