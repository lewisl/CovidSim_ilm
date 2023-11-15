####################################################################################################
#   accessing real COVID-19 data from Johns Hopkins provided data 
####################################################################################################

"""
 for Johns Hopkins US actual data
"""
struct Col_ref
    date::String
    col::Int64
end


"""
    function get_real_data(;series="confirmed")

You must clone the Johns Hopkinds COVID-19 tracking data repository on Github at:

```http://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data```

Place the cloned repository in the directories referenced in this function
"""
function get_real_data(;series="confirmed")
    us_confirmed = string("/Users/lewis/Dropbox/Online Coursework/Covid References/COVID-19/csse_covid_19_data/", 
                    "csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")
    us_dead = string("/Users/lewis/Dropbox/Online Coursework/Covid References/COVID-19/csse_covid_19_data/",
                    "csse_covid_19_time_series/time_series_covid19_deaths_US.csv")
    if series == "confirmed"
        dat = read_actual(datapath = us_confirmed, series=series)   # returns (dat = dat, first = firstdata, last=lastdata)
    elseif series == "dead"
        dat = read_actual(datapath = us_dead, series=series)
    else
        @warn "Series must be \"confirmed\" or \"dead\", got $series. Returning empty array."
        return
    end
    return (dat = dat, first=first, last=last)
end


function loc2df(;confdat=nothing, deaddat=nothing, loc=53033, days="all")
    if isnothing(confdat) && isnothing(deaddat)
        @warn "A datafile must be provided for confirmed infections (confdat=<array>), deaths (deaddat=<array>), or both."
        return
    end

    fips = 0
    county_state = ""
    infected = []
    dead = []
    confdays = deaddays = 0:0
    if isa(loc, Integer)  || typeof(loc) <: AbstractString
        search_item = loc
    else
        @warn "loc must be a valid FIPS number for a US county or a \"<county>, <state>, US\" string"
        return
    end

    if !isnothing(confdat)
        confdays = days == "all" ? (1:(size(confdat,2) - 2)) : days
        infected = select_row(search_item, confdat, confdays)
    end

    if !isnothing(deaddat)
        deaddays = days == "all" ? (1:(size(deaddat,2) - 2)) : days
        dead = select_row(search_item, deaddat, deaddays)
    end

    if !isempty(infected) && !isempty(dead)
        try
            df = DataFrame(infected = reshape(infected, :, 1)[:,1], dead = reshape(dead, :, 1)[:,1])
        catch
            @warn "number of days did not match for infected and dead data series"
            return
        end
    elseif !isempty(infected)
        df = DataFrame(infected = reshape(infected, :, 1)[:,1])
    elseif !isempty(dead)
        df = DataFrame(dead = reshape(dead, :, 1)[:,1])
    else
        @warn "neither data series could be created"
        return
    end

    return df
end


function select_row(search_item, dat, days)
    if typeof(search_item) <: Int  # search_item is a FIPS number
        founditem = search_item .== dat[:,1]  # 1 true, rest false if good
    elseif typeof(search_item) <: AbstractString  # search_item looks like "King, Washington"
        founditem = occursin.(search_item, dat[:,2])  # 1 true, rest false if good
    end
    num = count(founditem)  # 1 true, rest false if good
    if num == 1  
        return dat[founditem, days.start+2:days.stop+2]
    elseif num > 1
        @warn "more than one occurence of $search_item was found in data array"
        return []
    else
        @warn "$search_item was not found in data array"
        return []
    end                
end


function read_actual(;datapath="", series="")
    grab = readdlm(datapath, ',', header=true)
    if series == "confirmed"   # thanks to Johns Hopkins inconsistent file formats!
        dat = grab[1][:, [5,11,12:end...]]
        hdr = grab[2][1, [5,11,12:end...]]
    elseif series == "dead"  # the dead file uses column 12 for population!
        dat = grab[1][:, [5,11,13:end...]]
        hdr = grab[2][1, [5,11,13:end...]]
    end
    firstdata = Col_ref(hdr[3], 3)
    lastdata = Col_ref(hdr[end], length(hdr))
    return (dat=dat, first=firstdata, last=lastdata)
end
