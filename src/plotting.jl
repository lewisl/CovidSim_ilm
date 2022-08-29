#########################################################################################
# plotting.jl
#########################################################################################

gr()  # Plots backend for plotting


function cumplot(series, locale, plotcols=[:unexposed, :infectious, :recovered, :dead]; 
    days="all", geo=[], thm=:ggplot2)

    cumhist = series[locale].cum
    newhist = series[locale].new

    theme(thm, foreground_color_border=:black, 
          tickfontsize=9, gridlinewidth=1)

    !(typeof(plotcols) <: Array) && (plotcols = [plotcols])

    # the data is the 2d array cumseries
    n = length(cumhist)
    days = days == "all" ? (1:n) : days
    caldays = cumhist.calday[days] # simulation days as actual calendar dates
    cumseries = hcat(columns(getproperties(cumhist,Tuple(Symbol(plcol,"_","total") for plcol in plotcols)))...)

    # labels and annotations
    labels = [titlecase(string(col)) for col in plotcols]
    labels = reshape([labels...], 1, length(labels))
    people = if !isempty(geo)
                geo.pop[geo.fips .== locale]
             else # this will off by a tiny bit because of rounding
                getproperty(cumhist, :unexposed_total)[1] + getproperty(cumhist, :infectious_total)[1]
             end   
    # cityname = !isempty(geo) ? geo[geo[:,fips] .== locale, city][1] : ""
    died =  getproperty(cumhist, :dead_total)[end]   
    unexp = getproperty(cumhist, :unexposed_total)[end]
    infected = people - unexp    
    recovered = infected - died  
    co_pal = length(plotcols) == 2 ? [theme_palette(thm)[2], theme_palette(thm)[4]] : theme_palette(thm)
 
    xtickrange = caldays[10]:Day(180):caldays[length(caldays)-10]
    xtickfmt = Dates.format.(xtickrange, "yyyy-mm-dd")

    # the plot
    plot(   caldays, cumseries[days, :], 
            size = (700,500),
            label = labels, 
            lw=2.3,
            title = "Covid for $people people for $n days\nActive Cases for Each Day",
            xlabel = "Simulation Days",
            xticks = (xtickrange, xtickfmt),
            # xticks = caldays[10]:Day(180):caldays[length(caldays)-10],
            ylabel = "People",
            legendfontsize = 10,
            color_palette = co_pal,
            reuse = false,
            legend_position = :right,
            background_color_legend=nothing,
            foreground_color_legend=nothing
        )
    annotate!(caldays[1] + Day(6), 0.51 * ylims()[2],              # half_yscale,
            text("Died: $died\nInfected: $infected\nRecovered: $recovered\nUnexposed: $unexp", 
                11, :left))
end


function newplot(series, locale, plotcols=[:infectious]; days="all", geo=[], thm=:ggplot2)

    cumhist = series[locale].cum
    newhist = series[locale].new

    theme(thm, foreground_color_border=:black, 
            tickfontsize=9, gridlinewidth=1)

    !(typeof(plotcols) <: Array) && (plotcols = [plotcols])

    # the data and labels
    n = length(newhist)
    days = days == "all" ? (1:n) : days
    caldays = newhist.calday[days]
    newseries = hcat(columns(getproperties(newhist,Tuple(Symbol(plcol,"_","total") for plcol in plotcols)))...)

    labels = [titlecase(string(col)) for col in plotcols]
    labels = reshape([labels...], 1, length(labels))
    people = if !isempty(geo)
        geo.pop[geo.fips .== locale]
     else # this will off by a tiny bit because of rounding
        getproperty(cumhist, :unexposed_total)[1] + getproperty(cumhist, :infectious_total)[1]
     end   
     co_pal = length(plotcols) == 2 ? [theme_palette(thm)[2], theme_palette(thm)[4]] : theme_palette(thm)

     xtickrange = caldays[10]:Day(180):caldays[length(caldays)-10]
     xtickfmt = Dates.format.(xtickrange, "yyyy-mm-dd")
 

    # the plot
    plot(       caldays, newseries[days, :], 
                size = (700,500),
                label = labels, 
                lw=1.5,
                # bar_width=1,
                title = "Daily Change for $people people over $n days",
                xlabel = "Simulation Days",
                xticks = (xtickrange, xtickfmt),
                # xticks = caldays[10]:Day(180):caldays[length(caldays)-10],
                yaxis = ("People"),
                color_palette = co_pal,
                reuse =false,
                background_color_legend=nothing,
                foreground_color_legend=nothing
             )
end


function selpos!(tab::Table)
    for c in columns(tab)
        for i in eachindex(c)
            c[i] = c[i] > 0.0 ? c[i] : 0.0
        end
    end
end


function daily_cases_plot(series, popdat, locale, plotcols=[:total]; days="all", geo=[], thm=:ggplot2)
    cumhist = series[locale].cum
    newhist = series[locale].new
    locdat = popdat[locale]

    theme(thm, foreground_color_border=:black, 
        tickfontsize=9, gridlinewidth=1)    
    co_pal = length(plotcols) == 2 ? [theme_palette(thm)[2], theme_palette(thm)[4]] : theme_palette(thm)


    n = length(newhist)
    days = days == "all" ? (1:n) : days
    numcols = length(AGEGRPS) + 1
    numseries = length(plotcols)
    caldays = newhist.calday[days]

    # construct the data series
    maxtimes = maximum(length.(locdat.sickday)) - 1  # maximum no. of times anyone has gotten infected
    daily_cases_series = zeros(Int, n, numcols)
    for i in 2:maxtimes+1
        for age in AGEGRPS  # accumulate all days on which anyone got sick the 1st, 2nd, 3rd... time
            daygotsick = countmap(get.(locdat.sickday[locdat.agegrp .== age],i,0))
            for (k,v) in daygotsick
                if k == 0
                    continue  # ignore people who never got infected
                end
                daily_cases_series[k, Int(age)] += v
            end
        end
    end

    daily_cases_series[:, numcols] = sum(daily_cases_series, dims=2)  # sum the columns across each row

    # prepare plot series
    series_selector = Int[]
    labels = String[]
    for col in plotcols
        if col in AGEGRPS
            push!(series_selector, Int(col))
            push!(labels, titlecase(string(col)))
        elseif col === :total
            push!(series_selector, numcols)
            push!(labels, "Total")
        else
            throw(DomainError(col, "plotcols argument must be array of agegrp enums and/or :total"))
        end
    end

    # annotations and labels
    labels = reshape([labels...], 1, length(labels))
    people = if !isempty(geo)
                    geo.pop[geo.fips .== locale]
                else # this will off by a tiny bit because of rounding
                    getproperty(cumhist, :unexposed_total)[1] + getproperty(cumhist, :infectious_total)[1]
                end   
    died =  getproperty(cumhist, :dead_total)[end]   
    unexp = getproperty(cumhist, :unexposed_total)[end]
    infected = people - unexp    
    recovered = infected - died  

    xtickrange = caldays[10]:Day(180):caldays[length(caldays)-10]
    xtickfmt = Dates.format.(xtickrange, "yyyy-mm-dd")



    # the plot
    plot(   caldays, daily_cases_series[days, series_selector], 
            label=labels,
            size = (700,500),
            lw = 1.5,
            title = "Daily Cases",
            xlabel = "Simulation Days",
            xticks = (xtickrange, xtickfmt),
            # xticks = caldays[10]:Day(180):caldays[length(caldays)-10],
            # yaxis = ("People"),
            ylabel = "People",
            color_palette = co_pal,
            reuse =false,
            # background_color_legend=nothing,
            # foreground_color_legend=nothing
         )

    annotate!(caldays[1] + Day(6), 0.51 * ylims()[2],              # half_yscale,
         text("Died: $died\nInfected: $infected\nRecovered: $recovered\nUnexposed: $unexp", 
             11, :left))
           
end