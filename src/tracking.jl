#########################################################################################
# tracking.jl
#########################################################################################

plotly()

# for Johns Hopkins US actual data
struct Col_ref
    date::String
    col::Int64
end



#################################################################################
#  Simulation Stats -- very preliminary
#################################################################################
"""
    function stat1(series, locale)

Returns a table of status outcomes for the simulation by age group.
"""
function stat1(series, locale)
    cumhist = series[locale].cum
    newhist = series[locale].new

    stat1 = Table(item=Symbol[], total=Int[], age0_19=Int[], age20_39=Int[], age40_59=Int[], age60_79=Int[], age80_up=Int[])
    calc_cols = [:total, :age0_19, :age20_39, :age40_59, :age60_79, :age80_up]

    stat1 = push!(stat1, (item=:pop,
                    total    = getproperty(cumhist, :unexposed_total)[1] + getproperty(cumhist, :infectious_total)[1],
                    age0_19  = getproperty(cumhist, :unexposed_age0_19)[1] + getproperty(cumhist, :infectious_age0_19)[1],
                    age20_39 = getproperty(cumhist, :unexposed_age20_39)[1] + getproperty(cumhist, :infectious_age20_39)[1],
                    age40_59 = getproperty(cumhist, :unexposed_age40_59)[1] + getproperty(cumhist, :infectious_age40_59)[1],
                    age60_79 = getproperty(cumhist, :unexposed_age60_79)[1] + getproperty(cumhist, :infectious_age60_79)[1],
                    age80_up = getproperty(cumhist, :unexposed_age80_up)[1] + getproperty(cumhist, :infectious_age80_up)[1])
                )
                
    push!(stat1, (item=:died,
                  total    = getproperty(cumhist, :dead_total)[end],
                  age0_19  = getproperty(cumhist, :dead_age0_19)[end],
                  age20_39 = getproperty(cumhist, :dead_age20_39)[end],
                  age40_59 = getproperty(cumhist, :dead_age40_59)[end],
                  age60_79 = getproperty(cumhist, :dead_age60_79)[end],
                  age80_up = getproperty(cumhist, :dead_age80_up)[end])
                  )

    push!(stat1, (item=:unexposed,
                  total    = getproperty(cumhist, :unexposed_total)[end],
                  age0_19  = getproperty(cumhist, :unexposed_age0_19)[end],
                  age20_39 = getproperty(cumhist, :unexposed_age20_39)[end],
                  age40_59 = getproperty(cumhist, :unexposed_age40_59)[end],
                  age60_79 = getproperty(cumhist, :unexposed_age60_79)[end],
                  age80_up = getproperty(cumhist, :unexposed_age80_up)[end])
                  )

    # inter-row calculation infected = pop - unexposed
    poprow = stat1[stat1.item .=== :pop]
    unexprow = stat1[stat1.item .=== :unexposed]
    newrow = zeros(Int, 6)
    for (i, col) in enumerate(calc_cols)
        newrow[i] = getproperty(poprow, col)[] - getproperty(unexprow, col)[]
    end

    push!(stat1, (item = :ever_infected,
        total = newrow[1],
        age0_19 = newrow[2],
        age20_39 = newrow[3],
        age40_59 = newrow[4],
        age60_79 = newrow[5],
        age80_up = newrow[6])
        )

    # inter-row calculation recovered = infected - died
    diedrow = stat1[stat1.item .=== :died]
    infectedrow = stat1[stat1.item .=== :ever_infected]
    newrow[:] .= 0
    for (i, col) in enumerate(calc_cols)
        newrow[i] = (getproperty(infectedrow, col)[] - getproperty(diedrow, col)[] 
                    - getproperty(cumhist, Symbol(:infectious, "_", col))[end])
    end

    push!(stat1, (item=:recovered,
                  total    = newrow[1],
                  age0_19  = newrow[2],
                  age20_39 = newrow[3],
                  age40_59 = newrow[4],
                  age60_79 = newrow[5],
                  age80_up = newrow[6])
                  )

    return stat1

end


# outcomes per agegrp  THIS IS REALLY JUST THE PCT CALCULATION:  PUT IT IN STAT1
function virus_outcome(series, locale; agegrp=totalcol, base=:infected)  # denom in (:infected, :pop, :none)

    n = size(series[locale].cum, 1)
    outcomes = Dict{Symbol, Float64}()  # TODO should we have integer outcomes for totals when base=:none?
    # agegrp = Int(agegrp)
    if agegrp == totalcol
        agegrp = :total
    end
    
    # each denominator for data summary
    # total_pop = series[locale].cum[1, map2series[:unexposed][agegrp]] + series[locale].cum[1, map2series[:infectious][agegrp]]
    total_pop = getproperty(series[locale].cum, Symbol(unexposed, "_", agegrp))[1] + getproperty(series[locale].cum, Symbol(infectious, "_", agegrp))[1]
    # total_infected = series[locale].cum[end, map2series[:totinfected][agegrp]]
    total_infected = getproperty(series[locale].cum, Symbol(:totinfected, "_", agegrp))[end]

    denom = if base == :pop 
                total_pop 
            elseif base == :infected
                total_infected
            else  # :none or wrong entry
                1
            end

    for cond in statuses
        stsym = Symbol(cond)
        outcomes[stsym] = getproperty(series[locale].cum, Symbol(stsym, "_", agegrp))[n] / denom
    end

    return outcomes
end


function onecond(series, locale, cond; case=:new, agegrp=totalcol, filt=:pos)
    map2series = series[locale].cols

    datacol = getproperty(series[locale], case)[:,map2series[Symbol(cond)][Int(agegrp)]]
    if filt == :pos
        datacol[datacol .> 0]
    else
        datacol
    end
end

###########################################################################################
#  Plotting
###########################################################################################


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
    caldays = cumhist.calday[days]
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
    died =  getproperty(cumhist, :dead_total)[end]   #      series.data[locale].cum[end, map2series.dead[totalcol]]
    unexp = getproperty(cumhist, :unexposed_total)[end]
    infected = people - unexp    #series.data[locale].cum[end, map2series.unexposed[totalcol]]
    recovered = infected - died  # series.data[locale].cum[end, map2series.recovered[totalcol]]
    co_pal = length(plotcols) == 2 ? [theme_palette(thm)[2], theme_palette(thm)[4]] : theme_palette(thm)
 

    # the plot
    plot(   caldays, cumseries[days,1:end], 
            size = (700,500),
            label = labels, 
            lw=2.3,
            title = "Covid for $people people for $n days\nActive Cases for Each Day",
            xlabel = "Simulation Days",
            xticks = caldays[10]:Day(180):caldays[length(caldays)-10],
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


    # the plot
    bar(        caldays, newseries[days, 1:end], 
                size = (700,500),
                label = labels, 
                lw=0.2,
                bar_width=1,
                title = "Daily Change for $people people over $n days",
                xlabel = "Simulation Days",
                xticks = caldays[10]:Day(180):caldays[length(caldays)-10],
                yaxis = ("People"),
                color_palette = co_pal,
                reuse =false,
                background_color_legend=nothing,
                foreground_color_legend=nothing
             )
end

