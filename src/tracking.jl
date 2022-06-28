#########################################################################################
# tracking.jl
#########################################################################################

plotly()  # Plots backend for plotting

# for Johns Hopkins US actual data
struct Col_ref
    date::String
    col::Int64
end



#################################################################################
#  Simulation Stats -- very preliminary
#################################################################################

gt(dict, key) = get(dict, key, 0)


"""
    function stat1(series, locale)

Returns a table of status outcomes for the simulation by age group.
"""
function stat_cond(series, locale)
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
            total = newrow[1], age0_19 = newrow[2], age20_39 = newrow[3], 
            age40_59 = newrow[4], age60_79 = newrow[5], age80_up = newrow[6])
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
                    total    = newrow[1], age0_19  = newrow[2], age20_39 = newrow[3],
                    age40_59 = newrow[4], age60_79 = newrow[5], age80_up = newrow[6])
                  )


    @show sum(stat1.total)

    # add pct of population columns
    statpct = Table(item = stat1.item, total = stat1.total, total_pct = stat1.total ./ stat1.total[1],   # 
            age0_19 = stat1.age0_19, age0_19_pct = stat1.age0_19 ./ stat1.age0_19[1], 
            age20_39=stat1.age20_39, age20_39_pct = stat1.age20_39 ./ stat1.age20_39[1],
            age40_59 = stat1.age40_59, age40_59_pct = stat1.age40_59 ./ stat1.age40_59[1], 
            age60_79=stat1.age60_79, age60_79_pct = stat1.age60_79 ./ stat1.age60_79[1], 
            age80_up=stat1.age80_up, age80_up_pct = stat1.age80_up ./ stat1.age80_up[1])

    return statpct

end


function stat_vax(popdat, locale)
    td = @Select(vaxrcvd, agegrp, status)(popdat[locale]) # view with 2 columns on popdat
    vaxes = td.vaxrcvd

    stat_vax = Table(item=Symbol[], total=Int[], age0_19=Int[], age20_39=Int[], age40_59=Int[], age60_79=Int[], age80_up=Int[])
    calc_cols = [:total, :age0_19, :age20_39, :age40_59, :age60_79, :age80_up]
    
    col_age0_19  = countmap(last.(vaxes[(td.agegrp .== age0_19) .& (td.status .!= dead)]))
    col_age20_39 = countmap(last.(vaxes[(td.agegrp .== age20_39) .& (td.status .!= dead)]))
    col_age40_59 = countmap(last.(vaxes[(td.agegrp .== age40_59) .& (td.status .!= dead)]))
    col_age60_79 = countmap(last.(vaxes[(td.agegrp .== age40_59) .& (td.status .!= dead)]))
    col_age80_up = countmap(last.(vaxes[(td.agegrp .== age80_up) .& (td.status .!= dead)]))
    col_total    = countmap(last.(vaxes))
    
    for k in (:none, :Moderna, :Pfizer, :JnJ)
        push!(stat_vax, 
              (item=k, total=gt(col_total,k), age0_19=gt(col_age0_19,k), age20_39=gt(col_age20_39,k), age40_59=gt(col_age40_59,k),
               age60_79=gt(col_age60_79,k), age80_up=gt(col_age80_up,k))
               )
    end

    stat_pct = @Select(item, total, total_pct = $total ./ sum($total),
                age0_19,  age0_19_pct  = $age0_19  ./ sum($age0_19),
                age20_39, age20_39_pct = $age20_39 ./ sum($age20_39),
                age40_59, age40_59_pct = $age40_59 ./ sum($age40_59),
                age60_79, age60_79_pct = $age60_79 ./ sum($age60_79),
                age80_up, age80_up_pct = $age80_up ./ sum($age80_up))(stat_vax)

    return stat_pct

end


function stat_repeat(popdat, locale)
    thisdat = @Select(variant, agegrp)(popdat[locale]) # view with 2 columns on popdat
    variants = thisdat.variant

    stat_count = Table(item=Symbol[], total=Int[], age0_19=Int[], age20_39=Int[], age40_59=Int[], age60_79=Int[], age80_up=Int[])

    col_total    = sort(countmap(length.(variants)))
    col_age0_19  = sort(countmap(length.(variants[thisdat.agegrp .== age0_19])))
    col_age20_39 = sort(countmap(length.(variants[thisdat.agegrp .== age20_39])))
    col_age40_59 = sort(countmap(length.(variants[thisdat.agegrp .== age40_59])))
    col_age60_79 = sort(countmap(length.(variants[thisdat.agegrp .== age40_59])))
    col_age80_up = sort(countmap(length.(variants[thisdat.agegrp .== age80_up])))

    for k in 0:10
        push!(stat_count, 
              (item=Symbol("Times_",k), total=gt(col_total,k), age0_19=gt(col_age0_19,k), age20_39=gt(col_age20_39,k), age40_59=gt(col_age40_59,k),
               age60_79=gt(col_age60_79,k), age80_up=gt(col_age80_up,k))
               )
    end

    stat_pct = @Select(item, total, total_pct = $total ./ sum($total),
                        age0_19,  age0_19_pct  = $age0_19  ./ sum($age0_19),
                        age20_39, age20_39_pct = $age20_39 ./ sum($age20_39),
                        age40_59, age40_59_pct = $age40_59 ./ sum($age40_59),
                        age60_79, age60_79_pct = $age60_79 ./ sum($age60_79),
                        age80_up, age80_up_pct = $age80_up ./ sum($age80_up))(stat_count)

    return stat_pct
end

function stat_breakout(popdat, locale)

    @Select(vaxday, sickday, agegrp)(locdat[((locdat.status .== infectious) .| (locdat.status .== recovered)) .& 
    (locdat.vaxstatus .!= :none) .& (last.(locdat.sickday) .> first.(locdat.vaxday))])

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

