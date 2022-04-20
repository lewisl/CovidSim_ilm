#########################################################################################
# tracking.jl
#########################################################################################

# gr()   # initialize plotting backend for Plots


# for debugging simulations: daily outcome entries as named tuples
const spreadq = []
const transq = []
const tntq = []
const r0q = []


# for Johns Hopkins US actual data
struct Col_ref
    date::String
    col::Int64
end


# tracking statistics

function showq(qname)
    for item in qname
        println(item)
    end
end


function reviewdays(q=spreadq)
    for it in q
        println(it)
        print("\nPress enter to continue, q enter to quit.> ");
        ans = chomp(readline())
        if ans == "q"
            break
        end
    end
end

function reviewdays(df::DataFrame)
    for it in eachrow(df)
        display(it)
        print("\nPress enter to continue, q enter to quit.> ");
        ans = chomp(readline())
        if ans == "q"
            break
        end
    end
end


#################################################################################
#  Epidemiological Stats -- very preliminary
#################################################################################

# outcomes per agegrp
function virus_outcome(series, locale; agegrp=totalcol, base=:infected)  # denom in (:infected, :pop, :none)
    map2series = series[locale].cols

    n = size(series[locale].cum, 1)
    outcomes = Dict{Symbol, Float64}()  # TODO should we have integer outcomes for totals when base=:none?
    agegrp = Int(agegrp)
    
    # each denominator for data summary
    total_pop = series[locale].cum[1, map2series[:unexposed][agegrp]] + series[locale].cum[1, map2series[:infectious][agegrp]]
    total_infected = series[locale].cum[end, map2series[:totinfected][agegrp]]

    denom = if base == :pop 
                total_pop 
            elseif base == :infected
                total_infected
            else  # :none or wrong entry
                1
            end

    for cond in statuses
        ssym = Symbol(cond)
        outcomes[ssym] = series[locale].cum[n, map2series[ssym][agegrp]] / denom
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
    days="all", geo=[], thm=:wong2)

    cumhist = series[locale].cum
    newhist = series[locale].new

    # theme(:ggplot2, foreground_color_border =:black, reuse = false)
    theme(thm, foreground_color_border=:black, 
          tickfontsize=9, gridlinewidth=1)

    !(typeof(plotcols) <: Array) && (plotcols = [plotcols])

    # the data is the 2d array cumseries
    n = size(cumhist, 1)
    days = days == "all" ? (1:n) : days
    cumseries = hcat(columns(getproperties(cumhist,Tuple(Symbol(plcol,"_","total") for plcol in plotcols)))...)

    # labels and annotations
    labels = [titlecase(string(col)) for col in plotcols]
    labels = reshape([labels...], 1, length(labels))
    people = if !isempty(geo)
                geo[geo[:,fips] .== locale, popsize][1]
             else # this will off by a tiny bit because of rounding
                getproperty(cumhist, :unexposed_total)[1] + getproperty(cumhist, :infectious_total)[1]
             end   
    cityname = !isempty(geo) ? geo[geo[:,fips] .== locale, city][1] : ""
    died =  getproperty(cumhist, :dead_total)[end]   #      series.data[locale].cum[end, map2series.dead[totalcol]]
    unexp = getproperty(cumhist, :unexposed_total)[end]
    infected = people - unexp    #series.data[locale].cum[end, map2series.unexposed[totalcol]]
    recovered = infected - died  # series.data[locale].cum[end, map2series.recovered[totalcol]]
    co_pal = length(plotcols) == 2 ? [theme_palette(thm)[2], theme_palette(thm)[4]] : theme_palette(thm)
 

    # the plot
    plot(   days, cumseries[1:end,1:end], 
            size = (700,500),
            label = labels, 
            lw=2.3,
            title = "Covid for $people people in $cityname over $n days\nActive Cases for Each Day",
            xlabel = "Simulation Days",
            ylabel = "People",
            legendfontsize = 10,
            color_palette = co_pal,
            reuse = false,
            legend_position = :right
        )
    plot!(annotate = ((6, 0.51 * ylims()[2],              # half_yscale,
            Plots.text("Died: $died\nInfected: $infected\nRecovered: $recovered\nUnexposed: $unexp", 
                11, :left))))
end


function newplot(series, locale, plcols=[:infectious]; days="all")

    # pyplot()
    theme(:ggplot2, foreground_color_border =:black)

    !(typeof(plcols) <: Array) && (plcols = [plcols])

    map2series = series.mapper

    # the data and labels
    n = size(series.data[locale].new, 1)
    days = days == "all" ? (1:n) : days
    newseries = series.data[locale].new[days, [getproperty(map2series, i)[totalcol] for i in plcols]]
    labels = [titlecase(string(col)) for col in plcols]
    labels = reshape([labels...], 1, length(labels))
    people = series.data[locale].cum[1, map2series.unexposed[totalcol]] + series.data[locale].cum[1, map2series.infectious[totalcol]]

    # the plot
    groupedbar( days, newseries[1:end, 1:end], 
                size = (700,500),
                label = labels, 
                lw=0.2,
                bar_width=1,
                title = "Covid Daily Change for $people people over $n days",
                xlabel = "Simulation Days",
                yaxis = ("People"),
                reuse =false
        )
    # gui()
end


function day2df(spreadq::Array)
    spreadseries = DataFrame(spreadq)

    spreadseries[!, :cuminfected] .= zeros(Int, size(spreadseries,1))
    spreadseries[1, :cuminfected] = copy(spreadseries[1,:infected])
    for i = 2:size(spreadseries,1)
       spreadseries[i,:cuminfected] = spreadseries[i-1,:cuminfected] + spreadseries[i,:infected]
    end

    return spreadseries
end


function dayplot(spreadq, plseries=[])
    dayplot(DataFrame(spreadq), plseries)
end


function dayplot(spreadseries::DataFrame, plseries=[])
    
    theme(:ggplot2, foreground_color_border =:black)
    
    pl = bar(   spreadseries[!,:day], spreadseries[!,:infected],label="Infected", 
            lw=0.2,
            bar_width=1,
            size = (700,300),
            dpi=180,
            xlabel="Simulation Days", 
            ylabel="People", 
            title="Daily Spread of Covid",
            bg_legend=:white)
    
    for addlseries in plseries
        lbl = titlecase(string(addlseries))
        plot!(spreadseries[!,:day], spreadseries[!,addlseries],label=lbl, lw=2)
    end
    # gui()  # force instant plot window
    return pl
end


function day_animate2(spreadseries)
    n = size(spreadseries,1)
    # daymat = Matrix(spreadseries)

    xd = spreadseries[1:5,:]

    topy = max(maximum(spreadseries[!,:spreaders]),maximum(spreadseries[!,:contacts]),
                maximum(spreadseries[!,:touched]),maximum(spreadseries[!,:infected]) )

    @df xd plot(:day, [:spreaders :contacts :touched :infected], color=^([:red :blue :green :orange]),
                labels=^(["Spreaders" "Contacts" "Touched" "Infected"]),dpi=200, lw=2,ylim=(0,topy))

    for i = 5:2:n
        xd = spreadseries[i-2:i,:]

        @df xd plot!(:day, [:spreaders :contacts :touched :infected], color=^([:red :blue :green :orange]),
                 labels=false, dpi=200, lw=2, ylim=(0,3e4))
        gui()

        if i < round(Int, n/4)
            sleep(0.3)
        elseif i < round(Int,n/2)
            sleep(0.1)
        else
            sleep(.001)
        end
        # print("\nPress enter to continue, q enter to quit.> ");
        # ans = chomp(readline()) 
        # if ans == "q"
        #     break
        # end    
    end
end

function catplot()
    # groupedbar(datpct', bar_position=:stack,label=labels)
    # plot!(xticks=(1:3,["one", "two", "three"]))
    #= julia> datpct'
                3×5 LinearAlgebra.Adjoint{Float64,Array{Float64,2}}:
                 0.32036   0.30348    0.141939  0.155273  0.0789479
                 0.199046  0.413894   0.201274  0.054734  0.131053
                 0.252381  0.0677606  0.26775   0.161113  0.250996
    =#
end

# Plots.AnimatedGif("/var/folders/mf/73qj_8c91dzg4sw459_7mchm0000gn/T/jl_Js4px6.gif")