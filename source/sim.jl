
####################################################################################
#   simulation runner: ILM Model
####################################################################################


function buildsim(ndays, locales;
    dovax = false,
    paramdir = "../parameters",
    geofilename = "../data/geo2data.csv", 
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    variantfilename = "variants.yml")

    locales = locales isa Int ? [locales] : locales

    model = setup(ndays, locales; 
        dovax=dovax, 
        paramdir=paramdir,
        geofilename=geofilename, 
        socialfilename=socialfilename,
        vaccinefilename=vaccinefilename,
        variantfilename=variantfilename
        )

    return model
end


function runsim(model; 
            runcases=[], 
            showr0 = false, 
            silent=true, 
            dovariant=false,
            dovax=false
            )

    empty_all_caches!() # from previous runs

    # split up some members of model
    ndays = model.ndays
    locales = model.locales
    transitionset = model.transitionset  # transition arrays
    trvec = model.trvec # preallocated small vector
    popdat = deepcopy(model.dat["popdat"])   # Copy the population data so model can be reused!!!
    agegrp_idx = model.dat["agegrp_idx"]   # first key is locale
    cumhistmx = deepcopy(model.dat["cumhistmx"])   # first key is locale
    newhistmx = deepcopy(model.dat["newhistmx"])   # first key is locale
    geodf = model.geo
    infectset = model.infectset
    socialparams = model.social
    vaxset = model.vaxset
    vaxschedset = model.vaxschedset

    # restart the day counter to zero
    reset!(day_ctr, :day)  # return and reset key to 0 :day leftover from prior runs

    locales = locales   # force local scope to be visible in the loop

    sdcases = Dict{Symbol, Spreadcase}()  # hold definitions of spreadcases


    # execution timing
    vaxtime = 0
    sprtime = 0
    trtime = 0
    idxtime = 0
    histtime = 0


    ######################
    # simulation loop
    ######################
    for i = 1:ndays
        inc!(day_ctr, :day)  # increment the simulation day counter
        silent || println("simulation day: ", day_ctr[:day])

        for loc in locales     

            silent || println("Simulation starting for location $loc")
            
            locdat = popdat[loc]
            age_idx_loc = agegrp_idx[loc]  # indices by agegrp
            
            density_factor = geodf[geodf[!, :fips] .== loc, :density_factor][]  # TODO not a good place for this
            
            for case in runcases
                case(loc, popdat, socialparams, infectset, sdcases, age_idx_loc; startofday=true)  # TODO extend ages to be any filter for 
            end                                                 # who participates in a given case

            # filter for key people
            idxtime += @elapsed begin
                # @bp
                infect_idx = findall(locdat.status .== infectious)
                contactable_idx = findall(locdat.status .!= dead)
            end

            # if dovax give shots
            dovax && (vaxtime += @elapsed vaccinate!(locdat, vaxschedset, contactable_idx, vaxset))
            # @show size(contactable_idx), vaxset

            # two fundamental steps of the simulation: spread! and transition!
            sprtime += @elapsed spread!(locdat, infect_idx, contactable_idx, sdcases, socialparams, 
                                        infectset, vaxset, density_factor, dovax, dovariant)   
            trtime += @elapsed transition!(locdat, infect_idx, infectset, transitionset, vaxset, dovax, dovariant; trvec=trvec)                        

            for case in runcases
                case(loc, popdat, socialparams, infectset, sdcases, age_idx_loc; startofday=false)  # TODO extend ages to be any filter for 
            end                                                 # who participates in a given case

            # r0 displayed every 10 days
            if showr0 && (mod(day_ctr[:day],10) == 0)   # do we ever want to do this by locale -- maybe
                current_r0 = r0_sim(locdat, age_dist=age_dist, dectree=dectree, socialparams=socialparams, infectparams=infectparams, sdcases=sdcases)
                println("day $(day_ctr[:day]), locale $loc: rt = $current_r0")
            end

        end

        # accumulate simulation statistics
        histtime += @elapsed do_history!(locales, popdat, cumhistmx, newhistmx, agegrp_idx)

        silent || println("Simulation completed for $(day_ctr[:day]) days.")
    end

    #######################

    # simulation history series for plotting: arrays NOT dataframes
    series = Dict(loc=>Dict(:cum=>cumhistmx[loc], :new=>newhistmx[loc]) for loc in locales)

    # sum agegrps to total for all series groups (by agegrp)
    hist_total_agegrps!(series, locales)

    for loc in locales
        add_totinfected_series!(series, loc)
        # we might need a totvaccinated series
    end

    @show idxtime, vaxtime, sprtime, trtime, histtime

    return popdat, series
end



################################################################################
#  Update daily history series
################################################################################

@views function do_history!(locales, popdat, cumhist, newhist, agegrp_idx)
    thisday = day_ctr[:day]
    for loc in locales
        dat = popdat[loc]  # source
        cumdat = cumhist[loc]   # sink
        newdat = newhist[loc]   # sink

        @inbounds for age in instances(agegrp)
            int_age = Int(age)

            # get the source data: status
            dat_age = dat[agegrp_idx[loc][age]]
            status_today = countmap(dat_age.status)    # cumulative position for thisday, keys are Enum status

            # get the source data: conditions in (nil, mild, sick, severe)
            filt_infectious = findall(dat_age.status .== infectious)
            if length(filt_infectious) > 0
                sick_today = countmap(dat_age.cond[filt_infectious])  # keys are enum condition
            else   # there can be days when no one is infected
                sick_today = Dict()
            end

            # get the source data: vaccination
            filt_vaccinated = findall(last.(dat_age.vaxrcvd) .!= :none)
            if length(filt_vaccinated) > 0
                vax_today = countmap(last.(dat_age.vaxrcvd[filt_vaccinated]))  # keys are symbol
            else
                vax_today = Dict()
            end


            #
            # cumulative data
            #
            for st in statuses 
                cumdat[thisday, getindex(map2series, Symbol(st))[int_age]] = get(status_today, st, 0)
            end

            # insert condition into sink: cum
            for co in infectious_cases # 5:8
                cumdat[thisday, getindex(map2series, Symbol(co))[int_age]] = get(sick_today, co, 0)
            end

            # insert vaccination into sink: cum
            for vax in vaxlist
                cumdat[thisday, getindex(map2series, vax)[int_age]] = get(vax_today, vax, 0)
            end


            #
            # new (each day) data
            #
            for ac in allconds   # status and conds
                if thisday == 1
                    newdat[thisday, getindex(map2series, Symbol(ac))[int_age]] = get(status_today, Int(ac), 0)
                else  # on all other days
                    newdat[thisday, getindex(map2series, Symbol(ac))[int_age]] = (
                        cumdat[thisday, getindex(map2series,Symbol(ac))[int_age]] 
                        - cumdat[thisday - 1, getindex(map2series, Symbol(ac))[int_age]]
                        )         
                end
            end

            for vax in vaxlist
                if thisday == 1
                    newdat[thisday, getindex(map2series, vax)[int_age]] = get(vax_today, vax, 0)
                else  # on all other days
                    newdat[thisday, getindex(map2series, vax)[int_age]] = (
                        cumdat[thisday, getindex(map2series, vax)[int_age]] 
                        - cumdat[thisday - 1, getindex(map2series, vax)[int_age]]
                        )         
                end
            end                

        end # for age in agegrps

    end # for loc in locales

end # function


function hist_total_agegrps!(series, locales)
    for loc in locales
        for kind in [:cum, :new]
            for cond in Int.(allconds)  # infectious cases and statuses
                series[loc][kind][:,map2series[cond][totalcol]] = sum(series[loc][kind][:,map2series[cond][collect(Int.(agegrps))]],dims=2)
            end
            for vax in vaxlist
                series[loc][kind][:,map2series[vax][totalcol]] = sum(series[loc][kind][:,map2series[vax][collect(Int.(agegrps))]],dims=2)
            end
        end
    end
end


function review_history(histmx)
    for i in 1:size(histmx, 3)
        println("   *** Day $i ***")
        display(hcat(histmx[:,:,i], [:Unexposed, :Infectious, :Recovered, :Dead, :Nil, :Mild, :Sick, :Severe]))
        print("Press enter or type q and enter to quit> "); resp = chomp(readline())
        if resp == "q"; break; end
    end
end


# a single locale that already has both new and cum series
function add_totinfected_series!(series, locale)
    if !(haskey(series[locale], :cum) && haskey(series[locale], :new))
        error("locale series must contain both :cum and :new series")
        return
    end
    # for new
    @views begin
        n = size(series[locale][:new],1)
        series[locale][:new] = hcat(series[locale][:new], zeros(Int, n, 6))
        series[locale][:new][:,map2series.totinfected] = ( (series[locale][:new][:,map2series.unexposed] .< 0 ) .*
                                                          abs.(series[locale][:new][:,map2series.unexposed]) ) 
        # for cum
        series[locale][:cum] = hcat(series[locale][:cum], zeros(Int, n, 6))
        cumsum!(series[locale][:cum][:,map2series.totinfected], series[locale][:new][:,map2series.totinfected], dims=1)  
    end
    return
end



#####################################################################################
#  other functions used in simulation
#####################################################################################


function cleanup_stash(stash)
    for k in keys(stash)
        delete!(stash, k)
    end
end


function empty_all_caches!()
    # empty tracking queues
    !isempty(spreadq) && (deleteat!(spreadq, 1:length(spreadq)))   
    !isempty(transq) && (deleteat!(transq, 1:length(transq)))   
    !isempty(tntq) && (deleteat!(tntq, 1:length(tntq)))   
    !isempty(r0q) && (deleteat!(r0q, 1:length(r0q)))  
end


"""
    optfindall(p, X, maxlen=0)

Returns indices to X where p, a filter, is true.
Filters should be anonymous functions.
For maxlen=0, the length of the temporary vector is length(x).
For maxlen=n, the length of the temporary vector is n.
For maxlen=0.x, the length of temporary vector is 0.x * length(x) and
x should be in (0.0, 1.0).
"""
function optfindall(p, X, maxlen=1)
    if maxlen==1
        out = Vector{Int}(undef, length(X))
    elseif isa(maxlen, Int)
        out = Vector{Int}(undef, maxlen)
    else
        out = Vector{Int}(undef, floor(Int, maxlen * length(X)))
    end
    ind = 0
    @inbounds for (i, x) in pairs(X)
        if p(x)
            out[ind+=1] = i
        end
    end
    resize!(out, ind)
    return out
end


#######################################################################################
#  probability
#######################################################################################


# discrete integer histogram
function bucket(x; vals)
    if isempty(vals)
        vals = range(minimum(x), stop = maximum(x))
    end
    [count(x .== i) for i in vals]
end


# range counts to discretize PDF of continuous outcomes
function histo(x)
    big = ceil(maximum(x))
    bins = Int(big)
    sm = floor(minimum(x))
    ret = zeros(Int, bins)
    binbounds = collect(1:bins)
    @inbounds for i = 1:bins
        n = count(x -> i-1 < x <= i,x)
        ret[i] = Int(n)
    end
    return ret, binbounds
end


"""
Returns continuous value that represents gamma outcome for a given
approximate center point (scale value of gamma).  We can interpret this
as a funny sort of probability or as a number outcome from a gamma
distributed sample.
1.2 provides a good shape with long tail right and big clump left
"""
function gamma_prob(target; shape=1.0)
    @assert 0.0 <= target <= 99.0 "target must be between 0.0 and 99.0"
    dgamma = Gamma(shape,target)
    pr = rand(dgamma, 1)[1] / 100.0
end


"""
Returns a single number of successes for a
sampled outcome of cnt tries with the input pr of success.
"""
function binomial_one_sample(cnt, pr)::Int
    return rand.(Binomial.(cnt, pr))
end



"""
    categorical_sim(prs::Vector{Float64})
    categorical_sim(prs::Vector{Float64}, n::Int, do_assert=true)

Approximates sampling from a categorical distribution.
prs is an array of floats that must sum to 1.0.
do_assert determines if an assert tests this sum. For a single trial, this runs
in 10% of the time of rand(Categorical(prs)). For multiple trials, the 
second method runs in less than 50% of the time.

The second method generates results for n trials. 
The assert test is done only once if do_assert is true.
"""
function categorical_sim(prs)
    if !isapprox(sum(prs), 1.0)
        return 0
    end
    x = rand()
    cumpr = 0.0
    i = 0
    for pr in prs
        cumpr += pr
        i += 1
        if x <= cumpr 
            break
        end
    end
    i
end

function categorical_sim(prs, n::Int, do_assert=true)
    do_assert && @assert isapprox(sum(prs), 1.0)
    ret = Vector{Int}(undef, n)
    
    @inbounds for i in 1:n
        ret[i] = categorical_sim(prs, false)
    end
    ret
end
