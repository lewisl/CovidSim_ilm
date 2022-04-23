
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

        #=
        model = (ndays=ndays, locales=locales, dat=datadict, series=series, geo=geodata, 
                transitionset=transitionset, vaxset=vaxset, vaxschedset=vaxschedset, infectset=infectset, 
                social=socialparams, trvec=trvec)  
        =#    

    return model
end


function runsim(model; 
            runcases=[], 
            showr0 = false, 
            silent=true, 
            dovax=false
            )

    empty_all_caches!() # from previous runs

    # split up  members of model and initialize
        ndays = model.ndays
        locales = model.locales
        transitionset = model.transitionset  # transition arrays
        trvec = model.trvec # preallocated small vector
        popdat = deepcopy(model.dat["popdat"])   # Copy the population data so model can be reused!!!
        agegrp_idx = model.dat["agegrp_idx"]   # first key is locale
        series = deepcopy(model.series)  # contains series.mapper and series.data, which is a dict of locales, each local includes .cum and .new
        geodf = model.geo
        infectset = model.infectset
        socialparams = model.social
        vaxset = model.vaxset
        vaxschedset = model.vaxschedset
        for sched in values(vaxschedset) 
            for vax in values(sched.vaxesincluded) 
                vax.doses = vax.starting_doses   
            end
        end
        contact_vector = zeros(Int, 10)    

    # restart the day counter to zero
    reset!(day_ctr, :day)  # return and reset key to 0 :day leftover from prior runs

    sdcases = Dict{Symbol, Spreadcase}()  # hold definitions of spreadcases


    # execution timers
    vaxtime = 0     # vaccinate
    sprtime = 0     # spread infection
    trtime = 0      # transition infected population through stages of illness
    idxtime = 0     # calculate indices for infectious and susceptible
    histtime = 0    # update history time series
    misc_time = 0
    totalsimtime = 0


    ######################
    # simulation loop
    ######################
    totalsimtime += @elapsed begin

    for loc in locales     

        silent || println("Simulation starting for location $loc")
        
        # this should be the first and only place to deref the locale (as loc)
        locdat = popdat[loc]  
        newhist = series[loc].new
        cumhist = series[loc].cum
        age_idx_loc = agegrp_idx[loc]  # indices by agegrp
        density_factor = geodf[geodf[!, :fips] .== loc, :density_factor][]  # for the loc

        # Deref columns once per locale and not in the deeper loops. Pass needed columns to spread! and transition!
        c_cond       = locdat.cond
        c_status     = locdat.status
        c_agegrp     = locdat.agegrp
        c_sickday    = locdat.sickday
        c_sdcomply   = locdat.sdcomply
        c_variant    = locdat.variant
        c_vaxstatus  = locdat.vaxstatus
        c_variant    = locdat.variant
        c_recovday   = locdat.recovday
        c_vaxrcvd    = locdat.vaxrcvd
        c_vaxday     = locdat.vaxday
        c_deadday    = locdat.deadday

        # other per locale initialization
        poprange = 1:length(locdat)

        # day loop
        for i = 1:ndays  
            inc!(day_ctr, :day)  # increment the simulation day counter
            thisday = day_ctr[:day]
            silent || println("simulation day: ", thisday)

            for case in runcases  # cases that run at the beginning of the day
                case(locdat, socialparams, infectset, sdcases, age_idx_loc; day=thisday, startofday=true, locale=loc)  # TODO extend ages to be any filter for 
            end                                                 # who participates in a given case

            # filter for key people
            idxtime += @elapsed begin
                infect_idx = findall(locdat.status .== infectious)
                contactable_idx = findall(locdat.status .!= dead)
            end

            # if dovax vaccinate (e.g., give shots)
            dovax && (
                        vaxtime += @elapsed vaccinate!(locdat, vaxschedset, contactable_idx, vaxset)
                        )

            # person loop
            for p in infect_idx    
                    
                # is this person ACTIVELY infectious
                spr_sickday = @inbounds c_sickday[p]
                spr_variant = @inbounds c_variant[p][end]
                sendrisk = @inbounds infectset[spr_variant].sendrisk[spr_sickday]
                if sendrisk > 0.0     
                                                    
                    sprtime += @elapsed (
                        spread!(p, contact_vector, sdcases,  socialparams,   
                                    infectset, vaxset, density_factor, dovax, poprange, 
                                    c_cond,
                                    c_status,
                                    c_agegrp,
                                    c_sickday,
                                    c_sdcomply,
                                    c_variant,
                                    c_vaxstatus,
                                    c_recovday,
                                    c_vaxrcvd,
                                    c_vaxday,
                                    ))        
                end

                trtime += @elapsed (
                    transition!(p, infectset, transitionset, vaxset, dovax, noop, trvec,   
                                    c_cond,
                                    c_status,
                                    c_agegrp,
                                    c_sickday,
                                    c_sdcomply,
                                    c_variant,
                                    c_vaxstatus,
                                    c_recovday,
                                    c_vaxrcvd,
                                    c_vaxday,
                                    c_deadday
                                    ))
            end # people loop         
            
            for case in runcases  # cases that run at the end of the day
                case(locdat, socialparams, infectset, sdcases, age_idx_loc; day=day_ctr[:day], startofday=false, locale=loc)  # TODO extend ages to be any filter for 
            end                                                 # who participates in a given case

            # r0 displayed every 10 days
            if showr0 && (mod(day_ctr[:day],10) == 0)   # do we ever want to do this by locale -- maybe
                current_r0 = r0_sim(locdat, age_dist=age_dist, dectree=dectree, socialparams=socialparams, infectparams=infectparams, sdcases=sdcases)
                println("day $(day_ctr[:day]), locale $loc: rt = $current_r0")
            end

            # accumulate simulation statistics in series for plotting: arrays NOT dataframes
            histtime += @elapsed do_history!(locdat, newhist, cumhist, age_idx_loc)

        end # day loop


        histtime += @elapsed begin
            hist_total_agegrps!(newhist, cumhist) # sum agegrps to total for all series groups (by agegrp)
            add_totinfected_series!(newhist, cumhist) 
            add_totvaccinated_series!(newhist, cumhist)
        end

        silent || println("Simulation completed for $(day_ctr[:day]) days for locale $loc.")

    end # locale loop
    end # for totalsimtime

    print_timings(idxtime, vaxtime, sprtime, trtime, histtime, totalsimtime, misc_time)

    return popdat, series
end



################################################################################
#  Update daily history series
################################################################################

@inline function do_history!(locdat, newhist, cumhist, age_idx_loc)  
    thisday = day_ctr[:day]

    @inbounds for age in agegrps

        dat_age = locdat[age_idx_loc[age]]

        # get the source data: status
        status_today = countmap(dat_age.status)    # keys are Enum status
        update_series!(cumhist, newhist, statuses, status_today, age, thisday)

        # get the source data: conditions in (nil, mild, sick, severe)
        filt_infectious = findall(dat_age.status .== infectious)
        if length(filt_infectious) > 0
            sick_today = countmap(dat_age.cond[filt_infectious])  #         keys are enum condition
            update_series!(cumhist, newhist, infectious_cases, sick_today, age, thisday)
        end   

        # get the source data: vaccination
        filt_vaccinated = findall(last.(dat_age.vaxrcvd) .!= :none)
        if length(filt_vaccinated) > 0
            vax_today = countmap(last.(dat_age.vaxrcvd[filt_vaccinated]))  #         keys are symbol
            update_series!(cumhist, newhist, vaxlist, vax_today, age, thisday)
        end

        # get the source data: variants: use filt_infectious from above...
        if length(filt_infectious) > 0
            variant_today = countmap(last.(dat_age.variant[filt_infectious]))    #  keys are symbol
            update_series!(cumhist, newhist, variantlist, variant_today, age, thisday)
        end
        
    end # for age in agegrps

end 


@inline function update_series!(cumhist, newhist, categories, countsdict, age, thisday)

    @inbounds for item in categories
        colname = Symbol(item, "_", age)
        itemcount = get(countsdict, item, 0)
        if itemcount == 0
            continue
        end
        if thisday == 1
            getproperty(cumhist, colname)[thisday] = itemcount
            getproperty(newhist, colname)[thisday] = itemcount
        else
            getproperty(cumhist, colname)[thisday] = itemcount
            getproperty(newhist, colname)[thisday] = (itemcount -  
                    getproperty(cumhist, colname)[thisday-1])
        end
    end
end


@inline function hist_total_agegrps!(newhist, cumhist)
        # runs once per locale
    for item in seriesgroups
        for age in instances(agegrp)
            getproperty(newhist, Symbol(item, "_", "total"))[:] .+= getproperty(newhist, Symbol(item, "_", age))
            getproperty(cumhist, Symbol(item, "_", "total"))[:] .+= getproperty(cumhist, Symbol(item, "_", age))
        end
    end
    
end


# a single locale that already has both new and cum series
@inline function add_totinfected_series!(newhist, cumhist)

    for cond in infectious_cases
        for age in instances(agegrp)
            getproperty(newhist, Symbol(:totinfected, "_", age))[:] .+= getproperty(newhist, Symbol(cond, "_", age))
            getproperty(cumhist, Symbol(:totinfected, "_", age))[:] .+= getproperty(cumhist, Symbol(cond, "_", age))

        end
        getproperty(newhist, Symbol(:totinfected, "_", "total"))[:] .+= getproperty(newhist, Symbol(cond, "_", "total"))
        getproperty(cumhist, Symbol(:totinfected, "_", "total"))[:] .+= getproperty(cumhist, Symbol(cond, "_", "total"))
    end
    
end


@inline function add_totvaccinated_series!(newhist, cumhist)
    for vax in vaxlist
        for age in instances(agegrp)
            getproperty(newhist, Symbol(:totvaccinated, "_", age))[:] .+= getproperty(newhist, Symbol(vax, "_", age))
            getproperty(cumhist, Symbol(:totvaccinated, "_", age))[:] .+= getproperty(cumhist, Symbol(vax, "_", age))

        end
        getproperty(newhist, Symbol(:totvaccinated, "_", "total"))[:] .+= getproperty(newhist, Symbol(vax, "_", "total"))
        getproperty(cumhist, Symbol(:totvaccinated, "_", "total"))[:] .+= getproperty(cumhist, Symbol(vax, "_", "total"))
    end
end


function setx(series, x, cols, rows)
    for col in cols
        getproperty(series, col)[rows] .= x
    end
end


#####################################################################################
#  other functions used in simulation
#####################################################################################

function print_timings(idxtime, vaxtime, sprtime, trtime, histtime, totalsimtime, misc_time)
    println("\nExecution Times")
    @printf "Indexing    %.3f\n" idxtime
    @printf "Vaccination %.3f\n" vaxtime
    @printf "Spread      %.3f\n" sprtime
    @printf "Transition  %.3f\n" trtime
    @printf "History     %.3f\n" histtime
    @printf "Misc Time   %.3f\n" misc_time
    @printf "Total       %.3f\n" totalsimtime
end


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
