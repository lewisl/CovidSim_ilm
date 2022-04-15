
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

    # split up  members of model and initialize
        ndays = model.ndays
        locales = model.locales
        transitionset = model.transitionset  # transition arrays
        trvec = model.trvec # preallocated small vector
        popdat = deepcopy(model.dat["popdat"])   # Copy the population data so model can be reused!!!
        agegrp_idx = model.dat["agegrp_idx"]   # first key is locale
        series = deepcopy(model.series)
        seriescols = series.cols
        geodf = model.geo
        infectset = model.infectset
        socialparams = model.social
        vaxset = model.vaxset
        vaxschedset = model.vaxschedset
        for sched in values(vaxschedset) # k1 is name of a schedule, v1 is instance of struct Vaxsched,
                                    #   vaxesincluded in a field in Vaxsched, which is a dict
            for vax in values(sched.vaxesincluded) # k2 is the key for a vaccine, v2 is the value= an instance of struct Vaxinclude
                vax.doses = vax.starting_doses   # fields of Vaxinclude
            end
        end


    # restart the day counter to zero
    reset!(day_ctr, :day)  # return and reset key to 0 :day leftover from prior runs

    locales = locales   # force local scope to be visible in the loop

    sdcases = Dict{Symbol, Spreadcase}()  # hold definitions of spreadcases


    # execution timers
    vaxtime = 0     # vaccinate
    sprtime = 0     # spread infection
    trtime = 0      # transition infected population through stages of illness
    idxtime = 0     # calculate indices for infectious and susceptible
    histtime = 0    # update history time series
    totalsimtime = 0


    ######################
    # simulation loop
    ######################
    
    # timing begin block
    totalsimtime += @elapsed begin

    for loc in locales     

        silent || println("Simulation starting for location $loc")
        
        # this should be the first and only place to deref the locale (as loc)
        locdat = popdat[loc]  
        newhist = series.new[loc]
        cumhist = series.cum[loc]
        age_idx_loc = agegrp_idx[loc]  # indices by agegrp
        density_factor = geodf[geodf[!, :fips] .== loc, :density_factor][]  # for the loc


        # Deref columns once per locale and not in the deeper loops. Pass these columns to spread! and transition!
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


        for i = 1:ndays
            inc!(day_ctr, :day)  # increment the simulation day counter
            silent || println("simulation day: ", day_ctr[:day])

                for case in runcases  # cases that run at the beginning of the day
                    case(locdat, socialparams, infectset, sdcases, age_idx_loc; day=day_ctr[:day], startofday=true, locale=loc)  # TODO extend ages to be any filter for 
                end                                                 # who participates in a given case

                # filter for key people
                idxtime += @elapsed begin
                    # @bp
                    infect_idx = findall(locdat.status .== infectious)
                    contactable_idx = findall(locdat.status .!= dead)
                end

                # if dovax vaccinate (e.g., give shots)
                dovax && (
                            vaxtime += @elapsed vaccinate!(locdat, vaxschedset, contactable_idx, vaxset)
                         )

                # two fundamental steps of the simulation: spread! and transition!
                sprtime += @elapsed spread!(infect_idx, contactable_idx, sdcases,  socialparams, 
                                            infectset, vaxset, density_factor, dovax, dovariant;
                                            c_cond       = c_cond,
                                            c_status     = c_status,
                                            c_agegrp     = c_agegrp,
                                            c_sickday    = c_sickday,
                                            c_sdcomply   = c_sdcomply,
                                            c_variant    = c_variant,
                                            c_vaxstatus  = c_vaxstatus,
                                            c_recovday   = c_recovday,
                                            c_vaxrcvd    = c_vaxrcvd,
                                            c_vaxday     = c_vaxday
                                            )   
                trtime += @elapsed transition!(infect_idx, infectset, transitionset, vaxset, dovax, dovariant; trvec=trvec,
                                               c_cond       = c_cond,
                                               c_status     = c_status,
                                               c_agegrp     = c_agegrp,
                                               c_sickday    = c_sickday,
                                               c_sdcomply   = c_sdcomply,
                                               c_variant    = c_variant,
                                               c_vaxstatus  = c_vaxstatus,
                                               c_recovday   = c_recovday,
                                               c_vaxrcvd    = c_vaxrcvd,
                                               c_vaxday     = c_vaxday,
                                               c_deadday    = c_deadday
                                                )                        

                for case in runcases  # cases that run at the end of the day
                    case(locdat, socialparams, infectset, sdcases, age_idx_loc; day=day_ctr[:day], startofday=false, locale=loc)  # TODO extend ages to be any filter for 
                end                                                 # who participates in a given case

                # r0 displayed every 10 days
                if showr0 && (mod(day_ctr[:day],10) == 0)   # do we ever want to do this by locale -- maybe
                    current_r0 = r0_sim(locdat, age_dist=age_dist, dectree=dectree, socialparams=socialparams, infectparams=infectparams, sdcases=sdcases)
                    println("day $(day_ctr[:day]), locale $loc: rt = $current_r0")
                end

                # accumulate simulation statistics in series for plotting: arrays NOT dataframes
                histtime += @elapsed do_history!(locdat, newhist, cumhist, age_idx_loc, seriescols)

        end # for i in 1:n_days

        histtime += @elapsed begin
            hist_total_agegrps!(newhist, cumhist, seriescols) # sum agegrps to total for all series groups (by agegrp)
            add_totinfected_series!(newhist, cumhist, seriescols) 
            add_totvaccinated_series!(newhist, cumhist, seriescols)
        end

        silent || println("Simulation completed for $(day_ctr[:day]) days for locale $loc.")

    end # for loc in locales
    end # totalsimtime

    print_timings(idxtime, vaxtime, sprtime, trtime, histtime, totalsimtime)

    return popdat, series
end



################################################################################
#  Update daily history series
################################################################################

@views function do_history!(locdat, newhist, cumhist, age_idx_loc, seriescols)  # cumhist, newhist,

    @inbounds for age in instances(agegrp)

        # get the source data: status
        dat_age = locdat[age_idx_loc][age]
        status_today = @inbounds countmap(dat_age.status)    # cumulative position for thisday, keys are Enum status

        # get the source data: conditions in (nil, mild, sick, severe)
        filt_infectious = findall(dat_age.status .== infectious)
        if length(filt_infectious) > 0
            sick_today = @inbounds countmap(dat_age.cond[filt_infectious])  # keys are enum condition
        else   # there can be days when no one is infected
            sick_today = Dict()
        end

        # get the source data: vaccination
        filt_vaccinated = findall(last.(dat_age.vaxrcvd) .!= :none)
        if length(filt_vaccinated) > 0
            vax_today = @inbounds countmap(last.(dat_age.vaxrcvd[filt_vaccinated]))  # keys are symbol
        else
            vax_today = Dict()
        end

        # get the source data: variants: use filt_infectious from above...
        if length(filt_infectious) > 0
            variant_today = @inbounds countmap(last.(dat_age.variant[filt_infectious]))
        else
            variant_today = Dict()
        end
        
        #
        # cumulative and new series
        #

        int_age = Int(age)
        thisday = day_ctr[:day]
    
        saveseries!(cumhist, newhist, statuses, status_today, int_age, seriescols, thisday)
        saveseries!(cumhist, newhist, infectious_cases, sick_today, int_age, seriescols, thisday)
        saveseries!(cumhist, newhist, vaxlist, vax_today, int_age, seriescols, thisday)
        saveseries!(cumhist, newhist, variantlist, variant_today, int_age, seriescols, thisday)
        
    end # for age in agegrps

end # function


function saveseries!(cumdat, newdat, categories, countsdict, int_age, seriescols, thisday)

    @inbounds for item in categories
        seriescol = seriescols[Symbol(item)][int_age]
        itemcount = get(countsdict, item, 0)
        if thisday == 1
            cumdat[thisday, seriescol] = itemcount
            newdat[thisday, seriescol] = itemcount  # initialize 1st day of new
        else
            cumdat[thisday, seriescol] = itemcount
            newdat[thisday, seriescol] = (    # day 2... do cum(day n) - cum(day n-1)
                cumdat[thisday, seriescol]
                - cumdat[thisday - 1, seriescol]
                )
        end
    end
end


function hist_total_agegrps!(newhist, cumhist, seriescols)
    cols = seriescols

    @inbounds for cond in allconds  # infectious cases and statuses
        colgroup = cols[Symbol(cond)]
        newhist[:, colgroup[totalcol]] = sum(newhist[:, colgroup[collect(Int.(agegrps))]], dims=2)
        cumhist[:, colgroup[totalcol]] = sum(cumhist[:, colgroup[collect(Int.(agegrps))]], dims=2)
    end

    @inbounds for vax in vaxlist
        colgroup = cols[Symbol(vax)]
        newhist[:, colgroup[totalcol]] = sum(newhist[:, colgroup[collect(Int.(agegrps))]], dims=2)
        cumhist[:, colgroup[totalcol]] = sum(cumhist[:, colgroup[collect(Int.(agegrps))]], dims=2)
    end

    @inbounds for variant in variantlist
        colgroup = cols[Symbol(variant)]
        newhist[:, colgroup[totalcol]] = sum(newhist[:, colgroup[collect(Int.(agegrps))]], dims=2)
        cumhist[:, colgroup[totalcol]] = sum(cumhist[:, colgroup[collect(Int.(agegrps))]], dims=2)
    end
    
end


# a single locale that already has both new and cum series
function add_totinfected_series!(newhist, cumhist, seriescols)
    cols = seriescols
    # for new
    @views begin
        n = size(newhist,1)
        newhist[:,cols[:totinfected]] = ( (newhist[:,cols[:unexposed]] .< 0 ) .*
                                                          abs.(newhist[:,cols[:unexposed]]) ) 
        cumsum!(cumhist[:, cols[:totinfected]], newhist[:, cols[:totinfected]], dims=1)  
    end
    
end


function add_totvaccinated_series!(newhist, cumhist, seriescols)
    cols = seriescols
    @views begin
        n = size(newhist, 1)
        newhist[:, cols[:totvaccinated]] = (newhist[:, cols[:JnJ]] .+ newhist[:, cols[:Pfizer]] .+ 
                                                    newhist[:, cols[:Moderna]])
        cumhist[:, cols[:totvaccinated]] = (cumhist[:, cols[:JnJ]] .+ cumhist[:, cols[:Pfizer]] .+ 
                                                    cumhist[:, cols[:Moderna]])
    end
end



#####################################################################################
#  other functions used in simulation
#####################################################################################

function print_timings(idxtime, vaxtime, sprtime, trtime, histtime, totalsimtime)
    println("\nExecution Times")
    @printf "Indexing    %.3f\n" idxtime
    @printf "Vaccination %.3f\n" vaxtime
    @printf "Spread      %.3f\n" sprtime
    @printf "Transition  %.3f\n" trtime
    @printf "History     %.3f\n" histtime
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
