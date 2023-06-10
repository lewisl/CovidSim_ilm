
####################################################################################
#   simulation runner: ILM Model
####################################################################################


function buildsim(ndays, locales;
    day1 = Date("2020-01-01", "yyyy-mm-dd"),    # first calendar day of simulation
    dovax = false,                              # vaccinations for people
    paramdir = "../sample_parameters",          # a directory of required parameters
    geofilename = "../data/geo2data.csv", 
    socialfilename = "socialparams.yml",
    vaccinefilename = "vaccines.yml",
    scheddir="vaccine_schedule",
    variantfilename = "variants.yml")

    locales = locales isa Int ? [locales] : locales

    model = setup(ndays, locales; 
        day1=day1,
        dovax=dovax, 
        paramdir=paramdir,
        geofilename=geofilename, 
        socialfilename=socialfilename,
        vaccinefilename=vaccinefilename,
        scheddir=scheddir,
        variantfilename=variantfilename,
        )

        #=
        model = (ndays=ndays, day1=day1, locales=locales, dat=dat, series=series, geo=geodata, 
                progressionset=progressionset, vaxset=vaxset, vaxschedset=vaxschedset, infectset=infectset, 
                social=socialparams, trvec=trvec, vaxlist=vaxlist, variantlist=variantlist, 
                seriescolnames=seriescolnames)  
        =#    

    return model
end


function buildsim(yaml_model)
    model = setup(yaml_model)
end


function runsim(model; 
            runcases=[], 
            showr0 = false, 
            silent=true, 
            dovax=false,
            vaxscheds=:none
            )

    # split up  members of model
    ndays = model.ndays
    day1 = model.day1
    locales = model.locales
    progressionset = model.progressionset  # progression arrays
    trvec = model.trvec # preallocated small vector
    popdat = deepcopy(model.dat.popdat)   # Copy the population data so model can be reused!!!
    agegrp_idx = model.dat.agegrp_idx   # first key is locale
    series = deepcopy(model.series)  # dict of locales => namedtuple(.cum, .new), TypedTable of history columns
    geodf = model.geo
    infectset = model.infectset
    variantlist = model.variantlist
    socialparams = model.social
    vaxset = model.vaxset
    vaxlist = model.vaxlist
    indoor_seq = model.indoor_seq
    seriescolnames = model.seriescolnames
    vaxschedset = model.vaxschedset

    # initialize some factors
    for sched in values(vaxschedset) 
        for vax in values(sched.vaxesincluded) 
            vax.doses = vax.starting_doses   
        end
    end

    sdcases = Dict{Symbol, SpreadCase}()  # hold definitions of spreadcases

    if showr0
        r0sim_output = IOBuffer()  # bullshit to create correct output in VS Code notebooks
    end

    # restart the day counter to zero
    reset!(DAY_CTR, :day)  # return and reset key to 0 :day leftover from prior runs

    # execution timers
    vaxtime = 0     # vaccinate
    sprtime = 0     # spread infection
    trtime = 0      # progression infected population through stages of illness
    idxtime = 0     # calculate indices for infectious and susceptible
    histtime = 0    # update history time series
    totalsimtime = 0


    ######################
    # simulation loop
    ######################
    totalsimtime += @elapsed begin

    for loc in locales     

        silent || println("Simulation starting for location $loc")
        
        # first and only place to deref the locale (as loc)
        locdat = popdat[loc]  
        newhist = series[loc].new
        cumhist = series[loc].cum
        age_idx_loc = agegrp_idx[loc]  # indices by agegrp
        density_factor = geodf.density_factor[geodf.fips .== loc][1]
        caldays = cumhist.caldays
        indoor_seq = indoor_seq[loc]

        # Deref columns once per locale and not in the deeper loops. Pass needed columns to spread! and progression!
        c_cond       = locdat.cond
        c_status     = locdat.status
        c_agegrp     = locdat.agegrp
        c_duration   = locdat.duration
        c_sdcase     = locdat.sdcase
        c_variant    = locdat.variant
        c_vaxstatus  = locdat.vaxstatus
        c_sickday    = locdat.sickday
        c_variant    = locdat.variant
        c_recovday   = locdat.recovday
        c_vaxrcvd    = locdat.vaxrcvd
        c_vaxday     = locdat.vaxday
        c_deadday    = locdat.deadday

        # other per locale initialization
        poprange = 1:length(locdat)

        # day loop:  simulation time step is one day
        for i = 1:ndays  
            inc!(DAY_CTR, :day)  # increment the simulation day counter
            today = DAY_CTR[:day]
            silent || println("simulation day: ", today)

            for case in runcases  # cases that run at the beginning of the day
                case(locdat, socialparams, infectset, sdcases, age_idx_loc; day=today, startofday=true, locale=loc) 
            end    ## TODO extend ages to be any filter for who participates in a given case

            # filter for key people
            idxtime += @elapsed infect_idx = findall(locdat.status .== :infectious) # all the sick and maybe infectious
            
            # if dovax vaccinate (e.g., give shots)
            dovax && begin
                    vaxtime += @elapsed vaccinate!(vaxschedset, vaxset, vaxscheds,
                                        c_status,
                                        c_agegrp,
                                        c_vaxstatus,
                                        c_recovday,
                                        c_vaxrcvd,
                                        c_vaxday
                                    )
                     end

            # person loop
            @inbounds for p in infect_idx    # p is an infected person who potentially spreads virus

                sprtime += @elapsed begin
                    # is this person ACTIVELY infectious
                    spr_duration = c_duration[p]  # duration determines if spreader is really able to spread the virus
                    spr_variant = c_variant[p][end]
                    
                    sendrisk = infectset[spr_variant].sendrisk[spr_duration]
                    
                    if sendrisk > 0.0     
                        # transmission of the virus
                        spread!(p, today, sdcases,  socialparams, 
                                infectset, vaxset, density_factor, indoor_seq, poprange, 
                                    # columns of simulation data table, rows are persons
                                    c_cond,
                                    c_status,
                                    c_agegrp,
                                    c_duration,
                                    c_sdcase,
                                    c_sickday,
                                    c_variant,
                                    c_vaxstatus,
                                    c_recovday,
                                    c_vaxrcvd,
                                    c_vaxday)        
                    end
                end  # sprtime

                trtime += @elapsed begin
                    # progression of the disease for each infected person
                    progression!(p, infectset, progressionset, vaxset, dovax, trvec,   
                                    # columns of simulation data table, rows are persons
                                    c_cond,
                                    c_status,
                                    c_agegrp,
                                    c_duration,
                                    c_sdcase,
                                    c_variant,
                                    c_vaxstatus,
                                    c_recovday,
                                    c_vaxrcvd,
                                    c_vaxday,
                                    c_deadday
                                    ) 
                    end
            end # people loop         
            
            for case in runcases  # cases that run at the end of the day
                case(locdat, socialparams, infectset, sdcases, age_idx_loc; day=today, startofday=false, locale=loc) 
            end   # TODO extend ages to be any filter for who participates in a given case

            # r0 displayed every 10 days
            if showr0 && (mod(DAY_CTR[:day],10) == 0)   # do we ever want to do this by locale -- maybe
                current_r0 =  r0_sim(locdat, progressionset, trvec, infectset, vaxset, 
                        socialparams, dovax, :base, density_factor, 3)
                write(r0sim_output, "Day: $(DAY_CTR[:day])")
                write(r0sim_output, " Locale: $loc")
                write(r0sim_output, " Current r(t): $current_r0 \n"); 
            end

            histtime += @elapsed collect_history!(locdat, newhist, cumhist, age_idx_loc, today, vaxlist, variantlist, seriescolnames)

        end # day loop

        # calculate history totals by agegroup, infected for all conditions, vaccinated for all vaccines
        histtime += @elapsed begin
            update_total_agegrps!(newhist, cumhist, seriescolnames) # sum agegrps to total for all series groups (by agegrp)
            update_totinfected_series!(newhist, cumhist, seriescolnames) 
            update_totvaccinated_series!(newhist, cumhist, vaxlist, seriescolnames)
        end

        silent || println("Simulation completed for $(DAY_CTR[:day]) days for locale $loc.")

    end # locale loop
    end # for totalsimtime

    if showr0                  # bullshit for correct output in VS Code notebooks
        flush(r0sim_output)
        printthis = String(take!(r0sim_output))
        close(r0sim_output)
    end

    flush(stdout); print_timings(idxtime, vaxtime, sprtime, trtime, histtime, totalsimtime); flush(stdout)
    showr0 && begin; flush(stdout); println(printthis); end

    return popdat, series  # final state of population matrix, history data for simulation run by day
end



################################################################################
#  Update daily history series
################################################################################

@inline function collect_history!(locdat, newhist, cumhist, age_idx_loc, today, vaxlist, variantlist, seriescolnames)  

    statuscol = locdat.status
    condcol = locdat.cond
    vaxcol = locdat.vaxrcvd
    variantcol = locdat.variant
    scn = seriescolnames

    @inbounds @fastmath for age in AGEGRPS

        age_idx = age_idx_loc[age]
        # dat_age = sourcedat[age_idx]   

        # vectors of count of outcome by trait column for status, cond, vax, and variant
        status_today = zeros(Int, 4)
        sick_today = zeros(Int, 4)
        vax_today = zeros(Int, 3)
        variant_today = zeros(Int, 6)

        # iterate through each person p in the age group
        for p in age_idx

            # get the source data: status
            status_today[mapstatus(statuscol[p])] += 1

            # get the source data: conditions in (nil, mild, sick, severe)
            # get the source data: variant in <list of variants.
            if statuscol[p] == :infectious
                sick_today[mapcondition(condcol[p])] += 1
                variant_today[variantdict[last(variantcol[p])]] += 1
            end
 
            # get the source data: vaccination
            vax_of_p = last(vaxcol[p])
            if vax_of_p != :none
                vax_today[vaxdict[vax_of_p]] += 1
            end
                
        end # for p in age_idx

        # put the counts into the history series
        update_series!(cumhist, newhist, scn, STATUSES, status_today, age, today, group=:statuscols, intmapper=mapstatus)
        update_series!(cumhist, newhist, scn, INFECTIOUS_CASES, sick_today, age, today, group=:condcols, intmapper=mapcondition)
        update_series!(cumhist, newhist, scn, vaxlist, vax_today, age, today, group=:vaxcols, mapdict=vaxdict)
        update_series!(cumhist, newhist, scn, variantlist, variant_today, age, today, group=:variantcols, mapdict=variantdict)

    end # for age in AGEGRPS

    # :unexposed special case:  no new people on day 1
    if today == 1  
        for colname in values(seriescolnames[:statuscols][:unexposed])
            getproperty(newhist, colname)[today] = 0
        end
    end

end 


@inline function update_series!(cumhist, newhist, scn, categories, countsvec, age, today; group, intmapper=mapviadict, mapdict=Dict())

    @inbounds @fastmath for item in categories
        seriescol = scn[group][Symbol(item)][Symbol(age)]
        itemcount = isempty(mapdict) ? countsvec[intmapper(item)] : countsvec[intmapper(mapdict, item)]
        if itemcount == 0
            continue
        end
        if today == 1
            getproperty(cumhist, seriescol)[today] = itemcount
            getproperty(newhist, seriescol)[today] = itemcount
        else
            getproperty(cumhist, seriescol)[today] = itemcount
            getproperty(newhist, seriescol)[today] = (itemcount -  
                    getproperty(cumhist, seriescol)[today-1])
        end
    end
end


"""
Sum all the series columns for all ages for all groups and items into a :total column by group and item.
"""
@inline function update_total_agegrps!(newhist, cumhist, scn)
    # runs once per locale
    @fastmath @inbounds for group in keys(scn)
        for item in keys(scn[group])
            totalcol = scn[group][item][:total]
            sumcols = Tuple(scn[group][item][age] for age in AGEGRPS)
            getproperty(newhist, totalcol)[:] .= .+(columns(getproperties(newhist, sumcols))...)  # a tuple of column names
            getproperty(cumhist, totalcol)[:] .= .+(columns(getproperties(cumhist, sumcols))...)
        end
    end
end


"""
Sum all the series columns for all ages and total across ages for all infectious_cases into :totinfected series group of columns.
"""
@inline function update_totinfected_series!(newhist, cumhist, scn)
    @fastmath @inbounds for age in AGENAMES  # for each age and "total"
        totalcol = scn[:condcols][:totinfected][age]
        sumcols = Tuple(scn[:condcols][cond][age] for cond in Symbol.(INFECTIOUS_CASES)) # tuple of all of condition column names
        # not the most obvious below, but faster than looping one column at a time!
        getproperty(newhist, totalcol)[:] .= .+(columns(getproperties(newhist, sumcols))...)  # sum all of the infectious_cases columns
        getproperty(cumhist, totalcol)[:] .= .+(columns(getproperties(cumhist, sumcols))...)
    end 
end


"""
Sum all the series columns for all ages and total across ages for all vaccines into :totvaccinated series group of columns.
"""
@inline function update_totvaccinated_series!(newhist, cumhist, vaxlist, scn)
    @fastmath @inbounds for age in AGENAMES  # for each age and "total"  # Tuple(scn[:vaxcols][vax][age] for vax in vaxlist)
        totalcol = scn[:vaxcols][:totvaccinated][age]
        sumcols = Tuple(scn[:vaxcols][vax][age] for vax in vaxlist) # tuple of all of the vax column names
        getproperty(newhist, totalcol)[:] .= .+(columns(getproperties(newhist, sumcols))...)        # sum all of the vax columns
        getproperty(cumhist, totalcol)[:] .= .+(columns(getproperties(cumhist, sumcols))...)
    end
end


#####################################################################################
#  other functions used in simulation
#####################################################################################

function print_timings(idxtime, vaxtime, sprtime, trtime, histtime, totalsimtime)
    println("\nExecution Times")
    @printf "Indexing     %.3f\n" idxtime
    @printf "Vaccination  %.3f\n" vaxtime
    @printf "Spread       %.3f\n" sprtime
    @printf "Progression  %.3f\n" trtime
    @printf "History      %.3f\n" histtime
    @printf "Total        %.3f\n" totalsimtime
end


function cleanup_stash(stash)
    for k in keys(stash)
        delete!(stash, k)
    end
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
