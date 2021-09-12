################################
# spread.jl for ilm model
#    social distancing cases
#    spreading the infection
################################


#######################################################################
# social distancing cases
#      struct to hold parameters for defining the case
#      implement the case: 
#           - set social distance compliance for each person
#           - define the contactfactors and touchfactors for the case
#######################################################################          
# mod_90 = sd_gen(start=90,cf=(.2,1.5), tf=(.18,.6),comply=.85)
# str_45 = sd_gen(start=45, comply=.90, cf=(.2,1.0), tf=(.18,.3))
# str_55 = sd_gen(start=55, comply=.95, cf=(.2,1.0), tf=(.18,.3))


Base.@kwdef struct Infectparams
    sendrisk::Vector{Float64}
    recvrisk::Vector{Float64}
end


"""
Method for converting a dict loaded from YAML to this struct
"""
Infectparams(id::Dict) = 
    (Infectparams(
        sendrisk = id[:sendrisk],
        recvrisk = id[:recvrisk]
        )
    )


Base.@kwdef struct Socialparams
    gammashape::Float64
    contactfactors::Dict{Enum, Dict{Enum, Float64}}
    touchfactors::Dict{Enum, Dict{Enum, Float64}}
end


Base.@kwdef struct Spreadcase                 # Base.@kwdef -> use keyword arguments in constructor
    name::Symbol
    day::Int
    cfdelta::Tuple{Float64,Float64}  
    tfdelta::Tuple{Float64,Float64}  
    comply::Float64             # compliance fraction
    cfcase::Dict{Enum, Dict{Enum, Float64}}
    tfcase::Dict{Enum, Dict{Enum, Float64}}
end

function sd_gen(;startday::Int, comply::Float64, cf::Tuple{Float64, Float64},
    tf::Tuple{Float64, Float64}, name::Symbol, include_ages=[])
    function runcase(locale, dat, socialparams, spreadset, sdcases, ages; startofday)   
        s_d_seed!(dat, sdcases, startday, comply, cf, tf, name, include_ages, locale, socialparams, spreadset, ages;
                    startofday=startofday)
    end
end


@inline function s_d_seed!(dat, sdcases, startday, comply, cf, tf, name, include_ages, locale, socialparams, spreadset, ages; startofday)
    @assert 0.0 <= comply <= 1.0  "comply must be floating point in 0.0 to 1.0 inclusive"
    
    startofday || return

    if startday == day_ctr[:day]
        locdat = dat[locale]

        if comply == 0.0  # magic signal: if comply is zero turn off this case for include_ages
            cancel_sd_case!(locdat, sdcases, name, include_ages, ages)
            return
        end

        # create the Spreadcase in sdcases
        sdcases[name] = Spreadcase(
                        name    = name,   # TODO  if we never use this get rid of it
                        day     = startday,   
                        cfdelta = cf,         
                        tfdelta = tf,         
                        comply  = comply,     
                        cfcase  = shifter(socialparams.contactfactors, cf...),  
                        tfcase  = shifter(socialparams.touchfactors, tf...)     
                        )

        # load the sdcomply column of the population table
        # filter1 is everyone who is unexposed, recovered or sick: nil or mild
        filter1 = findall(((locdat.status .== unexposed) .| (locdat.status .== recovered)) .| 
                ((locdat.cond .== nil) .| (locdat.cond .== mild)))
        if (comply == 1.0)   # include everyone in filter1 in this case
            complyfilter = filter1
        else
            complyfilter = sample(filter1, round(Int, comply*length(filter1)), replace=false)
        end
        
        if isempty(include_ages)   # include all include_ages
            locdat.sdcomply[complyfilter] .= name
        else
            byage_idx = intersect(complyfilter, union((ages[i] for i in include_ages)...))
            locdat.sdcomply[byage_idx] .= name
        end
    end
end


function cancel_sd_case!(locdat, sdcases, name, include_ages, ages)
    # filter on who is in this case now
    incase_idx = findall(locdat.sdcomply .== name)

    if isempty(include_ages)   # include all ages
        locdat.sdcomply[incase_idx] .= :none
        delete!(sdcases, name)  # there is no one left in this case...
    else  # only turn it off for some ages
        byage_idx = intersect(incase_idx, union((ages[i] for i in include_ages)...))
        locdat.sdcomply[byage_idx] .= :none
    end    

end

###################################################################
# basic functions for the default definition of spread
###################################################################


"""
    numcontacts(density_factor, gammashape, agegrp, cond, contactfactors)::Int

Returns the number of contacts that someone spreading the disease will make on a day. This
method uses the default contactfactors for the current spreader.
"""
@inline function numcontacts(density_factor, gammashape, agegrp, cond, contactfactors)::Int 
    @inbounds @fastmath scale = density_factor * contactfactors[agegrp][cond]
    @fastmath round(Int,rand(Gamma(gammashape, scale)))
end

"""
    numcontacts(density_factor, gammashape, agegrp, cond, acase::Spreadcase)::Int

Returns the number of contacts that someone spreading the disease will make on a day. This 
method uses the spreadcase applicable to the current spreader.
"""
@inline function numcontacts(density_factor, gammashape, agegrp, cond, acase::Spreadcase)::Int
    @inbounds @fastmath scale = density_factor * acase.cfcase[agegrp][cond]  
    @fastmath round(Int,rand(Gamma(gammashape, scale)))
end


"""
    function istouched(agegrp, lookup, touchfactors)::Bool

Returns true if the contact made was significant to the recipient or false if not.
This method uses the default touchfactors for the current recipient.
"""
@inline function istouched(agegrp, lookup, touchfactors)::Bool
    return @inbounds @fastmath rand(Binomial(1, touchfactors[agegrp][lookup])) == 1
end


"""
    function istouched(agegrp, lookup, acase::Spreadcase)::Bool

Returns true if the contact made was significant to the recipient or false if not.
This method uses the spreadcase for the recipient.
"""
@inline function istouched(agegrp, lookup, acase::Spreadcase)::Bool
    return @inbounds @fastmath rand(Binomial(1, acase.tfcase[agegrp][lookup])) == 1
end


"""
    function isinfected(infectparams, spreadersickday, contactagegrp)::Bool

Returns true if the spreader infected the contact. 
"""
@inline function isinfected(spreadset, vaxset, spreader, contact, locdat)::Bool
    if locdat.vaxrcvd[spreader][end] == :none
        variant = locdat.variant[spreader]
        sendrisk  = spreadset[variant].sendrisk[locdat.sickday[spreader]]
        recvrisk  = spreadset[variant].recvrisk[Int(locdat.agegrp[contact])]
    else
        # sendrisk
            vaxtype = locdat.vaxrcvd[spreader][end]  # the most recent or last vaccine received
            sendrisk = vaxset[vaxtype].sendrisk[locdat.sickday[spreader]]
        # recvrisk
            vaxtype = locdat.vax[contact][end]
            recvrisk = vaxset[vaxtype].recvrisk[Int(locdat.agegrp[contact])]
    end

    @inbounds @fastmath prob = sendrisk * recvrisk            # TODO also vaccinated people will have partially unsusceptible
    return @fastmath rand(Binomial(1, prob)) == 1
end


"""
    spread!(locdat, infect_idx, contactable_idx, sdcases, socialparams, infectparams, density_factor)

Infectious people spread the virus to susceptible people for a single locale. Changes attribute
columns in the population table. Runs social distancing cases.
"""
@inline function spread!(locdat, infect_idx, contactable_idx, sdcases, socialparams,
     spreadset, vaxset, density_factor, dovax, dovariant)

    n_newly_infected = 0

    # retrieve params
    contactfactors = socialparams.contactfactors
    touchfactors   = socialparams.touchfactors
    gammashape     = socialparams.gammashape

    # column aliases as vector v_...
    v_cond     = locdat.cond
    v_status   = locdat.status
    v_agegrp   = locdat.agegrp
    v_sickday  = locdat.sickday
    v_sdcomply = locdat.sdcomply
    dovariant && (v_variant  = locdat.variant)
    dovax && (begin; v_vax      = locdat.vaxrcvd; v_vaxday = locdat.vaxday; end)

    # assign contacts, do touches, do new infections
    @inbounds for spr in infect_idx      # spr is the person who is the spreader
        # determine number of outbound contacts
        contact_param = v_sdcomply[spr] == :none ? contactfactors : sdcases[v_sdcomply[spr]]
        nc = numcontacts(density_factor, gammashape, v_agegrp[spr], v_cond[spr], contact_param)  
        
        # TODO we could keep track of contacts for contact tracing
        # target is the outbound contact reached by the spr (spreader)
        @inbounds @fastmath for target in sample(contactable_idx, nc, replace=true) # people can get contacted more than once
                      # combine a status or a condition value        
            # contactlookup = v_status[target] == infectious ?  v_cond[target] : v_status[target]  # unexposed or recovered
            if v_status[target] == infectious
                continue
            end          
            if v_status[target] == unexposed  # only condition that can get infected   TODO: handle reinfection of recovered
                touch_param = v_sdcomply[target] == :none ? touchfactors : sdcases[v_sdcomply[target]]
                touched = istouched(v_agegrp[target], unexposed, touch_param)   # contactlookup or v_status[target]

                # infection outcome
                if touched         # TODO some recovered people will become susceptible again
                    if isinfected(spreadset, vaxset, spr, target, locdat)
                        v_cond[target] = nil # nil === asymptomatic or pre-symptomatic
                        v_status[target] = infectious
                        v_sickday[target] = 1
                        n_newly_infected += 1
                        dovariant && (v_variant[target] = v_variant[spr])
                    end
                end  # if (touched ...)
            end  # if contactstatus
        end  # for target in sample(...)
    end  # for p in infect_idx

    return n_newly_infected # n_contacts, n_touched, n_newly_infected
end


function make_sick!(dat; cnt, fromage, tocond, tovariant, tosickday=1)

    @assert size(cnt, 1) == size(fromage, 1)

    filt_unexp = optfindall(==(unexposed), dat.status, 1) # must be unexposed

    for i in 1:size(fromage, 1)  # by target age groups

        filt_age = dat.agegrp[filt_unexp] .== fromage[i] # age of the unexposed
        rowrange = 1:cnt[i]
        filt_all = filt_unexp[filt_age][rowrange]

        dat.status[filt_all] .= infectious
        dat.cond[filt_all] .= tocond
        dat.sickday[filt_all] .= tosickday
        dat.variant[filt_all] .= tovariant
    end
end


"""
    r0_sim(; pop=200_000, age_dist=age_dist, dectree=dectree, socialparams=socialparams, infectparams=infectparams, density_factor=1.0, scale=5)
    r0_sim(locdat; age_dist=age_dist, dectree=dectree, socialparams=socialparams, infectparams=infectparams, sdcases=sdcases, density_factor=1.0, scale=5)

Simulates r0 or rt. The first method creates a population and tracks how many infections
are caused by first generation spreaders and NOT spreaders who were infected by the
first generation. The simulates r0

The second method simulates r at time t given the characteristics of the simulation
you are running. This shows how r, reproduction rate, is affected by public health
measures and the characteristics of the population over time. This simulates r(t).
"""
function r0_sim(; pop=200_000, age_dist=age_dist, dectree=dectree, socialparams=socialparams, infectparams=infectparams, density_factor=1.0, scale=5)
    # create simulation population
    r0pop = pop_data(pop)

    # seed spreaders in each age group proportional to age distribution
    cnt_accessible = count(r0pop.status .!= dead)
    cnt_by_agedist = round.(Int, age_dist ./ minimum(age_dist))
    scale = set_by_level(cnt_accessible)
    cnt_by_agedist .*= scale # update with scale

    cnt_spreaders = sum(age_relative)  # COMPARE TO GEN1_INFECTED

    for i in agegrps
        idx = findall(r0pop.agegrp .== i) 

        for j = 1:cnt_by_agedist[Int(i)]
            r0pop.status[idx] = infectious
            r0pop.cond[idx] = nil
            r0pop.sickday[idx] = 1
            idx += 1
        end
    end

    # set infect_idx based on seeding: never update so we measure only 1st gen. spreaders
    gen1_infect_idx = findall(r0pop.status .== infectious)
    gen1_infected = length(gen1_infect_idx)

    sdcases = []   # TODO MAYBE this should be an input based on current context of simulation
    r0_infected = 0

    for i = 1:sickdaylim        
        contactable_idx = findall(r0pop.status .!= dead)
        n_newly_infected = spread!(r0pop, gen1_infect_idx, contactable_idx,  sdcases, socialparams, infectparams, density_factor)  
        infect_idx = findall(r0pop.status .== infectious)
        r0_infected += n_newly_infected
        transition!(r0pop, infect_idx, dectree) 
        gen1_infect_idx = filter(x -> r0pop.status[x] == infectious, gen1_infect_idx)
    end

    r0 =  r0_infected / gen1_infected   # n_newly_infected / cnt_spreaders
    return r0

end


function r0_sim(locdat; age_dist=age_dist, dectree=dectree, socialparams=socialparams, infectparams=infectparams, sdcases=sdcases, density_factor=1.0, scale=5)
    # create simulation population
    r0pop = deepcopy(locdat)

    ignore_idx = optfindall(==(infectious), r0pop.status, 0.5)
    # the following only works because we treat recovered as if they are immune
    r0pop.status[ignore_idx] .= recovered # can't catch what they already have; won't spread for calc of r0

    cnt_accessible = count(r0pop.status .!= dead)
    age_relative = round.(Int, age_dist ./ minimum(age_dist)) # counts by agegrp
    scale = set_by_level(cnt_accessible)
    age_relative .*= scale # update with scale
    cnt_spreaders = sum(age_relative)

    for i in agegrps  # set the spreaders for the r0 simulation
        idx = findall((r0pop.agegrp .== i) .& (r0pop.status .== unexposed))
        for j = 1:age_relative[Int(i)]
            spr = idx[j]
            r0pop.status[spr] = infectious
            r0pop.cond[spr] = nil
            r0pop.sickday[spr] = 1
        end
    end     

    r0_infected = 0 
    for i = 1:sickdaylim      
        infect_idx = findall((r0pop.status .== infectious) .& (r0pop.sickday .> 0))
        contactable_idx = findall(r0pop.status .!= dead)
        # spread!(locdat, infect_idx, contactable_idx, sdcases, socialparams, infectparams, density_factor)                                    
        r0_infected += spread!(r0pop, infect_idx, contactable_idx, sdcases, socialparams, infectparams, density_factor)  

        transition!(r0pop, infect_idx, dectree) 

        # eliminate the new spreaders so we only track the original spreaders
        newsick_idx = findall(r0pop.sickday .== 1)

        # r0pop.status[newsick_idx] .= unexposed
        r0pop.status[newsick_idx] .= recovered # only works because infectious and recovered are treated as immune
    end

    r0 =  r0_infected / cnt_spreaders   # n_newly_infected / cnt_spreaders
    return r0
end


function set_by_level(x, levels=[[1, 300_000], [5, 500_000], [10, 10_000_000_000]])
    ret = 0
    for lvl in levels
        if x <= lvl[2]
            ret = lvl[1]
            break
        end
    end
    return ret
end


function r0_table(n=6, cfstart = 0.9, tfstart = 0.3; socialparams=socialparams, infectparams=infectparams, dt=dt)
    tbl = zeros(n+1,n+1)
    cfiter = [cfstart + (i-1) * .1 for i=1:n]
    tfiter = [tfstart + (i-1) * 0.05 for i=1:n]
    for (j,cf) in enumerate(cfiter)
        for (i,tf) = enumerate(tfiter)
            tbl[i+1,j+1] = r0_sim(socialparams=socialparams, infectparams=infectparams, dt=dt, decpoints=decpoints, shift_contact=(0.2,cf), shift_touch=(.18,tf)).r0
        end
    end
    tbl[1, 2:n+1] .= cfiter
    tbl[2:n+1, 1] .= tfiter
    tbl[:] = round.(tbl, digits=2)
    display(tbl)
    return tbl
end

#=
approximate r0 values from model
using default age distribution
model selects a c_f based on age and infectious case
model selects a t_f based on age and condition (includes unexposed and recovered)
r0 depends on the selection of both c_f and t_f
Note: simulation uses samples so generated values will vary

           c_f
  tf       1.1   1.2      1.3   1.4     1.5   1.6    1.7    1.8    1.9    2.0
           ----------------------------------------------------------
     0.18 | 0.38| 0.38 | 0.42 | 0.46 | 0.49 | 0.51 | 0.55 | 0.57 | 0.59 | 0.64
     0.23 | 0.47| 0.47 | 0.49 | 0.55 | 0.64 | 0.65 | 0.68 | 0.68 | 0.73 | 0.77
     0.28 | 0.53| 0.61 | 0.62 | 0.65 | 0.69 | 0.73 | 0.79 | 0.82 | 0.83 | 0.88
     0.33 | 0.61| 0.66 | 0.7  | 0.79 | 0.8  | 0.83 | 0.9  | 0.95 | 0.99 | 1.04
     0.38 | 0.7 | 0.74 | 0.85 | 0.84 | 0.94 | 0.98 | 1.04 | 1.08 | 1.11 | 1.17
     0.43 | 0.8 | 0.85 | 0.89 | 0.93 | 1.03 | 1.11 | 1.16 | 1.2  | 1.27 | 1.34
     0.48 | 0.88| 0.91 | 0.99 | 1.03 | 1.16 | 1.23 | 1.26 | 1.32 | 1.42 | 1.47
     0.53 | 0.97| 1.06 | 1.08 | 1.18 | 1.26 | 1.27 | 1.42 | 1.47 | 1.52 | 1.61
     0.58 | 1.01| 1.09 | 1.17 | 1.25 | 1.33 | 1.43 | 1.52 | 1.52 | 1.68 | 1.76
     0.63 | 1.11| 1.2  | 1.25 | 1.38 | 1.42 | 1.5  | 1.65 | 1.75 | 1.78 | 1.95


=#|
