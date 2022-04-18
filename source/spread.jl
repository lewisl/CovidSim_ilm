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
    recovery_immunity::Dict{Symbol, Float64}
    immunehalflife::Int64
    basemultiplier::Float64
end

        """
        Method for converting a dict loaded from YAML to this struct
        """
        function Infectparams(indict::Dict{Symbol, Any})
            Infectparams(
                sendrisk = indict[:sendrisk],
                recvrisk = indict[:recvrisk],
                recovery_immunity = indict[:recovery_immunity],
                immunehalflife = indict[:immunehalflife],
                basemultiplier = indict[:basemultiplier]
                )
        end


Base.@kwdef struct Socialparams
    gammashape::Float64
    contactfactors::Matrix{Float64}     
    touchfactors::Matrix{Float64}     
end


Base.@kwdef struct Spreadcase                 # Base.@kwdef -> use keyword arguments in constructor
    name::Symbol
    day::Int
    cfdelta::Tuple{Float64,Float64}  
    tfdelta::Tuple{Float64,Float64}  
    comply::Float64             # compliance fraction
    cfcase::Matrix{Float64}
    tfcase::Matrix{Float64}
end

function sd_gen(;startday::Int, comply::Float64, cf::Tuple{Float64, Float64},
                tf::Tuple{Float64, Float64}, name::Symbol, include_ages=[])
    function caserunner(locdat, socialparams, infectset, sdcases, age_idx_loc; startofday)   
        s_d_seed!(locdat, sdcases, startday, comply, cf, tf, name, include_ages, socialparams, infectset, age_idx_loc;
                    startofday=startofday)
    end
end


@inline function s_d_seed!(locdat, sdcases, startday, comply, cf, tf, name, include_ages, socialparams, infectset, age_idx_loc; startofday)
    @assert 0.0 <= comply <= 1.0  "comply must be floating point in 0.0 to 1.0 inclusive"
    
    startofday || return

    if startday == day_ctr[:day]

        if comply == 0.0  # magic signal: if comply is zero turn off this case for include_ages
            cancel_sd_case!(locdat, sdcases, name, include_ages, age_idx_loc)
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
            byage_idx = intersect(complyfilter, union((age_idx_loc[i] for i in include_ages)...))
            locdat.sdcomply[byage_idx] .= name
        end
    end
end


function cancel_sd_case!(locdat, sdcases, name, include_ages, age_idx_loc)
    # filter on who is in this case now
    incase_idx = findall(locdat.sdcomply .== name)

    if isempty(include_ages)   # include all ages
        locdat.sdcomply[incase_idx] .= :none
        delete!(sdcases, name)  # there is no one left in this case...
    else  # only turn it off for some ages
        byage_idx = intersect(incase_idx, union((age_idx_loc[i] for i in include_ages)...))
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
    @inbounds @fastmath scale = density_factor * contactfactors[mapcondition(cond), mapagegrp(agegrp)]
    @fastmath round(Int,rand(Gamma(gammashape, scale)))
end

"""
    numcontacts(density_factor, gammashape, agegrp, cond, acase::Spreadcase)::Int

Returns the number of contacts that someone spreading the disease will make on a day. This 
method uses the spreadcase applicable to the current spreader.
"""
@inline function numcontacts(density_factor, gammashape, agegrp, cond, acase::Spreadcase)::Int
    @inbounds @fastmath scale = density_factor * acase.cfcase[mapcondition(cond), mapagegrp(agegrp)]  
    @fastmath round(Int,rand(Gamma(gammashape, scale)))
end


"""
    function istouched(agegrp, lookup, touchfactors)::Bool
    function istouched(agegrp, lookup, acase::Spreadcase)::Bool

Returns true if the contact made was significant to the recipient or false if not.
First method uses the default touchfactors for the current recipient.
Second method uses the spreadcase for the recipient.
"""
@inline function istouched(agegrp, lookup, touchfactors)::Bool
    return @inbounds @fastmath rand(Binomial(1, touchfactors[maptouch(lookup), mapagegrp(agegrp)])) == 1
end



@inline function istouched(agegrp, lookup, acase::Spreadcase)::Bool
    return @inbounds @fastmath rand(Binomial(1, acase.tfcase[maptouch(lookup), mapagegrp(agegrp)])) == 1
end


"""
    function isinfected(infectparams, spreadersickday, contactagegrp)::Bool

Returns true if the spreader infected the contact. 
"""
@inline function isinfected(risk)::Bool
    return @fastmath rand(Binomial(1, risk)) == 1
end

#############################################################################
#
#  effect of immunity from recovery and vaccination
#
#############################################################################

# decay functions for immunity for decline from 1.0 to lower positive limit of function
# multiply times max immunity if less than 1.0        



"""
    riseup(days_since_shot, delay_days, lower, upper)
    riseup14(t)     

Vaccine infectreduce ramps up after receiving a shot.
This returns a value between 0.0 and 1.0 which 
must be multiplied times the vaccine specific infectreduce 
because this function only represents the time-based change.

The second method is curried to only require days_since_shot as
input. The other arguments to this method are: 14, 0.0, 1.0
"""
@fastmath function riseup(days_since_shot, delay_days, lower, upper)
    clamp(days_since_shot * (1.0 / delay_days), lower, upper)
end

riseup14(t) = riseup(t, 14, 0.0, 1.0)  # curried to only input the day as time t



# gradual decay of vaccine infectreduce based on assumed half-life

@fastmath function lindecay(t, h, lower)
    y = 0.5 ./ -h * t  + 1.0
    y = y < lower ? lower : y
end

expdecay(t,h) = exp(-(log(2)/h) * t)  

sigdecay(t, h; csig=5.0) = 1.0 / (1.0 + exp.((t - h)/(t / csig + (h / csig))))    

tbrk(h, lower) = 2.0 * h - (2.0 * h * lower)

intercept(t, hl, lower) = -0.3 * t / hl + 1.0

function lindecay2(t,hl,lower1, lower2)
    f1 = -t * 0.5 / hl + 1.0
    if  f1 >= lower1
        f1
    else
        clamp(-t * 0.2 / hl + intercept(tbrk(hl, lower1), hl, lower1), lower2, 1.0)
    end
end


function lindecayarr(t::AbstractVector{T} where T, hl, lower1, lower2)
    arr = zeros(size(t,1))
    icept = intercept(tbrk(hl, lower1), hl, lower1)
    for i = eachindex(arr)
        f1 = -t[i] * 0.5 / hl + 1.0
        if  f1 >= lower1
            arr[i] = f1
        else
            arr[i] = clamp(-t[i] * 0.2 / hl + icept, lower2, 1.0)
        end
    end
    return arr
end


@inline @fastmath function vaxmodifier(full_effect_days, today, lastshotday, halflife; rise_lower=0.5, decay_lower=0.05)
    # combines the effect of the rise to full infectreduce post shot with
        # the decay in infectreduce over time: based on current date
    @assert today >= lastshotday "today's date must be >= to day of most recent shot"

    rise = riseup(today - lastshotday, full_effect_days, rise_lower, 1.0)

    days_after_full_effect = today - (lastshotday + full_effect_days)     #clamp(today - (lastshotday + full_effect_days), 0, Int)
    decay = lindecay(days_after_full_effect, halflife, decay_lower)
    return rise * decay
end


function sigmoidshift(x; risk_discount=0.2)
    sigmoid(
            shifter(
                    clamp(x, 0.0, 1.0 + risk_discount),
                    0.0, 1.0, -4.0, 4.0
                    )
            )
end

@inline function simpleclamp(x)
    clamp(x, 0.0, 1.0)
end

squashfunc = simpleclamp


"""
    spr_vaxeffect(infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday)

Immunity from vaccination for a single person.
"""
@inline @fastmath function spr_vaxeffect(infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday)

    today = day_ctr[:day]
    oneshotfactor = 0.85   # TODO yet another parameter to put somewhere...!

    # person's vaccine conditions
    days_after_vax = today - vaxday
    @assert today >= days_after_vax "today's date must be >= to day of most recent shot"

    # vaccine characteristics
    halflife = vaxset[vaxrcvd].halflife
    full_effect_days = vaxset[vaxrcvd].full_effect_days
    infectfactor = vaxset[vaxrcvd].infectfactor
    vaxeffect = @inbounds vaxset[vaxrcvd].infectreduce[vaxstatus][spr_variant]

    # rise and decay of vaccine effectiveness
    rise_lower=0.5    # TODO need to make these inputs somewhere
    decay_lower=0.05
    rise = riseup(today - days_after_vax, full_effect_days, rise_lower, 1.0)
    days_after_full_effect = today - (days_after_vax + full_effect_days)     #clamp(today - (lastshotday + full_effect_days), 0, Int)
    decay = lindecay(days_after_full_effect, halflife, decay_lower)
    vaxmod = rise * decay

    factor = 1.0 - (vaxmod * vaxeffect * infectfactor) 

    return factor
end


"""
    spr_recoveffect(recovday, targ_variant, spr_variant, infectset)

Immunity from recovery for a single person.
"""
@inline function spr_recoveffect(recovday, targ_variant, spr_variant, infectset)

        today = day_ctr[:day]
        days_post_recov = today - recovday 

        if days_post_recov > 0   # TODO should be an assert: does this run day of or day after recovery?
            # get the max immunity
            immstrength = infectset[targ_variant].recovery_immunity[spr_variant]

            # get the declined value
            immhalflife = infectset[targ_variant].immunehalflife
            immdecline = lindecay(days_post_recov, immhalflife, 0.05)

            factor = 1.0 - (immdecline * immstrength)
        else
            factor = 1.0
        end

    return factor
end



@inline @fastmath function infectrisk(infectset, spr_variant, spr_sickday, targ_agegrp, recovfactor, vaxfactor)

    # spreader person characteristics
    sendrisk = @inbounds infectset[spr_variant].sendrisk[spr_sickday]

    # target person characteristics
    recvrisk = @inbounds infectset[spr_variant].recvrisk[Int(targ_agegrp)]

    combinedfactor = recvrisk * sendrisk * min(recovfactor, vaxfactor)
    riskfactor = squashfunc(combinedfactor)  
end


# impact of differential infectiousness for vaccine or variant
spreadin(risk) = 8.0 * risk - 4.0  # rescale risk to -4.0, 4.0
sigmoid(x) = 1.0 / (1.0 + exp(-x))  # smoosh input to 0.0, 1.0

altrisk(risk) = sigmoid(spreadin(risk))


"""
    spread!(locdat, infect_idx, contactable_idx, sdcases, socialparams, infectparams, density_factor)

Infectious people spread the virus to susceptible people for a single locale. Changes attribute
columns in the population table. Runs social distancing cases.
"""
@inline function spread!(spr::Int, contact_vector::Vector{Int}, sdcases, socialparams,
     infectset, vaxset, density_factor, dovax, dovariant,     
        c_pid,
        c_cond,
        c_status,
        c_agegrp,
        c_sickday,
        c_sdcomply,
        c_variant,
        c_vaxstatus,
        c_recovday,
        c_vaxrcvd,
        c_vaxday
     )

    # retrieve params
    contactfactors = socialparams.contactfactors
    touchfactors   = socialparams.touchfactors
    gammashape     = socialparams.gammashape

    # initialize
    tocond = nil
    n_newly_infected = 0
    tovariant = :base
    targ_agegrp = age0_19
    history_changes = ()
    target_status = unexposed

    # determine number of outbound contacts 
    contact_param = c_sdcomply[spr] == :none ? contactfactors : sdcases[c_sdcomply[spr]]
    nc = @inbounds numcontacts(density_factor, gammashape, c_agegrp[spr], c_cond[spr], contact_param)  
    sample!(c_pid, contact_vector)
    sel = 0
    for i = 1:nc
        if sel >= length(contact_vector)
            sample!(c_pid, contact_vector)  # draw another sample
            sel = 1
        else
            sel += 1
        end
    
        target = @inbounds contact_vector[sel]
        target_status = @inbounds c_status[target]

        @inbounds if in(target_status, (unexposed, recovered))  # only conditions that can get infected   
            # choose the touch_param for the social distancing case or the input social parameters
            touch_param = c_sdcomply[target] == :none ? touchfactors : sdcases[c_sdcomply[target]]
            touched = @inbounds istouched(c_agegrp[target], unexposed, touch_param)   # is the contact significant?

            # infection outcome
            if touched  # if the contact is consequential
                # gather characteristics of target and spreader
                    recovday = @inbounds c_recovday[target][end]
                    targ_variant = @inbounds c_variant[target][end]
                    spr_variant = @inbounds c_variant[spr][end]
                    recovfactor = target_status == recovered ? spr_recoveffect(recovday, targ_variant, spr_variant, infectset) : 1.0

                    vaxstatus = @inbounds c_vaxstatus[target]
                    vaxrcvd = @inbounds c_vaxrcvd[target][end]
                    vaxday = @inbounds c_vaxday[target][end]
                    vaxfactor = vaxstatus != :none ? spr_vaxeffect(infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday) : 1.0

                    spr_sickday = @inbounds c_sickday[spr]
                    targ_agegrp = @inbounds c_agegrp[target]

                risk = infectrisk(infectset, spr_variant, spr_sickday, targ_agegrp, recovfactor, vaxfactor)

                if isinfected(risk)
                    tovariant = @inbounds dovariant ? c_variant[spr][end] : :base
                    tocond = nil
                    tostatus = infectious
                    @inbounds set_infected!(target, c_cond, tocond, c_status, tostatus, c_sickday, 1, c_variant, tovariant)
                    # history_changes = ((tocond, 1), (tovariant, 1), (infectious, 1), (target_status, -1))
                end
            end  # if (touched ...)
        end  # if contactstatus
    end # for i = 1:nc

    return  # targ_agegrp, history_changes # n_contacts, n_touched, n_newly_infected
end       

function set_infected!(target, condcol, condval, statcol, statval, sickdaycol, sickdayval, variantcol, variantval)
    # wrapping args in some containers and deref'ing the containers will take too much time in the hottest loop of the simulation
    @inbounds begin
        condcol[target] = condval
        statcol[target] = statval
        sickdaycol[target] = sickdayval
        push!(variantcol[target], variantval)
    end
end


# simple make_sick! for a single person. Assumes that caller doesn't invoke structure of population data
function make_sick!(locdat, target::Int; cond, variant, sickday)
    @inbounds locdat.condition[target] = cond
    @inbounds locdat.status[target] = infectious
    push!(locdat.variant, variant)
    @inbounds locdat.sickday[target] = sickday
end

# complex make sick
function make_sick!(dat; cnt, ages, tocond, tovariant, tosickday=1) 

    @assert size(cnt, 1) == size(ages, 1) "size(cnt, 1) = $(size(cnt,1)) not equal size(ages, 1) = $(size(ages,1))"

    filt_unexp = optfindall(==(unexposed), dat.status, 1) # must be unexposed

    for i in 1:size(ages, 1)  # by target age groups

        filt_age = dat.agegrp[filt_unexp] .== ages[i] # age of the unexposed
        rowrange = 1:cnt[i]
        do_filt = filt_unexp[filt_age][rowrange]

        if size(do_filt, 1) == 0
            continue
        end

        dat.status[do_filt] .= infectious
        dat.cond[do_filt] .= tocond
        dat.sickday[do_filt] .= tosickday
        push!.(dat.variant[do_filt], tovariant)

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
model selects a contact_factor, c_f, based on age and infectious case
model selects a touch_factor, t_f, based on age and condition (includes unexposed and recovered)
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
