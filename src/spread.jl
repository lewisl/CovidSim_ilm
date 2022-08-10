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
    function caserunner(locdat, socialparams, infectset, sdcases, age_idx_loc; day, startofday, locale)   
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
    how_many_contacts(density_factor, gammashape, agegrp, cond, contactfactors)::Int

Returns the number of contacts that someone spreading the disease will make on a day. This
method uses the default contactfactors for the current spreader.
"""
@inline function how_many_contacts(density_factor, gammashape, agegrp, cond, contactfactors)::Int 
    @inbounds @fastmath scale = density_factor * contactfactors[mapcondition(cond), mapagegrp(agegrp)]
    @fastmath round(Int,rand(Gamma(gammashape, scale)))
end

"""
    how_many_contacts(density_factor, gammashape, agegrp, cond, acase::Spreadcase)::Int

Returns the number of contacts that someone spreading the disease will make on a day. This 
method uses the spreadcase applicable to the current spreader but with contactfactors set by
a spreadcase.
"""
@inline function how_many_contacts(density_factor, gammashape, agegrp, cond, acase::Spreadcase)::Int
    @inbounds @fastmath scale = density_factor * acase.cfcase[mapcondition(cond), mapagegrp(agegrp)]  
    @fastmath round(Int,rand(Gamma(gammashape, scale)))
end


"""
    function istouched(agegrp, lookup, touchfactors)::Bool

Returns true if the contact made was significant to the recipient or false if not.
First method uses the default touchfactors for the current recipient.
Second method uses the spreadcase for the recipient.
"""
@inline function istouched(agegrp, lookup, touchfactors)::Bool
    return @inbounds @fastmath rand(Binomial(1, touchfactors[maptouch(lookup), mapagegrp(agegrp)])) == 1
end



"""
    function isinfected(infectparams, spreaderduration, contactagegrp)::Bool

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
    effect_rise(days_since; mineff=0.65, delay_days=14)
  
Immunity effectiveness from vaccination or recovery ramps up.
Returns a value between mineff and 1.0. Linear increase.
"""
@inline @fastmath function effect_rise(days_since; mineff=0.65, delay_days=14)
    if days_since >= delay_days
        1.0
    else
        mineff + (days_since/delay_days * (1.0 - mineff))
    end
end


# gradual decay of vaccine effectiveness based on assumed half-life

@inline @fastmath function lindecay(t, h, lower)
    y = 0.5 ./ -h * t  + 1.0
    y = y < lower ? lower : y
end

expdecay(t,h) = exp(-(log(2)/h) * t)  

sigdecay(t, h; csig=5.0, decay_lower=0.1) = max(1.0 / (1.0 + exp.((t - h)/(t / csig + (h / csig)))), decay_lower)    

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
    @inbounds for i = eachindex(arr)
        f1 = -t[i] * 0.5 / hl + 1.0
        if  f1 >= lower1
            arr[i] = f1
        else
            arr[i] = clamp(-t[i] * 0.2 / hl + icept, lower2, 1.0)
        end
    end
    return arr
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



function vax_recov1(vaxfactor, recovfactor)
    x = vaxfactor * recovfactor
    x * exp(0.2 - x)
end

function vax_recov2(vaxfactor, recovfactor)
    min(vaxfactor, recovfactor)  # each factor is 1 - immunity_effect: small is good because risk = infectrisk * combined factor
end


# choice of simple factor adjustments
vax_recov_combo = vax_recov2

# the squashfunc must keep the product of ALL combinations of sendrisk and recvrisk between 0.0 and 1.0 inclusive
squashfunc = simpleclamp


"""
    spr_vaxeffect(today, infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday; csig=6.0, decay_lower=0.15)

Immunity from vaccination for a single person.
"""
@inline @fastmath function vaxeffect(today, infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday; mode=:spread, csig=6.0, decay_lower=0.15)

    # vaccine characteristics
    @inbounds begin
    vs               = vaxset[vaxrcvd]
    halflife         = vs.halflife
    vaxeffect        = vs.effectiveness[vaxstatus][spr_variant]
    mineff           = vs.day1_effect
    full_effect_days = vs.full_effect_days

        if mode == :spread
            infectfactor     = vaxset[vaxrcvd].infectfactor[spr_variant]
        elseif mode == :transition
            infectfactor     = 1.0
        else
            throw(DomainError(mode, "Argument must be :spread or :transition"))
        end
    end

    # person's vaccine conditions
    days_after_vax = max(today - vaxday, 0)
    days_after_full_effect = max(days_after_vax - full_effect_days, 0)     #clamp(today - (lastshotday + full_effect_days), 0, Int)

    rise = effect_rise(days_after_vax; mineff=mineff, delay_days=full_effect_days)
    decay =  sigdecay(days_after_full_effect, halflife, csig=csig, decay_lower=decay_lower)     #   lindecay(days_after_full_effect, halflife, decay_lower)
    time_mod = rise * decay

    factor = max(1.0 - (time_mod * vaxeffect * infectfactor), 0.0)

    return factor
end


"""
    recoveffect(recovday, targ_variant, spr_variant, infectset)

Immunity from recovery for a single person.
"""
@inline function recoveffect(today, recovday, targ_variant, spr_variant, infectset; csig=6.0, decay_lower=0.15)

        days_post_recov = today - recovday 

        @inbounds if days_post_recov >= 0   # TODO should be an assert: does this run day of or day after recovery?
            # get the max immunity
            immstrength = infectset[targ_variant].recovery_immunity[spr_variant]

            # get the declined value
            immhalflife = infectset[targ_variant].immunehalflife

            # immdecline = lindecay(days_post_recov, immhalflife, decay_lower)
            decay = sigdecay(days_post_recov, immhalflife, csig=csig, decay_lower=decay_lower)
            rise = effect_rise(days_post_recov)
            time_mod = rise * decay

            factor = 1.0 - (time_mod * immstrength)
        else
            factor = 1.0
        end

    return factor
end



@inline @fastmath function infectrisk(infectset, spr_variant, spr_duration, targ_agegrp, recovfactor, vaxfactor)

    # spreader person characteristics
    sendrisk = @inbounds infectset[spr_variant].sendrisk[spr_duration]

    # target person characteristics
    recvrisk = @inbounds infectset[spr_variant].recvrisk[Int(targ_agegrp)]

    combinedfactor = recvrisk * sendrisk * vax_recov_combo(vaxfactor, recovfactor)
    risk = squashfunc(combinedfactor)                 # this is required because combinedfactor could exceed 1.0
end


# impact of differential infectiousness for vaccine or variant
spreadin(risk) = 8.0 * risk - 4.0  # rescale risk to -4.0, 4.0
sigmoid(x) = 1.0 / (1.0 + exp(-x))  # smoosh input to 0.0, 1.0

altrisk(risk) = sigmoid(spreadin(risk))


"""
Infectious people spread the virus to susceptible people for a single locale. Changes attribute
columns in the population table. Runs social distancing cases.
"""
@inline function spread!(spr::Int, thisday::Int, sdcases, socialparams,   
     infectset, vaxset, density_factor, poprange,    
        c_cond,
        c_status,
        c_agegrp,
        c_duration,
        c_sdcomply,
        c_sickday,
        c_variant,
        c_vaxstatus,
        c_recovday,
        c_vaxrcvd,
        c_vaxday
     )

    today = thisday

    # retrieve params
    contactfactors = socialparams.contactfactors
    touchfactors   = socialparams.touchfactors
    gammashape     = socialparams.gammashape

    targets = social_model(spr, poprange, contactfactors, touchfactors, sdcases, density_factor, gammashape,
                           c_sdcomply, c_agegrp, c_cond, c_status)

    
    infection_model!(spr, targets, today, infectset, vaxset,  
                    c_recovday, c_variant, c_vaxstatus, c_vaxrcvd, c_vaxday,   
                    c_duration, c_agegrp, c_sickday, c_cond, c_status)        

    return nothing 
end       


"""
    Who has been touched by a spreader and might later become infected?
"""
function social_model(spr, poprange, contactfactors, touchfactors, sdcases, density_factor, gammashape,
                      c_sdcomply, c_agegrp, c_cond, c_status)

    @inbounds contact_param = c_sdcomply[spr] === :none ? contactfactors : sdcases[c_sdcomply[spr]]
    numcontacts = @inbounds @fastmath how_many_contacts(density_factor, gammashape, c_agegrp[spr], c_cond[spr], contact_param)  

    targets =  @fastmath [target for target in [rand(poprange) for i in 1:numcontacts] if (
                        @inbounds begin 
                            target_status = c_status[target]
                            if (target_status == unexposed) | (target_status == recovered)  # only conditions that can get infected   
                                touch_param = c_sdcomply[target] === :none ? touchfactors : sdcases[c_sdcomply[target].tfcase]
                                istouched(c_agegrp[target], target_status, touch_param)   # is the contact significant? returns true or false
                            else
                                false
                            end
                        end
                        )
                ]
    return targets
end


"""
    Of the targets who have been touched by a spreader, update trait columns for those who become infected...
"""
function infection_model!(spr, targets, today, infectset, vaxset,  
            c_recovday, c_variant, c_vaxstatus, c_vaxrcvd, c_vaxday,   # trait columns
            c_duration, c_agegrp, c_sickday, c_cond, c_status)         # trait columns

    for target in targets
        recovday = c_recovday[target][end]
        spr_variant = c_variant[spr][end]

        recovfactor = if c_status[target] == recovered
                            targ_variant = c_variant[target][end]
                            recoveffect(today, recovday, targ_variant, spr_variant, infectset)
                        else 
                            1.0
                        end

        vaxstatus = c_vaxstatus[target]
        vaxfactor = if vaxstatus === :none
                        1.0 
                    else
                        vaxrcvd = c_vaxrcvd[target][end]
                        vaxday = c_vaxday[target][end]
                        vaxeffect(today, infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday; mode=:spread)
                    end

        spr_duration = c_duration[spr]
        targ_agegrp = c_agegrp[target]

        risk = infectrisk(infectset, spr_variant, spr_duration, targ_agegrp, recovfactor, vaxfactor)

        if isinfected(risk)
            push!(c_variant[target], spr_variant)  # first of possibly several infections...  c_variant[spr][end]
            push!(c_sickday[target], today)
            c_duration[target] = 1
            c_cond[target] = nil
            c_status[target] = infectious
        end
    end
end


# simple make_sick! for a single person. Assumes that caller doesn't invoke structure of population data
function make_sick!(locdat, target::Int; cond, variant, duration)
    @inbounds locdat.condition[target] = cond
    @inbounds locdat.status[target] = infectious
    @inbounds push!(c_sickday[target], day_ctr[:day])
    push!(locdat.variant, variant)
    @inbounds locdat.duration[target] = duration
end

# complex make sick
function make_sick!(dat; cnt, ages, tocond, tovariant, toduration=1) 

    @assert size(cnt, 1) == size(ages, 1) "size(cnt, 1) = $(size(cnt,1)) not equal size(ages, 1) = $(size(ages,1))"

    filt_unexp = optfindall(==(unexposed), dat.status, 1) # must be unexposed

    @inbounds for i in 1:size(ages, 1)  # by target age groups

        filt_age = dat.agegrp[filt_unexp] .== ages[i] # age of the unexposed
        rowrange = 1:cnt[i]
        do_filt = filt_unexp[filt_age][rowrange]

        if size(do_filt, 1) == 0
            continue
        end

        dat.status[do_filt] .= infectious
        dat.cond[do_filt] .= tocond
        dat.duration[do_filt] .= toduration
        push!.(dat.variant[do_filt], tovariant)

    end
end


