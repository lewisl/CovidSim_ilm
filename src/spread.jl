################################
# spread.jl for ilm model
#    social distancing cases
#    spreading the infection
################################

"""
Infectious people spread the virus to susceptible people for a single locale on thisday. 
Changes attribute columns in the population table. Runs social distancing cases.
"""
@inline function spread!(spr::Int, thisday::Int, sdcases, socialparams, 
     infectset, vaxset, density_factor, indoor_seq, poprange,    
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
        c_vaxday
     )

    # initialize return value
    num_infected = 0
     
    # retrieve params
    contactfactors = socialparams.contactfactors
    touchfactors   = socialparams.touchfactors
    gammashape     = socialparams.gammashape
    indoor_factor  = indoor_seq[thisday]

    # columns to modify in make_sick function
    sickcolumns = tuple(c_cond, c_status, c_duration, c_variant, c_sickday)

    # how many contacts does the infected person have?
    @inbounds contact_param = c_sdcase[spr] === :none ? contactfactors : sdcases[c_sdcase[spr]]
    numcontacts = @inbounds @fastmath how_many_contacts(density_factor, indoor_factor, gammashape, 
                                                        c_agegrp[spr], c_cond[spr], contact_param) 
    
    contacts = rand(poprange, numcontacts)
    for contact in contacts
        
        contact_status = c_status[contact]

        # does this contact experience a meaningful touch by the spreader?
        touched =   if (contact_status == :unexposed) | (contact_status == :recovered)  # only conditions that can get infected   
                        touch_param = c_sdcase[contact] === :none ? touchfactors : sdcases[c_sdcase[contact]].tfcase
                        istouched(c_agegrp[contact], contact_status, indoor_factor, touch_param)   # returns true or false
                    else
                        false
                    end

        # will this contact get infected?
        if touched
            recovday = isempty(c_recovday[contact]) ? 0 : c_recovday[contact][end]
            spr_variant = isempty(c_variant[spr])  ? 0 : c_variant[spr][end]
    
            # effect on transmission based on how long ago a previously infected contact got over the disease
            recovfactor =   if c_status[contact] == :recovered
                                contact_variant = c_variant[contact][end]
                                recoveffect(thisday, recovday, contact_variant, spr_variant, infectset)
                            else 
                                1.0
                            end
    
            vaxstatus = c_vaxstatus[contact]
            # effect on transmission based on which vaccine the contact received, how many times, and how long ago
            vaxfactor = if vaxstatus === :none
                            1.0 
                        else
                            vaxrcvd = c_vaxrcvd[contact][end]
                            vaxday = c_vaxday[contact][end]
                            vaxeffect(thisday, infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday; mode=:spread)
                        end
    
            spr_duration = c_duration[spr]  # number of days spreader has been infected
            contact_agegrp = c_agegrp[contact]
    
            # binomial probability of the contact getting infected from the contact with this spreader
            risk = infectrisk(infectset, spr_variant, spr_duration, contact_agegrp, recovfactor, vaxfactor)

            if isinfected(risk)
                make_sick!(sickcolumns, contact, thisday, :nil, spr_variant)
                num_infected += 1
            end
        end

    end
    return num_infected
end       


###################################################################
# basic functions for the default definition of spread
###################################################################


# simple make_sick! for a single person. Assumes that caller doesn't invoke structure of population data
@inline function make_sick!(locdat, target::Int; cond, variant, duration)
    @inbounds locdat.condition[target] = cond
    @inbounds locdat.status[target] = :infectious
    @inbounds push!(c_sickday[target], DAY_CTR[:day])
    @inbounds push!(locdat.variant[target], variant)
    @inbounds locdat.duration[target] = duration
end

# make_sick! by column for a single person: documents what passing the sickcolumns tuple actually does
@inline function make_sick!(c_cond, c_status, c_duration, c_variant, c_sickday, target, thisday, cond, variant)
    push!(c_variant[target], variant)  # first of possibly several infections...  c_variant[spr][end]
    push!(c_sickday[target], thisday)
    c_duration[target] = 1
    c_cond[target] = cond
    c_status[target] = :infectious
end

# sickcolumnsle = tuple(c_cond, c_status, c_duration, c_variant, c_sickday)
@inline function make_sick!(sickcolumns, target, thisday, cond, variant)   # this whole thing appears to compile away so it's free
    make_sick!(sickcolumns[1], sickcolumns[2], sickcolumns[3], sickcolumns[4], sickcolumns[5], target, thisday, cond, variant)
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

        dat.status[do_filt] .= :infectious
        dat.cond[do_filt] .= tocond
        dat.duration[do_filt] .= toduration
        push!.(dat.variant[do_filt], tovariant)

    end
end


"""
    how_many_contacts(density_factor, gammashape, agegrp, cond, contactfactors)::Int

Returns the number of contacts that someone spreading the disease will make on a day. This
method uses the default contactfactors for the current spreader.
"""
@inline function how_many_contacts(density_factor, indoor_factor, gammashape, agegrp, cond, contactfactors)::Int64 
    # indoor_factor is in [1.0, 1.4]. greater than 1.0 increases scale factor for gamma distribution
    @inbounds @fastmath scale = density_factor * indoor_factor * contactfactors[mapcondition(cond), mapagegrp(agegrp)]
    @fastmath round(Int,rand(Gamma(gammashape, scale)))
end


"""
    how_many_contacts(density_factor, gammashape, agegrp, cond, acase::SpreadCase)::Int

Returns the number of contacts that someone spreading the disease will make on a day. This 
method uses the spreadcase applicable to the current spreader but with contactfactors set by
a spreadcase.
"""
@inline function how_many_contacts(density_factor, indoor_factor, gammashape, agegrp, cond, acase::SpreadCase)::Int64
    # indoor_factor is in [1.0, 1.4]. greater than 1.0 increases scale factor for gamma distribution
    @inbounds @fastmath scale = density_factor * indoor_factor * acase.cfcase[mapcondition(cond), mapagegrp(agegrp)]  
    @fastmath round(Int,rand(Gamma(gammashape, scale)))
end


"""
    function istouched(agegrp, lookup, touchfactors)::Bool

Returns true if the contact made was significant to the recipient or false if not.
First method uses the default touchfactors for the current recipient.
Second method uses the spreadcase for the recipient.
"""
@inline function istouched(agegrp, lookup, indoor_factor, touchfactors)::Bool
    touchprob = (   indoor_factor == 1.0 ? touchfactors[maptouch(lookup), mapagegrp(agegrp)] : 
                    # squash multiplicative factor to stay under 1.0
                    simpleclamp(indoor_factor * touchfactors[maptouch(lookup), mapagegrp(agegrp)]) # or tanh--much slower
                    )
    return @inbounds @fastmath rand(Binomial(1, touchprob)) == 1
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
#  effect of immunity from prior recovery and vaccination
#
#############################################################################

# decay functions for immunity for decline from 1.0 to lower positive limit of function
# multiply times max immunity if less than 1.0        


"""
    effect_rise(days_since; mineff=0.65, delay_days=14)
  
Immunity effectiveness from vaccination or recovery ramps up.
Returns a value between mineff and 1.0. Linear increase.
"""
@inline @fastmath function effect_rise(days_since; mineff=0.65, delay_days=14)::Float64
    if days_since >= delay_days
        1.0
    else
        mineff + (days_since/delay_days * (1.0 - mineff))
    end
end


# gradual decay of vaccine effectiveness based on assumed half-life

@inline @fastmath function lindecay(t, h, lower)::Float64
    y = 0.5 ./ -h * t  + 1.0
    y = y < lower ? lower : y
end

expdecay(t,h) = exp(-(log(2)/h) * t)  

sigdecay(t, h; csig=5.0, decay_lower=0.1) = max(1.0 / (1.0 + exp.((t - h)/(t / csig + (h / csig)))), decay_lower)    

tbrk(h, lower) = 2.0 * h - (2.0 * h * lower)

intercept(t, hl, lower) = -0.3 * t / hl + 1.0

function lindecay2(t,hl,lower1, lower2)::Float64
    f1 = -t * 0.5 / hl + 1.0
    if  f1 >= lower1
        f1
    else
        clamp(-t * 0.2 / hl + intercept(tbrk(hl, lower1), hl, lower1), lower2, 1.0)
    end
end


function lindecayarr(t::AbstractVector{T} where T, hl, lower1, lower2)::Float64
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


function sigmoidshift(x; risk_discount=0.2)::Float64
    sigmoid(
            shifter(
                    clamp(x, 0.0, 1.0 + risk_discount),
                    0.0, 1.0, -4.0, 4.0
                    )
            )
end

@inline function simpleclamp(x; bot=0.0, top=0.97)::Float64
    clamp(x, bot, top)
end



function vax_recov_1(vaxfactor, recovfactor)::Float64
    x = vaxfactor * recovfactor
    x * exp(0.2 - x)
end

function vax_recov_2(vaxfactor, recovfactor)::Float64
    min(vaxfactor, recovfactor)  # each factor is 1 - immunity_effect: small is good because risk = infectrisk * combined factor
end


# choice of simple factor adjustments
vax_recov_combo = vax_recov_2  # this seems to cause type instability at line 355 or 356

# the squashfunc must keep the product of ALL combinations of sendrisk and recvrisk between 0.0 and 1.0 inclusive
squashfunc = simpleclamp   # This also causes type instability at line 358


"""
    spr_vaxeffect(today, infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday; csig=6.0, decay_lower=0.15)

Immunity from vaccination for a single person.
"""
@inline @fastmath function vaxeffect(today, infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday; 
    mode=:spread, csig=6.0, decay_lower=0.15)::Float64

    # vaccine characteristics
    @inbounds begin
    vs               = vaxset[vaxrcvd]
    halflife         = vs.halflife
    vaxeffect        = vs.effectiveness[vaxstatus][spr_variant]
    mineff           = vs.day1_effect
    full_effect_days = vs.full_effect_days

        if mode == :spread
            infectfactor     = vaxset[vaxrcvd].infectfactor[spr_variant]
        elseif mode == :progression
            infectfactor     = 1.0
        else
            throw(DomainError(mode, "Argument must be :spread or :progression"))
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
@inline function recoveffect(today, recovday, targ_variant, spr_variant, infectset; csig=6.0, decay_lower=0.15)::Float64

        days_post_recov = today - recovday 

        @inbounds if days_post_recov >= 0   
            # get the max immunity for the variant that target recovered from against the variant of the spreader
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



@inline @fastmath function infectrisk(infectset, spr_variant, spr_duration, 
    targ_agegrp, recovfactor::Float64, vaxfactor::Float64)

    # spreader person characteristics
    sendrisk = @inbounds infectset[spr_variant].sendrisk[spr_duration]

    # target person characteristics
    recvrisk = @inbounds infectset[spr_variant].recvrisk[mapagegrp(targ_agegrp)]

    combinedfactor = recvrisk * sendrisk * vax_recov_2(vaxfactor, recovfactor)
    risk = simpleclamp(combinedfactor)                 # this is required because combinedfactor could exceed 1.0
end


# impact of differential infectiousness for vaccine or variant
spreadin(risk) = 8.0 * risk - 4.0  # rescale risk to -4.0, 4.0
sigmoid(x) = 1.0 / (1.0 + exp(-x))  # smoosh input to 0.0, 1.0

altrisk(risk) = sigmoid(spreadin(risk))



#######################################################################
# social distancing cases
#      struct SpreadCase holds parameters for the case
#      implement the case: 
#           - set social distance compliance for each person
#           - define the contactfactors and touchfactors for the case
#######################################################################          
# mod_90 = sd_gen(start=90,cf=(.2,1.5), tf=(.18,.6),comply=.85)
# str_45 = sd_gen(start=45, comply=.90, cf=(.2,1.0), tf=(.18,.3))
# str_55 = sd_gen(start=55, comply=.95, cf=(.2,1.0), tf=(.18,.3))




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

    if startday == DAY_CTR[:day]

        if comply == 0.0  # magic signal: if comply is zero turn off this case for include_ages
            cancel_sd_case!(locdat, sdcases, name, include_ages, age_idx_loc)
            return
        end

        # create the SpreadCase in sdcases
        sdcases[name] = SpreadCase(
                            name    = name,   # TODO  if we never use this get rid of it
                            day     = startday,   
                            cfdelta = cf,         
                            tfdelta = tf,         
                            comply  = comply,     
                            cfcase  = shifter(socialparams.contactfactors, cf...),  
                            tfcase  = shifter(socialparams.touchfactors, tf...)     
                            )

        # load the sdcase column of the population table
        # filter1 is everyone who is unexposed, recovered or sick: nil or mild
        filter1 = findall(((locdat.status .== unexposed) .| (locdat.status .== recovered)) .| 
                ((locdat.cond .== nil) .| (locdat.cond .== mild)))
        if (comply == 1.0)   # include everyone in filter1 in this case
            complyfilter = filter1
        else
            complyfilter = sample(filter1, round(Int, comply*length(filter1)), replace=false)
        end
        
        if isempty(include_ages)   # include all include_ages
            locdat.sdcase[complyfilter] .= name
        else
            byage_idx = intersect(complyfilter, union((age_idx_loc[i] for i in include_ages)...))
            locdat.sdcase[byage_idx] .= name
        end
    end
end


function cancel_sd_case!(locdat, sdcases, name, include_ages, age_idx_loc)
    # filter on who is in this case now
    incase_idx = findall(locdat.sdcase .== name)

    if isempty(include_ages)   # include all ages
        locdat.sdcase[incase_idx] .= :none
        delete!(sdcases, name)  # there is no one left in this case...
    else  # only turn it off for some ages
        byage_idx = intersect(incase_idx, union((age_idx_loc[i] for i in include_ages)...))
        locdat.sdcase[byage_idx] .= :none
    end    

end
