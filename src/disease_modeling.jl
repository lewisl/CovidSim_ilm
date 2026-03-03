#######################################################################################
#  disease modeling helper functions used in spread! and progression!
#######################################################################################


"""
    get_contacts(spreader, poprange, density_factor, indoor_factor, gammashape, contact_param)

Return a vector of contacts, which are indices to the population table for the current locale
"""
@inline function get_contacts(spreader, poprange, density_factor, indoor_factor, gammashape, contact_param)

    numcontacts = @inbounds @fastmath how_many_contacts(density_factor, indoor_factor, gammashape,
        spreader.agegrp, spreader.cond, contact_param)

    contacts = rand(poprange, numcontacts)
end


"""
    how_many_contacts(density_factor, gammashape, agegrp, cond, contactfactors)::Int

Return the number of contacts that someone spreading the disease will make on a day. This
method uses the default contactfactors for the current spreader.
"""
@inline function how_many_contacts(density_factor, indoor_factor, gammashape, spr_agegrp, spr_cond, contactfactors)::Int64
    # indoor_factor is in [1.0, 1.4]. greater than 1.0 increases scale factor for gamma distribution
    @inbounds @fastmath scale = density_factor * indoor_factor * contactfactors[mapcondition(spr_cond), mapagegrp(spr_agegrp)]
    @fastmath round(Int, rand(Gamma(gammashape, scale)))
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
    @fastmath round(Int, rand(Gamma(gammashape, scale)))
end

# TODO assuming touch only depends on recipient (the contact) may be BAD.
"""
    function istouched(agegrp, lookup, touchfactors)::Bool

Returns true if the contact made was significant to the recipient or false if not.
This assumes touch only depends on the recipient.    

Very tricky.  We are ignoring the 4 sick conditions because even if touched they can't get sick.
"""
@inline function istouched(contact, touch_param, indoor_factor)
    touched = if (contact.status == :unexposed) | (contact.status == :recovered)  # only conditions that can get infected   
        touchprob = (indoor_factor == 1.0 ? touch_param[maptouch(contact.status), mapagegrp(contact.agegrp)] :
                     # squash multiplicative factor to stay under 1.0
                     clamp(indoor_factor * touch_param[maptouch(contact.status), mapagegrp(contact.agegrp)], 0.0, 0.97)) # or tanh--much slower
        rand(Binomial(1, touchprob)) == 1
    else
        false
    end
end


"""
    isinfected(contact, spreader, vaxset, dovax, infectset, thisday)::Bool

Returns true if the spreader infected the contact or false if not. 
Considers partial immunity if contact has recovered from previous infection.
Considers vaccination status of the contact.
Considers the variant of the disease the spreader is carrying.
"""
@inline function isinfected(contact, spreader, vaxset, dovax, infectset, thisday)::Bool

    @inbounds spr_variant = isempty(spreader.variant) ? 0 : spreader.variant[end]

    # effect on transmission based on how long ago a previously infected contact got over the disease
    recovfactor = recoveffect(thisday, contact, spr_variant, infectset)

    # effect on transmission based on whether, when, and which vaccine contact received
    vaxfactor = dovax ? vaxeffect(thisday, contact, vaxset, infectfactor, spr_variant) : 1.0

    # binomial probability of the contact getting infected from the contact with this spreader
    risk = infectrisk(infectset, spr_variant, spreader.duration, contact.agegrp, recovfactor, vaxfactor)
    return @fastmath rand(Binomial(1, risk)) == 1
end


@inline @fastmath function infectrisk(infectset, spr_variant, spr_duration,
    targ_agegrp, recovfactor::Float64, vaxfactor::Float64)

    # spreader person characteristics
    sendrisk = @inbounds infectset[spr_variant].sendrisk[spr_duration]

    # target person characteristics
    recvrisk = @inbounds infectset[spr_variant].recvrisk[mapagegrp(targ_agegrp)]

    combinedfactor = recvrisk * sendrisk * vax_recov(vaxfactor, recovfactor)
    risk = clamp(combinedfactor, 0.0, 0.97)    # required because combinedfactor could exceed 1.0
end


"""
    vaxeffect(thisday, infectfactor, p_vax, vaxstatus, spr_variant, vaxday; 
                csig=6.0, decay_lower=0.15)

Immunity from vaccination for a single person. In function progression!, the vaxeffect changes the severity and/or duration of the 
disease. In function spread!, the vaxeffect reduces the likelihood of getting the disease.
"""
@inline @fastmath function vaxeffect(thisday, contact, vaxset, infectfactor, spr_variant;
    csig=6.0, decay_lower=0.15)   # TODO decide where these inputs come from and what values to use

    factor = 1.0 # default value
    vaxstatus = contact.vaxstatus

    if !(vaxstatus === :none)

        # vaccine characteristics shortcuts
        @inbounds begin
            vaxday = contact.vaxday[end]
            p_vax = vaxset[contact.vaxrcvd[end]]  # struct of characteristics of vaccine received by the contact
            infectfactor = p_vax.infectfactor[spr_variant]
            halflife = p_vax.halflife
            vaxeffect = p_vax.effectiveness[vaxstatus][spr_variant]
            mineff = p_vax.day1_effect
            full_effect_days = p_vax.full_effect_days
        end

        # person's vaccine conditions
        days_after_vax = max(thisday - vaxday, 0)
        days_after_full_effect = max(days_after_vax - full_effect_days, 0)    
        rise = effect_rise(days_after_vax; mineff=mineff, delay_days=full_effect_days)
        decay = sigdecay(days_after_full_effect, halflife, csig=csig, decay_lower=decay_lower)     #   lindecay(days_after_full_effect, halflife, decay_lower)
        time_mod = rise * decay

        factor = max(1.0 - (time_mod * vaxeffect * infectfactor), 0.0)
    end

    return factor
end



"""
    recoveffect(recovday, contact_varient, spr_variant, infectset)

Immunity from recovery for a single person. Also, affects progression through disease stages.
"""
@inline function recoveffect(thisday, contact, spr_variant, infectset; csig=6.0, decay_lower=0.15)::Float64

    factor = 1.0 # default return value

    if contact.status === :recovered
        @inbounds recovday = isempty(contact.recovday) ? 0 : contact.recovday[end]
        days_post_recov = thisday - recovday

        contact_varient = contact.variant[end]

        @inbounds if days_post_recov >= 0
            # get the max immunity for the variant that target recovered from against the variant of the spreader
            immstrength = infectset[contact_varient].recovery_immunity[spr_variant]

            # get the declined value
            immhalflife = infectset[contact_varient].immunehalflife

            # immdecline = lindecay(days_post_recov, immhalflife, decay_lower)
            decay = sigdecay(days_post_recov, immhalflife, csig=csig, decay_lower=decay_lower)
            rise = effect_rise(days_post_recov)
            time_mod = rise * decay

            factor = 1.0 - (time_mod * immstrength)
        end
    end

    return factor  # this isn't Julian, but it's more obvious to most people
end


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


function vax_recov(vaxfactor, recovfactor)::Float64
    min(vaxfactor, recovfactor)  # each factor is 1 - immunity_effect: small is good because risk = infectrisk * combined factor
end


sigmoid(x) = 1.0 / (1.0 + exp(-x))  # smoosh input to 0.0, 1.0--> not used at this time 6/12/2023


spreadin(risk, range) = range * risk - (0.5 * range)  # not used as of 6/12/2023
