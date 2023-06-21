#######################################################################################
#  disease modeling helper functions used in spread! and progression!
#######################################################################################

"""
    vaxeffect(today, infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday; 
                mode=:spread, csig=6.0, decay_lower=0.15)

Immunity from vaccination for a single person.  Argument mode can be :spread or :progression to use
in either spread! or progression!, respectively.
"""
@inline @fastmath function vaxeffect(today, infectset, vaxset, vaxstatus, spr_variant, vaxrcvd, vaxday;
    mode=:spread, csig=6.0, decay_lower=0.15)

    # vaccine characteristics
    @inbounds begin
        vs = vaxset[vaxrcvd]
        halflife = vs.halflife
        vaxeffect = vs.effectiveness[vaxstatus][spr_variant]
        mineff = vs.day1_effect
        full_effect_days = vs.full_effect_days

        if mode == :spread
            infectfactor = vaxset[vaxrcvd].infectfactor[spr_variant]
        elseif mode == :progression
            infectfactor = 1.0
        else
            throw(DomainError(mode, "Argument must be :spread or :progression"))
        end
    end

    # person's vaccine conditions
    days_after_vax = max(today - vaxday, 0)
    days_after_full_effect = max(days_after_vax - full_effect_days, 0)     #clamp(today - (lastshotday + full_effect_days), 0, Int)

    rise = effect_rise(days_after_vax; mineff=mineff, delay_days=full_effect_days)
    decay = sigdecay(days_after_full_effect, halflife, csig=csig, decay_lower=decay_lower)     #   lindecay(days_after_full_effect, halflife, decay_lower)
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


function vax_recov(vaxfactor, recovfactor)::Float64
    min(vaxfactor, recovfactor)  # each factor is 1 - immunity_effect: small is good because risk = infectrisk * combined factor
end

sigmoid(x) = 1.0 / (1.0 + exp(-x))  # smoosh input to 0.0, 1.0--> not used at this time 6/12/2023


spreadin(risk, range) = range * risk - (0.5 * range)  # not used as of 6/12/2023
