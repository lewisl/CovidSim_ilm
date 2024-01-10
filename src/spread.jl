################################
# spread.jl for ilm model
#    social distancing cases
#    spreading the infection
################################

"""
Hold factors that characterize how infectious the disease is.
Loaded by function setup from YAML parameter file.
"""
Base.@kwdef struct Infectparams
    sendrisk::Vector{Float64}
    recvrisk::Vector{Float64}
    recovery_immunity::Dict{Symbol, Float64}
    immunehalflife::Int64
    basemultiplier::Float64
end

        """
        Method for converting a dict loaded from YAML to this struct.
        Derefing a small struct is much faster than derefing a dict.
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

"""
Hold factors that describe social characteristics affecting spread of the disease.
Loaded by function setup from YAML parameter file.
"""
Base.@kwdef struct SocialParams
    gammashape::Float64
    indoor_uplift::Float64
    contactfactors::Matrix{Float64}     
    touchfactors::Matrix{Float64}     
end


"""
Hold parameters for social distancing cases used by callback function caserunner, below in file spread.jl.
"""
Base.@kwdef struct SpreadCase       # Base.@kwdef -> use keyword arguments and defaults in constructor
    name::Symbol
    day::Int
    cfdelta::Tuple{Float64,Float64}  
    tfdelta::Tuple{Float64,Float64}  
    comply::Float64             # compliance fraction
    cfcase::Matrix{Float64}
    tfcase::Matrix{Float64}
end


"""
Infectious people spread the virus to susceptible people for a single locale on thisday. 
Changes attribute columns in the population table. Runs social distancing cases.
"""
@inline function spread!(locdat, spr::Int, thisday::Int, sdcases, socialparams, 
     infectset, dovax, vaxset, density_factor, indoor_seq, poprange)

    spreader = locdat[spr]  # row of traits of the spreader person: a spreader "object"
     
    # retrieve params
    contactfactors = socialparams.contactfactors
    touchfactors   = socialparams.touchfactors
    gammashape     = socialparams.gammashape
    indoor_factor  = indoor_seq[thisday]

    @inbounds contact_param = spreader.sdcase === :none ? contactfactors : sdcases[spreader.sdcase]

    # how many contacts does the infected person have?
    contacts = get_contacts(spreader, poprange, density_factor, indoor_factor, gammashape, contact_param)

    for c in contacts
        contact = locdat[c] # row of traits of the contact person: a contact "object" 
        touch_param = contact.sdcase === :none ? touchfactors : sdcases[contact.sdcase].tfcase

        # is this a meaningful interaction?
        if istouched(contact, touch_param, indoor_factor)
            # will this contact get infected?
            if isinfected(contact, spreader, vaxset, dovax, infectset, thisday)
                @inbounds spr_variant = isempty(spreader.variant)  ? 0 : spreader.variant[end]
                make_sick!(contact, thisday, :nil, spr_variant)
            end
        end
    end
end       


###################################################################
# basic functions for the default definition of spread
###################################################################

# make_sick! for a single person. Assumes that caller doesn't invoke structure of population data
@inline function make_sick!(contact, thisday, cond, variant)  
    push!(contact.variant, variant)  # first of possibly several infections...  
    push!(contact.sickday, thisday)
    contact.duration = 1
    contact.cond = cond
    contact.status = :infectious
end


# make_sick! for multiple people
function make_sick!(dat; cnt, ages, tocond, tovariant, toduration=1) 

    @assert size(cnt, 1) == size(ages, 1) "size(cnt, 1) = $(size(cnt,1)) not equal size(ages, 1) = $(size(ages,1))"

    filt_unexp = optfindall(==(unexposed), dat.status, 1) # must be unexposed

    @inbounds for i in eachindex(ages)  # by target age groups

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


"""
Callback function returnned by function caserunner for social distancing cases.
"""
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
        filter1 = findall(((locdat.status .=== :unexposed) .| (locdat.status .=== :recovered)) .| 
                ((locdat.cond .=== :nil) .| (locdat.cond .=== :mild)))
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
