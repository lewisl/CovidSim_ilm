# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Julia 1.5.2
#     language: julia
#     name: julia-1.5
# ---

# %% [markdown]
# # Prototype Spread and Transition with ILM Approach

# %%
using CovidSim_ilm

# %%
using StatsBase
using DelimitedFiles
using Distributions
using YAML

# %%
ilmat = readdlm("../data/ilmtestdata.csv", ',',Int, header=true)[1]
ilmat = repeat(ilmat, 10_000)
refresh = copy(ilmat)

# %% [markdown]
# ## Spread

# %%
infectfilename="../parameters/infectparams.yml"
socialfilename="../parameters/socialparams.yml"
infectparams = CovidSim.read_infectparams(infectfilename)
socialparams = CovidSim.read_socialparams(socialfilename)
contact_factors = socialparams[:contact_factors]
touch_factors = socialparams[:touch_factors]
send_risk = infectparams[:send_risk]
recvrisk = infectparams[:recv_risk]
riskmx = CovidSim.send_risk_by_recv_risk(send_risk, recv_risk) 

# %%
filt_spread_bit = ilmat[:,1] .== 2
filt_spread_idx = findall(filt_spread_bit)
@show size(filt_spread_idx,1) == count(filt_spread_bit)
size(filt_spread_idx,1)

# %%
# for this test case lets have fewer people who are infected
# resetidx = sample(findall(spreadidx),72000, replace=false)
# ilmat[resetidx,1] .= 1
# index of all spreaders
filt_spread_idx = findall(ilmat[:, cpop_status] .== 2)
n_spreaders = size(filt_spread_idx, 1);

# %%
# no of contacts for each spreader
time_it = @elapsed begin
    n_spreaders = size(filt_spread_idx, 1);
#     @show n_spreaders

    contacts_per_spreader = zeros(Int, size(filt_spread_idx,1),2) # second column for duration of the spreader
    for i in 1:size(contacts_per_spreader, 1)  # cond, agegrp
        scale = contact_factors[ilmat[filt_spread_idx[i], cpop_cond]-4, ilmat[filt_spread_idx[i], cpop_agegrp]]
        contacts_per_spreader[i, 1] = round(Int,rand(Gamma(1.3,scale))) # assume density_factor = 1.0
        contacts_per_spreader[i, 2] = ilmat[filt_spread_idx[i], cpop_duration]   # duration of the spreader who made this contact
    end
    n_contacts = sum(contacts_per_spreader[:,1])
#     @show n_contacts

    # assign the contacts 
    filt_contactable_idx = findall(ilmat[:, cpop_status] .!= dead)
    n_contactable = size(filt_contactable_idx, 1)
#     @show n_contactable
    choose_contacts = sample(filt_contactable_idx, min(n_contacts, n_contactable), replace=false)
    n_contacts = size(choose_contacts,1)
#     @show n_contacts, "2nd time"

    # determine contacts are consequential touches and if newly infected
    n_touched = 0
    n_newly_infected = 0

    contact_ptr = 0
    for spr in 1:n_spreaders
        nc = contacts_per_spreader[spr,1]
        duration_spr = contacts_per_spreader[spr, 2]
        # who are the contacts for this spreader?
            contact_selector = (contact_ptr+1):(contact_ptr+nc)
            contact_ptr += nc
        contact_persons = choose_contacts[contact_selector]
#         spr > 50 && break
        for contact_person in contact_persons
#             @show spr,contact_person

            status = ilmat[contact_person, cpop_status]
            characteristic =  status in [1,3] ? [1,0,2][status] : ilmat[contact_person, cpop_cond]-2 # max(0,ilmat[person, cpop_cond]-2
            @debug characteristic < 1 && error("bad characteristic value")
            agegrp = ilmat[contact_person, cpop_agegrp]
            touched = rand(Binomial(1, touch_factors[characteristic, agegrp]))
            n_touched += touched
            if touched == 1 && characteristic == unexposed
                prob = riskmx[duration_spr, agegrp]
                newly_infected = rand(Binomial(1, prob))
                if newly_infected == 1
                    ilmat[contact_person, cpop_cond] = nil # nil === asymptomatic or pre-symptomatic
                    ilmat[contact_person, cpop_status] = infectious
                end
                n_newly_infected += newly_infected
            end
        end
    end
    
end

@show time_it
@show n_touched
@show n_newly_infected

# %%
spreadidx = findall(ilmat[:, cpop_status] .== 2)
size(spreadidx,1)

# %%
ilmat

# %% [markdown]
# ## Transition

# %%
std_file = "../parameters/dec_tree_all_25.yml"
treedict, decpoints = setup_dt(std_file)
treedict

# %%
treedict[1][[25,7]]

# %%
function precalc_agegrp_filt(dat)
    agegrp_filt_bit = Dict(agegrp => dat[:, cpop_agegrp] .== agegrp for agegrp in agegrps)
    agegrp_filt_idx = Dict(agegrp => findall(agegrp_filt_bit[agegrp]) for agegrp in agegrps)
    return agegrp_filt_bit, agegrp_filt_idx
end
agegrp_filt_bit, agegrp_filt_idx = precalc_agegrp_filt(ilmat);

# %%
function transition!(dt, all_decpoints, locale, dat)  

    # @assert (length(locale) == 1 || typeof(locale) <: NamedTuple) "locale must be a single integer or NamedTuple"
    
    for agegrp in agegrps
        tree = dt[agegrp]
        for node in keys(tree)
            nodeduration, fromcond = node
            for branch in tree[node]["branches"]
                filt = ( (dat[agegrp_filt_idx[agegrp],cpop_cond] .== fromcond) .& 
                         (dat[agegrp_filt_idx[agegrp],cpop_status] .== infectious) .&
                         (dat[agegrp_filt_idx[agegrp], cpop_duration] .== nodeduration)  )
                
                    # boolean filt = ( (dat[:,cpop_cond] .== fromcond) .& 
                    #                  (dat[:,cpop_status] .== infectious) .&
                    #                  (dat[:, cpop_duration] .== nodeduration) .&
                    #                  (agegrp_filt_bit[agegrp])  )

                filt = agegrp_filt_idx[agegrp][filt]
                cnt_transition = size(filt, 1)  
                
                if cnt_transition > 0
                    
                    probs = tree[node]["probs"]
                    outcomes = tree[node]["outcomes"]
                    
                    choices = rand(Categorical(probs), cnt_transition) 
                                        
                    for (idx, person) in enumerate(filt)   # findall(filt)                       
                        new_stat_cond = outcomes[choices[idx]]
                        if new_stat_cond in (dead, recovered)  # change status
                            dat[person, cpop_status] = new_stat_cond
                        else   # change disease condition
                            dat[person, cpop_cond] = new_stat_cond
                        end  # if/else
                    end  # for (idx, person)
                    
                end  # if cnt_transition
            end  # for branch
        end  # for node
    end  #for agegrp

    # bump everyone who is still infectious all at once in one go
    filt = (dat[:,cpop_status] .== infectious)
    dat[filt, cpop_duration] .+= 1   
end  
    

transition!(treedict, decpoints, 38015, ilmat);   
ilmat


# %%
count((ilmat - refresh) .!= 0)

# %% [markdown]
# # Seed

# %%
time_it = @elapsed CovidSim.make_sick!(ilmat; cnt=[3,3], fromage=[2,3], tocond=nil)

#=
function seed!(day, cnt, duration, cond, agegrps, dat)  #  locale,
    if day == day_ctr[:day] ]
        
        make_sick!(dat; cnt, fromage, tocond, toduration=1)
        
    end 
end
=#
@show time_it

# %%
count((ilmat - refresh) .!= 0)

# %%
