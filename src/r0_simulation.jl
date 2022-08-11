
# spread!(spr::Int, thisday::Int, sdcases, socialparams,   
# infectset, vaxset, density_factor, poprange,    
#    c_cond,
#    c_status,
#    c_agegrp,
#    c_duration,
#    c_sdcomply,
#    c_sickday,
#    c_variant,
#    c_vaxstatus,
#    c_recovday,
#    c_vaxrcvd,
#    c_vaxday
# )

"""
r0_sim(; pop=200_000, age_dist=AGE_DIST, dectree=dectree, socialparams=socialparams, infectparams=infectparams, density_factor=1.0, scale=5)
r0_sim(locdat; age_dist=AGE_DIST, dectree=dectree, socialparams=socialparams, infectparams=infectparams, sdcases=sdcases, density_factor=1.0, scale=5)

Simulates r0 or rt. The first method creates a population and tracks how many infections
are caused by first generation spreaders and NOT spreaders who were infected by the
first generation. The simulates r0

The second method simulates r at time t given the characteristics of the simulation
you are running. This shows how r, reproduction rate, is affected by public health
measures and the characteristics of the population over time. This simulates r(t).
"""
function r0_sim(; pop=200_000, age_dist=AGE_DIST, progressionset, infectset, vaxset, socialparams, density_factor=1.0, scale=3)
    if pop < 200_000
        @warn "Population size should be greater than 200,000 for r0_sim. Proceeding with 200,000."
        pop = 200_000
    end
                

    # create simulation population
    r0pop = pop_data(pop, age_dist=AGE_DIST )

    # seed spreaders in each age group proportional to age distribution
    cnt_by_agedist = round.(Int, age_dist ./ minimum(age_dist))
    scale = set_by_level(count(r0pop.status .!= dead))
    cnt_by_agedist .*= scale # update with scale

    for i in AGEGRPS
        idx = findall(r0pop.agegrp .== i) 

        for j = 1:cnt_by_agedist[Int(i)]
            r0pop.status[idx] = infectious
            r0pop.cond[idx] = nil
            r0pop.duration[idx] = 1
            idx += 1
        end
    end

    # set infect_idx based on seeding: never update so we measure only 1st gen. spreaders
    gen1_infect_idx = findall(r0pop.status .== infectious)
    gen1_infected = length(gen1_infect_idx)
    r0_infected = 0

    for i = 1:DURATIONLIM        
        contactable_idx = findall(r0pop.status .!= dead)

        # spread with only spreaders from the gen1 infected pool
        n_newly_infected = spread!(r0pop, gen1_infect_idx, contactable_idx,  sdcases, socialparams, infectparams, density_factor)  
        r0_infected += n_newly_infected

        # progression all who are currently infected
        all_infect_idx = findall(r0pop.status .== infectious)
        progression!(r0pop, all_infect_idx, progressionset) 

        # of the gen1 infected, who is still infected? (some will have progressioned to recovered or dead)
        gen1_infect_idx = filter(x -> r0pop.status[x] == infectious, gen1_infect_idx)
    end

    r0 =  r0_infected / gen1_infected   
    return r0

end

# progression!(p, infectset, progressionset, vaxset, dovax, riskshift!, transvec,
#        c_cond,
#        c_status,
#        c_agegrp,
#        c_duration,
#        c_sdcomply,
#        c_variant,
#        c_vaxstatus,
#        c_recovday,
#        c_vaxrcvd,
#        c_vaxday,
#        c_deadday
#        )


function r0_sim(locdat; age_dist=AGE_DIST, progressionset=progressionset, socialparams=socialparams, infectparams=infectparams, sdcases=sdcases, density_factor=1.0, scale=5)
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

    for i in AGEGRPS  # set the spreaders for the r0 simulation
    idx = findall((r0pop.agegrp .== i) .& (r0pop.status .== unexposed))
    for j = 1:age_relative[Int(i)]
        spr = idx[j]
        r0pop.status[spr] = infectious
        r0pop.cond[spr] = nil
        r0pop.duration[spr] = 1
    end
    end     

    r0_infected = 0 
    for i = 1:DURATIONLIM      
    infect_idx = findall((r0pop.status .== infectious) .& (r0pop.duration .> 0))
    contactable_idx = findall(r0pop.status .!= dead)
    r0_infected += spread!(r0pop, infect_idx, contactable_idx, sdcases, socialparams, infectparams, density_factor)  

    progression!(r0pop, infect_idx, progressionset, infectset) 

    progression!(p, infectset, progressionset, vaxset, dovax, riskshift!, transvec,
        c_cond,
        c_status,
        c_agegrp,
        c_duration,
        c_sdcomply,
        c_variant,
        c_vaxstatus,
        c_recovday,
        c_vaxrcvd,
        c_vaxday,
        c_deadday
        )



    # eliminate the new spreaders so we only track the original spreaders
    newsick_idx = findall(r0pop.duration .== 1)

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
