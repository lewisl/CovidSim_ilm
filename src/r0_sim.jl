function r0_sim(pop=200_000, age_dist=AGE_DIST, dectree, infectparams, density_factor=1.0; scale=5)
    # create simulation population
    r0pop = pop_data(pop, age_dist)

    # seed spreaders in each age group proportional to age distribution
    cnt_by_agedist = round.(Int, age_dist ./ minimum(age_dist))
    cnt_by_agedist .*= scale # update with scale
    for i in agegrps
        idx = findfirst(x->x==i, r0pop.agegrp)
        for j = 1:cnt_by_agedist[i]
            r0pop.status[idx] = infectious
            r0pop.cond[idx] = nil
            r0pop.duration[idx] = 1
            idx += 1
        end
    end

    # set infect_idx based on seeding: never update so we measure only 1st gen. spreaders
    gen1_infect_idx = findall(r0pop .== infectious)
    gen1_infected = length(gen1_infect_idx)
    sdcases = []   # TODO this should be an input based on current context of simulation
    r0_infected = 0

    for i = 1:DURATIONLIM        
        contactable_idx = findall(locdat.status .!= dead)
        n_newly_infected = spread!(r0pop, gen1_infect_idx, contactable_idx,  sdcases, socialparams, infectparams, density_factor)  
        infect_idx = findall(locdat.status .== infectious)
        r0_infected += n_newly_infected
        transition!(r0pop, infect_idx, dectree) 
    end

    r0 =  r0_infected / gen1_infected   # n_newly_infected / cnt_spreaders
    return r0

end