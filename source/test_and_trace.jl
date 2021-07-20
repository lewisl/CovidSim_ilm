"""
Quar_Loc is a type alias for 

```NamedTuple{(:locale, :quar_date),Tuple{Int64,Int64}}```

example:  (locale=53038, quar_date=50)

Used as a locale in an isolatedmx 3-d array. Quar_Loc locale provides
the geographical locale and the simulation ordinal date that a cohort of infected folks
entered quarantine. Locale values are US Census FIPS codes at the county level.
"""
const Quar_Loc = NamedTuple{(:locale, :quar_date),Tuple{Int64,Int64}}

"""
Test_Loc is a type alias for 

```NamedTuple{(:locale, :test_date),Tuple{Int64,Int64}}```

example:  (locale=53038, test_date=50)

Used as a locale in an population matrix 3-d array for tracking tests conducted by locale. 
Provides the geographical locale and the simulation ordinal date that of the test.
Locale values are US Census FIPS codes at the county level.
"""
const Test_Loc = NamedTuple{(:locale, :test_date),Tuple{Int64,Int64}}

const map2access = (unexposed= 1, infectious=-1, recovered=2, dead=-1, 
                  nil= 3, mild=  4, sick= -1, severe= -1)

# cache values needed during test and track
const tnt_stash = Dict{NamedTuple, Array}() # works for Quar_loc or Test_loc--and they are different


"""
    function t_n_t_case_gen(start_day, end_day;         
        tc_perday=400, sensitivity=.90, specificity=0.90, infect_prior=0.5, 
        test_pct=.95,  test_delay=1, generations=3, qdays=15) 

Returns a function to use in runcases input to function ```run_a_sim```.
"""
function t_n_t_case_gen(start_day, end_day;         # these args go into the returned t_n_t_case
    tc_perday=400, sensitivity=.90, specificity=0.90, infect_prior=0.5, test_pct=.95,  # optional
    q_comply=0.8, c_comply=0.9, breakout_pct=.3, test_delay=3, generations=3, qdays=15,
    target_test=false, past_contacts=false) 
    # args match runcases loop in run_a_sim
    function runcase(locale; opendat, isodat, testdat, infectparams)  # case loop in run_a_sim provides required args
        t_n_t_case(start_day, end_day; 
                   socialparams=socialparams, infectparams=infectparams, opendat=opendat, isodat=isodat,  
                   testdat=testdat, locale=locale, tc_perday=tc_perday, sensitivity=sensitivity, specificity=specificity, 
                   infect_prior=infect_prior, test_pct=test_pct, q_comply=q_comply, c_comply=c_comply,
				   breakout_pct=breakout_pct, test_delay=test_delay, generations=generations, qdays=qdays,
                   target_test=target_test, past_contacts=past_contacts)
    end
end


function t_n_t_case(start_date, end_date; 
                socialparams, infectparams, opendat, isodat, testdat, locale,     
                tc_perday=1000, sensitivity=.95, specificity=0.90, infect_prior=0.5, test_pct=.95,  
                q_comply=0.8, c_comply=0.9, breakout_pct=.3, test_delay=3, generations=3, qdays=15,
                target_test=false, past_contacts=false)

    thisday = day_ctr[:day]
    ret_conds = [unexposed, recovered, nil, mild, sick, severe] 

    
    # Actions that can occur after or during the case:
        # do we need to UNquarantine anyone today?  
        for q in keys(isodat)
            if typeof(q) <: Quar_Loc # then k is a Quar_Loc = (locale=locale, quar_date=thisday)
                if q.locale == locale
                    # quarantine breakouts
                    if q.quar_date < thisday < q.quar_date + qdays
                        day_of_q = thisday - q.quar_date
                        # are there breakouts? 1 element of breakout array
                        cnt = get(get(tnt_stash, q, [Int(0)]), day_of_q, Int(0)) 

                        println("type of cnt: ", typeof(cnt))

                        if cnt > 0
                            # println("  GOT HERE:  unquarantine breakouts  $cnt ")
                            t_n_t_unquarantine(cnt_2_array(cnt, isodat[q]), q, opendat=opendat, 
                                               isodat=isodat, infectparams=infectparams)
                            push!(tntq, (day=thisday, breakout=cnt))
                        end
                    # quarantines ending today
                    elseif q.quar_date + qdays == thisday # end of quarantine is today
                        cnt = grab(ret_conds, agegrps, sickdays, q, isodat)
                        t_n_t_unquarantine(cnt, q, opendat, isodat, infectparams)
                        push!(tntq, (day=day_ctr[:day], unquarantine=sum(cnt)))
                        delete!(isodat, q)  # remove dated locale
                        delete!(tnt_stash, q)  # remove dated locale from stash
                    end
                end
            end
        end

        # are any delayed test results coming back today, which will lead to quarantines?
        if test_delay != 0
            for t in keys(tnt_stash)
                if typeof(t) <: Test_Loc
                    if t.locale == locale
                        if t.test_date + test_delay == thisday
                            # quarantine folks with positive test results
                            qloc = (locale=locale, quar_date=thisday)
                            put_in = tnt_stash[t]
                            if breakout_pct != 0.0  # future breakouts from this quarantine cohort
                                breakout!(breakout_pct, put_in, qloc, qdays) # Int[] (14, )
                            end
                            t_n_t_quarantine(put_in, qloc::Quar_Loc, opendat, 
                                             isodat, infectparams)
                            push!(tntq, (day=day_ctr[:day], quarantine=sum(put_in)))
                            delete!(tnt_stash, t) # pop the stash (delete!)
                        end
                    end
                end
            end
        end

    # do more testing today
    test_and_trace(start_date, end_date; 
        socialparams=socialparams, infectparams=infectparams, opendat=opendat, isodat=isodat,    # from case loop
        testdat=testdat, locale=locale, tc_perday=tc_perday, sensitivity=sensitivity, specificity=specificity, # optional
        infect_prior=infect_prior, test_pct=test_pct, q_comply=q_comply, c_comply=c_comply, 
		breakout_pct=breakout_pct, test_delay=test_delay, generations=generations, qdays=qdays,
        target_test=target_test, past_contacts=past_contacts)
end


function test_and_trace(start_date, end_date; 
    socialparams, infectparams, opendat, isodat, testdat, locale,  # from case loop
    tc_perday=1000, sensitivity=.95, specificity=0.95, infect_prior=0.5, test_pct=.95, # optional
    q_comply=0.8, c_comply=0.9, breakout_pct=.3, test_delay=3, generations=3, qdays=15,
    target_test=false, past_contacts=false)
    
    thisday = day_ctr[:day]
    test_conds = [unexposed, recovered, nil, mild]

    if thisday  == 1    # TODO when is the right time?  what is the right cleanup?
        cleanup_stash(tnt_stash)
    end


    if start_date <= thisday < end_date

        # who gets tested?  
        avail_to_test = zeros(Int, 25,4,5)
        if target_test
            sel_age = get_next!(nxt) # circular cycle through 1:4
            if sel_age == 4; sel_age = 4:5; end  # only sampled age group is non-zero
            avail_to_test[:,:,sel_age] = grab(test_conds, sel_age, sickdays, locale, opendat)
        else
            avail_to_test[:] = grab(test_conds, agegrps, sickdays, locale, opendat)
        end

        if sum(avail_to_test) == 0
            println("on $thisday no one to test")
            return
        end

        n_conds = length(map2access[unexposed]:map2access[mild])
        to_test = Dict(i=>zeros(Int, sickdaylim, n_conds, agegrps) for i in 1:generations) # TODO pre-allocate?
        postests = Dict(i=>zeros(Int, sickdaylim, n_conds, agegrps) for i in 1:generations)
        poscontacts = Dict(i=>zeros(Int, sickdaylim, n_conds, agegrps) for i in 1:generations)
        postouched = Dict(i=>zeros(Int, sickdaylim, n_conds, agegrps) for i in 1:generations)

        qloc = (locale=locale, quar_date=thisday)

        # initialize new tracking locales and stash
            tstloc = (locale=locale, test_date=thisday)
            if !haskey(testdat, tstloc)  
                testdat[tstloc] = zeros(Int, sickdaylim, length(conditions), length(agegrps))
            end
            # holds postest people to be quarantined after delay getting test results
            if !haskey(tnt_stash, tstloc) 
                tnt_stash[tstloc] = zeros(Int, sickdaylim, length(test_conds), length(agegrps))
            end
        
        density_factor = env.geodata[env.geodata[:, fips] .== locale, density_fac][1]
        conducted = 0 # per gen
        perday_conducted = 0 # per day

        @inbounds for gen in 1:generations
            if gen == 1
                to_test[gen][:] = avail_to_test
            else
                to_test[gen][:] = postouched[gen-1]
                tc_perday -= conducted
            end

            # test
            postests[gen][:], all_tests = simtests(to_test[gen]; tc_perday=tc_perday, sensitivity=sensitivity, 
                            specificity=specificity, infect_prior=0.05, test_pct=test_pct, infectparams=infectparams)
            conducted = sum(all_tests)
            perday_conducted += conducted
            if sum(all_tests) != 0  # track the test cases
                plus!(all_tests, test_conds, agegrps, sickdays, tstloc, dat=testdat)  
            end
            

            # trace contacts  TODO THIS WILL NO LONGER WORK
            target_cf = repeat(view(infectparams.contactfactors,1,:)',4,1) # use unexposed for all rows
            poscontacts[gen][:] = how_many_contacts!(poscontacts[gen], 
                                    postests[gen],  # equivalent to spreaders in spread
                                    avail_to_test,
                                    target_cf, 
                                    density_factor, infectparams=infectparams)  								
			poscontacts[gen][:] = round.(Int, c_comply .* poscontacts[gen])

            # contacts lead to consquential touches that we count
            target_tf = view(infectparams.touchfactors,map2access[unexposed]:map2access[mild], agegrps)
            postouched[gen][:] = how_many_touched!(postouched[gen], poscontacts[gen], 
                                        avail_to_test, test_conds, 
                                        target_tf, infectparams=infectparams)
            if past_contacts # going back "5" pretend days
                up_multiple = floor(Int, sum(shifter(rand(5),0.4, 1.3)))
                postouched[gen][:] = postouched[gen] .* up_multiple
            end

            # isolate positives and stash breakout
            put_in = round.(Int, q_comply .* postests[gen])
            if test_delay > 0 # delay positives into quarantine. (all results delayed, only pos matter)
                tnt_stash[tstloc] .+= put_in
            else
                if breakout_pct != 0.0  # future breakouts from this quarantine cohort
                    breakout!(breakout_pct, put_in, qloc, qdays) 
                end
                t_n_t_quarantine(put_in, qloc::Quar_Loc; opendat=opendat, isodat=isodat, infectparams=infectparams)
                push!(tntq, (day=day_ctr[:day], quarantine=sum(put_in)))
            end
        end  # for gen 

        # statistics
        avail = reduce(+, map(sum,values(to_test)))
        sumtests = reduce(+, map(sum,values(postests)))
        tc_perday > avail && (@warn "Happy Day $(day_ctr[:day]): more tests available than people to test")
        push!(tntq,(day=thisday, avail=avail, 
                    conducted=perday_conducted, postests=sumtests, 
                    poscontacts=reduce(+, map(sum,values(poscontacts))), 
                    postouched=reduce(+, map(sum,values(postouched)))))

    end  # if start_day
end


function simtests(to_test; tc_perday=1000, sensitivity=.9, specificity=.9, infect_prior=.05, 
    test_pct=.95, infectparams=infectparams)

    # distribute the tests across disease conditions by age group and condition
           # we could do probabilistically but the probs are very small and the whole thing
           # is somewhat artificial: we can capture "randomness" by randomly ignoring x% of the results
    pos_results = zeros(Int, 25,4,5)

    if tc_perday <= 0  # earlier generations of test and trace used up the available tests today
        return zeros(Int,sickdaylim, 4, length(agegrps)), [Int(0)]
    end
    if sum(to_test) == 0
        return zeros(Int,sickdaylim, 4, length(agegrps)), [Int(0)]
    end        

    # today_tests = rand(Binomial(tc_perday, test_pct), 1)[1]

    # println("today $(day_ctr[:day]) tc_perday $tc_perday  today tests $today_tests")

    alloc_pct = reshape(to_test ./ sum(to_test), length(to_test))
    alloc_pct[isnan.(alloc_pct)] .= 0.0  # eliminate NaNs (underflow)
    @assert isapprox(sum(alloc_pct), 1.0) "probabilities must sum to 1.0"

    # dist_tests = round.(Int, (alloc_pct .- 1e-5) .* today_tests) 
    dist_tests = zeros(Int, size(to_test))
    dcat = Categorical(alloc_pct)
    x = rand(dcat, round(Int, test_pct * tc_perday))
    dist_tests[:] = reshape([count(x .== i) for i in 1:length(dcat.p)], size(dist_tests))
    dist_tests[:] = clamp.(dist_tests, Int(0), Int.(to_test))

    # test results  

        # specificity tests apply to actual true not-infected: unexposed recovered
        false_pos = rand.(Binomial.(dist_tests[:, 1:2, :], 1.0 - specificity)) # @views 

        # sensitivity tests apply to actual true infected: nil, mild
        true_pos = rand.(Binomial.(dist_tests[:, 3:4, :], sensitivity))  # @views 
        false_neg = dist_tests[:, 3:4,:] - true_pos   # TODO should report this for curiosity

        pos_results[:] = cat(false_pos,true_pos,dims=2)  # (25,4,5)

        # println(" day $(day_ctr[:day])  ")
        # println(" False Pos  quarantined even though not sick: ", sum(false_pos) ,", ", 
        #         round(sum(false_pos) / sum(dist_tests[:, 1:2, :]), digits=4))
        # println(" True Pos  quarantined because actually sick: ", sum(true_pos) ,", ", 
        #         round(sum(true_pos) / sum(dist_tests[:, 3:4, :]), digits=4))
        # println(" False Neg  not quarantined even though sick: ", sum(false_neg) ,", ",
        #         round(sum(false_neg) / sum(dist_tests[:, 3:4,:]), digits=4))
        # println(" Total positive results, total tests conducted ", sum(pos_results), ", ", 
            # sum(dist_tests))

    return pos_results, dist_tests
end


function bayes(sensitivity, specificity, pr_pop) 
    # Bayesian probs. better but public tests don't use Bayes interpretation!
    pr_pos_given_pos_test = (
                       sensitivity * pr_pop /
       (sensitivity * pr_pop + (1.0 - specificity) * (1.0 - pr_pop))
    ) 
end


function t_n_t_quarantine(postests, qloc::Quar_Loc; opendat, isodat, infectparams)
    if !haskey(isodat, qloc)  
        isodat[qloc] = zeros(Int, sickdaylim, length(conditions), length(agegrps))
    end

    test_conds = [unexposed, recovered, nil, mild] 
    isolate_by!(postests, test_conds, agegrps, sickdays, qloc, opendat, isodat)
end


function t_n_t_unquarantine(cnt, qloc::Quar_Loc, opendat, isodat, infectparams)
    ret_conds = [unexposed, recovered, nil, mild, sick, severe] 
    unisolate_by!(cnt, ret_conds, agegrps, sickdays, qloc, 
                  opendat, isodat, mode=:plus) # delete the qloc when unq all
end


function breakout!(breakout_pct, inq, qloc, qdays; 
    breakout_dist = [0.0, 0.0, 0.0, 0.0, 0.0, 0.03, 0.03,            # TODO fix to make sure 
                    0.03, 0.03, 0.03, 0.03, 0.05, 0.05, 0.05, 0.67]) # (15,) must match qdays

    # holds a breakout array (14,) from the quarantine
    if !haskey(tnt_stash, qloc) 
        tnt_stash[qloc] = zeros(Int, qdays-1)  # no breakout last day
    end

    m = breakout_pct / sum(breakout_dist[1:14])
    breakout_dist[1:14] .*= m
    breakout_dist[15] = 1.0 - sum(breakout_dist[1:14])

    dcat = Categorical(breakout_dist)
    outs = rand(dcat, sum(inq))
    outs = [count(outs .== i) for i in 1:length(breakout_dist)][1:length(breakout_dist)-1]
    tnt_stash[qloc][:] = tnt_stash[qloc] .+ Int(outs)
end  


"""
    function cnt_2_array(cnt, col_mat)

Distribute a count of people to be unquarantined across
the cells of a population matrix, limited to the number of
quarantined people in each cell of the quarantine population
matrix.
"""
function cnt_2_array(cnt, col_mat; ret_conds=[nil, mild, sick, severe, unexposed, recovered])
    (ls, _, as) = axes(col_mat)
    map2unq = (unexposed=1, infectious=-1, recovered=2,  dead=-1, nil=3, mild=4, sick=5, severe=6)
    new_mat = zeros(Int, length(ls), length(ret_conds), length(as))
    for l in ls, c in ret_conds, a in as
        if cnt <= 0; break; end
        put = clamp(cnt, 0 , col_mat[l,c,a])
        new_mat[l,map2unq[c],a] = Int(put)
        cnt -= put
    end
    return new_mat
end


const nxt = [1] # a ref to an integer used as a circular counter

function get_next!(cval) # update the ref value
    cval[] = cval[] < 4 ? cval[] + 1 : 1
end