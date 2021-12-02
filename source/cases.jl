####################################################################
# cases.jl
#       pass these cases to run_a_sim in kwarg runcases as a list
#       cases = [CovidSim.<funcname>, CovidSim.<funcname>, ...]  
#       then run_a_sim(geofilename, n_days, locales; runcases=cases)
####################################################################


####################################################################
# seeding cases
####################################################################

"""
Generate seeding cases.
inputs: day, cnt, sickday, cond, agegrp
Two of the inputs may refer to multiple items and must match in number of items.

Returns a function that can be used in runcases input to run_a_sim.
"""
function seed_case_gen(day, cnt, sickday, cond, variant, agegrp) # these args go into the returned seed! case
    # this gets returned; assign it a value at the cmdline; use as an input to run_a_sim
    function runcase(locale, dat, socialparams, spreadset, sdcases, ages; startofday)  # args must match runcases loop in run_a_sim
        seed!(day, cnt, sickday, cond, variant, agegrp, locale, dat; startofday=startofday)  # payload: this is what the function will do when run
    end
end


"""
    seed!(day, cnt, sickday, conds, agegrps, locale, dat)

This is the action function that setups and implements a seeding case all in one execution.
"""
function seed!(day, cnt, sickday, conds, variants, agegrps, locale, dat; startofday)

    startofday || return  # if true->don't return.  if false->then do return

    @assert length(sickday) == 1 "input only one sickday value"
    # @warn "Seeding is for testing and may result in case counts out of balance"
    if day == day_ctr[:day]
        println("*** seed day $(day_ctr[:day]): $(sum(cnt)) $conds to $locale")
        # @assert (cond in [nil, mild, sick, severe]) "Seed cases must have conditions of nil, mild, sick, or severe" 
        make_sick!(dat[locale]; cnt=cnt, fromage=agegrps, tocond=conds, tovariant=variants, tosickday=sickday)
    end
end


# some generated seed! cases-->these are global (in code)
# seed_6_12 = seed_case_gen(8, [0,6,6,0,0], 5, nil, :base, agegrps)
# seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 5, nil, :base, agegrps)


####################################################################
# isolation cases
####################################################################



