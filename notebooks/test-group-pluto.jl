# ## A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

begin
	# using CovidSim_ilm
	push!(LOAD_PATH, joinpath(homedir(), "Dropbox/Online Coursework/Covid/group-src")) 

	try
		using Revise
	catch e 
		@warn(e.msg)
	end
end

using CovidSim_group

seed_1_6 = seed_case_gen(1, [0,3,3,0,0], 1, nil, :default, agegrps)

result_dict, series = run_a_sim(180, 38015, showr0=false, silent=true, runcases=[seed_1_6]);

cumplot(series, 38015)

keys(result_dict)

result_dict["dat"]["openmx"][38015]
