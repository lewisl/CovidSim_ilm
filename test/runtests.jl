using Test
using CSV
using Dates
using Random
using LazyTables
using CovidSim_ilm

function indoor_seq_for_window(day1, ndays, indoor_start, indoor_end; indoor_lift=1.1)
    loc = 1
    caldays = range(day1, step=Day(1), length=ndays)
    geodata = LazyTable(fips=[loc], indoor_st=[indoor_start], indoor_end=[indoor_end])
    series = Dict(loc => (cum=(caldays=caldays,),))

    seq = CovidSim_ilm.build_indoor_seq([loc], ndays, geodata, series, indoor_lift)[loc]
    return seq, collect(caldays)
end

@testset "build_indoor_seq" begin
    @testset "cross-year window follows calendar dates through leap day" begin
        seq, caldays = indoor_seq_for_window(
            Date(2020, 1, 1),
            160,
            Date(1, 9, 15),
            Date(2, 5, 30),
        )

        @test seq[findfirst(==(Date(2020, 2, 29)), caldays)] == 1.1
        @test seq[findfirst(==(Date(2020, 5, 30)), caldays)] == 1.1
        @test seq[findfirst(==(Date(2020, 5, 31)), caldays)] == 1.0
    end

    @testset "same-year window uses inclusive month-day bounds" begin
        seq, caldays = indoor_seq_for_window(
            Date(2020, 3, 1),
            110,
            Date(1, 3, 15),
            Date(1, 6, 15),
        )

        @test seq[findfirst(==(Date(2020, 3, 14)), caldays)] == 1.0
        @test seq[findfirst(==(Date(2020, 3, 15)), caldays)] == 1.1
        @test seq[findfirst(==(Date(2020, 6, 15)), caldays)] == 1.1
        @test seq[findfirst(==(Date(2020, 6, 16)), caldays)] == 1.0
    end
end

@testset "infectrisk clamp" begin
    infectset = Dict(
        1 => (sendrisk=[1.0], recvrisk=[1.0],),
    )

    @test CovidSim_ilm.infectrisk(infectset, 1, 1, :age0_19, 1.0, 1.0) == 1.0
    @test CovidSim_ilm.infectrisk(infectset, 1, 1, :age0_19, 0.5, 0.5) == 0.5
end

@testset "how_many_contacts matches C++ cap" begin
    Random.seed!(12345)
    contactfactors = fill(5.0, 4, 5)
    draws = [CovidSim_ilm.how_many_contacts(1.0, 1.0, 3.0, :age20_39, :nil, contactfactors)
             for _ in 1:200]

    @test all(0 .<= draws .<= 12)
    @test any(==(12), draws)
end

@testset "categorical_sim batch draws" begin
    Random.seed!(12345)
    draws = CovidSim_ilm.categorical_sim([0.0, 1.0, 0.0], 5)

    @test draws == fill(2, 5)
end

@testset "redistribute_probability! matches C++ logic" begin
    probvec = [0.1, 0.2, 0.1, 0.3, 0.2, 0.1]
    CovidSim_ilm.redistribute_probability!(probvec, 0.5, 10)
    @test probvec ≈ [0.2, 0.3, 0.2, 0.15, 0.1, 0.05]

    probvec = [0.1, 0.2, 0.1, 0.3, 0.2, 0.1]
    CovidSim_ilm.redistribute_probability!(probvec, 0.5, CovidSim_ilm.DURATIONLIM)
    @test probvec ≈ [0.4, 0.2, 0.1, 0.15, 0.1, 0.05]
end

@testset "runsim_seed_sweep_to_csv appends seed rows" begin
    model = buildsim(
        20, 38015;
        day1=Date("2020-01-01", "yyyy-mm-dd"),
        dovax=false,
        paramdir="sample_parameters",
        geofilename="geo2data.csv",
        socialfilename="socialparams.yml",
        vaccinefilename="vaccines.yml",
        scheddir="vaccine_100k",
        variantfilename="variants.yml",
    )

    seed20_39_day1 = maketraitseedfunc(; cond=:nil, variant=:base, duration=1,
        filter=[Term(:agegrp, :age20_39), Term(:status, :unexposed)],
        cnt=3, forlocale=0, triggerdate=1, forstartofday=true)
    seed40_59_day1 = maketraitseedfunc(; cond=:nil, variant=:base, duration=1,
        filter=[Term(:agegrp, :age40_59), Term(:status, :unexposed)],
        cnt=3, forlocale=0, triggerdate=1, forstartofday=true)

    mktempdir() do tmp
        rows = CovidSim_ilm.runsim_seed_sweep_to_csv(model, [111, 112];
            pathstr=tmp,
            basedir=:none,
            overwrite=true,
            locale=38015,
            runcases=[seed20_39_day1, seed40_59_day1],
            dovax=false,
            silent=true)

        @test length(rows) == 2

        filepath = joinpath(tmp, "seed_sweep.csv")
        @test isfile(filepath)

        result = CSV.File(filepath) |> collect
        @test length(result) == 2
        @test [row.seed for row in result] == [111, 112]
    end
end
