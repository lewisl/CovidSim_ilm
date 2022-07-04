
vaccine_template = """
# Include a complete vaccine section for each vaccine
Pfizer:                    # vaccine name will convert to symbol at load time
  halflife:         360
  reqdshots:          2
  delay2ndshot:      21
  full_effect_days:  14
  day1_effect:      0.65


  # reduce vaccine immunity effectiveness for infection risk
  #     compared to effectiveness for reducing severity
  # for variant being received by vaccinated person
  # multiply times effectiveness
  infectfactor:    
    base:     0.9
    alpha:    0.9
    delta:    0.85
    omicron_ba1:  0.75
    omicron_ba2: 0.75

  # effectiveness in reducing severity of illness 
  # shots: first, full, booster
  # we need data for each known variant
  # variants are the attacker variant
  effectiveness:
    first:
      base:     0.90
      alpha:    0.90
      delta:    0.70
      omicron_ba1:  0.65
      omicron_ba2:  0.65
    full:
      base:     0.94
      alpha:    0.94
      delta:    0.80
      omicron_ba1:  0.85
      omicron_ba2: 0.85
    booster:
      base:     0.94
      alpha:    0.94
      delta:    0.80
      omicron_ba1:  0.85
      omicron_ba2: 0.85


"""

variant_template = """
"""

vaxsched_template = """
vaxesincluded:
  Pfizer:
    mix:            0.45    # 0.55
    starting_doses: 70000    # 70000
    pct2ndshot:     0.9
    alternate:      [Moderna, JnJ]  # one of none, any, other shot
    booster:        false
  Moderna:
    mix:            0.45     # 0.30
    starting_doses: 60000      # 60000
    pct2ndshot:     0.9
    alternate:      [Pfizer, JnJ]
    booster:        false
  JnJ:
    mix:            0.1     # 0.15
    starting_doses: 20000      # 20000
    pct2ndshot:     0.0
    alternate:      [Pfizer, Moderna]
    booster:        false

dayrange:       [350, 700]
targetpct:      0.85    # 0.65
filterfunc:
shotmode:       all              # values in :first, :second, :all, :booster
pattern:        [0.0, .02, .05, .10, .15, .19, .21, .16, .08, .03, .01, 0.0]
spreadfunc:
"""

social_params_template = """
#
# gammashape: shape parameter of gamma distribution used to calculate contacts
#     - Float64, typically in [0.6, 1.5]
gammashape: 1.0
#
# contactfactors: gamma distribution size parameter 
#     - determines number of contacts per day for each spreader
#     - converted to Matrix{Float64} 4 x 5 (conditions converted to rows, agegrps to columns)
#
contactfactors:                                        
  age0_19:                                            
    {nil: 1.1, mild: 1.1, sick: 0.7, severe: 0.5}       
  age20_39:
    {nil: 2.1, mild: 2.0, sick: 1.0, severe: 0.6}
  age40_59:
    {nil: 2.1, mild: 2.0, sick: 1.0, severe: 0.6}
  age60_79:
    {nil: 1.7, mild: 1.6, sick: 0.7, severe: 0.5}
  age80_up:
    {nil: 1.0, mild: 0.9, sick: 0.6, severe: 0.5}
#
# touchfactors: probability input to binomial distribution 
#     - determines if target is "touched"
#     - must be between 0.0 and 1.0
#     - converted to Matrix{Float64} 5x6 when loaded
#
touchfactors:
  age0_19:                                                    
    {unexposed: 0.55, recovered: 0.55, nil: 0.55, mild: 0.55, sick: 0.28, severe: 0.18}
  age20_39:
    {unexposed: 0.63, recovered: 0.63, nil: 0.63, mild: 0.62, sick: 0.35, severe: 0.18}
  age40_59:
    {unexposed: 0.61, recovered: 0.61, nil: 0.61, mild: 0.58, sick: 0.30, severe: 0.18}
  age60_79:
    {unexposed: 0.41, recovered: 0.41, nil: 0.41, mild: 0.41, sick: 0.18, severe: 0.18}
  age80_up:
    {unexposed: 0.35, recovered: 0.35, nil: 0.35, mild: 0.28, sick: 0.18, severe: 0.18}
"""
