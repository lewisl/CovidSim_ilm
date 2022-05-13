# CovidSim

This is a classic SEIR (Susceptible, Exposed, Infected, Removed) simulation of the COVID outbreak of 2019-2022 with some new twists, written in the Julia programming language. [Look at some preliminary results...](https://github.com/lewisl/CovidSim/blob/master/reports/report%201/report%201.ipynb)

Beginning as group or compartment model, a shift was made to an agent-based model in this repository Covidilmsim. The compartment model contained 1000 groups before considering variants, vaccination, and social distancing. The compartment model ran very quickly but as more conditions and public health policies were added, the logic became very tangled. 

The agent-based model (or individual level model, "ilm", from now on) tracks each individual in a locale with an individual's specific traits and outcomes. Because a simulation doesn't know actual individual people, we use the same age groups, disease status conditions, and infection durations as in the group model. But, the ilm enables more complex policy scenarios to be simulated with more understandable logic. However, the ilm runs slower than the group model because each individual must be queried and updated. The rest of the readme describes the individual level model.

The statuses are:

- Unexposed
- Infectious (summary of the 4 disease conditions)
- Recovered
- Dead

The disease conditions are:

- Nil (infected and asymptomatic)
- Mild
- Sick
- Severe

The agegroups are:

- age0_19
- age20_39 
- age40_59 
- age60_79 
- age80_up

The basic processes of the simulation are:
- Spread

	The disease spreads from those who are infected to those who are not. Transmissibility varies with the number of days that an infected person has had it, with asymptomatic transmission assumed. Susceptibility varies by the age group of the recipient. Variants, vaccination, and recovery from prior infection all affect spread. Partial, temporary immunity conferred by recovering from an infection or being vaccinated declines over time.

- Transition

	A person who is sick with the virus transitions through stages from nil to either recovered or dead, based on user-defined transition arrays that vary by age group.

Basic tracking includes cumulative data series for each group, new daily values for each group, and detailed daily progression of spreading.  Charts are defined for cumulative data, daily data, and spreading progression.

There are many input parameters that control the behavior of the simulation. Key parameters that affect spreading are:
- contactfactors 

    Determine the number of people that infectious "spreaders" contact, on average, per day. These vary by age group and disease condition of the spreader.

- touchfactors 

    Determine the probability that a contact is consequential--significant enough to *potentially* transmit the virus 

- sendrisk and recvrisk 

    Determine the probability of actually transmitting the virus from sender, which varies by number of days the person has had the disease, and the probability of infecting the recipient, which varies by age group.

- r0 simulation

    The model is more complicated than assuming one r0 applies to the entire population. R0 is *not* an input; it is an outcome.  The factors above provide different effective transmission rates for different age groups, disease conditions, and stage of infection. The r0 simulation provides a sanity check on transmission to see the resulting r0 for a single cohort that includes all age groups and durations. The model defaults provide for an early stage R0 of roughly 1.8. (Early stage assumes that the infected group is small relative to the population so that transmission is *not* affected by a large group of non-susceptible people, who may be dead, recovered, or already infected).  The r0 simulation can be run "mid-stream" during a simulation to see how case scenarios and epidemic dynamics change shortrun r0, which is as much socially determined as biologically.

Transition of infected individuals (in the disease cell groups above) is controlled by input transition arrays:

- transition tree of arrays 

    1 per age group, determine when the condition of an infected person shifts from nil, to mild, to sick, to severe, to recovering or dying. The tree provides different paths through the stages of the illness to eventual recovering or dying.

- checkpoints

    Rather than assign probabilities at disease inception to determine the severity and duration of an infected individual's sickness, each infected person progresses through up to 25 days of infection. At set intervals--checkpoints--each infected person probabilistically transitions to a new infection condition, recovery, or death.

- sanity check on checkpoints and transition matrices

   Each decision tree (for an age group) must resolve all infected individuals to recovered or dead at the end of the maximum duration period (25 days by default). Total probability across all outcomes must sum to 1.0.  The sanity check can be quickly run on a set of decision trees for all 5 age groups. In addition to verifying that probabilities sum to 1, this reports the expected % of recovered and dead by age group, which can be compared to reported clinical outcomes. [Read more...](https://github.com/lewisl/CovidSim/blob/master/documentation/decision%20tree%20concept.md)

A benefit of the model is comparative ease for running a variety of test cases to examine the response of disease progression to events and potential policy interventions:

- seeding
  
    Seeding events can be defined to introduce infectious people to a locale.  This can occur on any day and introduce people of any condition or duration.  This enables "manually" causing travel of the disease to new locales when multiple locales are simulated. Multiple seeding events can easily be included in a single simulation run. Travel modeling was initially considered but as Covid progressed we saw that its high transmissibility resulted in community spread nearly everywhere, though in differing time waves. Other groups have done elaborate modeling of travel based on air flight schedules but events on the ground seemingly overwhelmed travel. A disease with different characteristics or at an early or late stage could certainly warrant more serious consideration of travel patterns. 

- isolation

    With a simple callback function approach, people can be isolated on a given day in a given locale and can be "un-isolated" later.

- social distancing

    The factors that drive spread of the virus can be changed with a complying and non-complying group. A subsequent "event" can change the degree of social distancing and the compliance to simulate varying degrees of "opening up" or the impact of mask usage.

- test, trace and isolate

    A group-based SEIR model cannot trace individual testing and outcomes but we can distribute tests for breadth, determine outcomes for the tested group, determine contacts, isolate those with positive test results and repeat through multiple generations of contacts. An individual level model enables easier treatment of test and trace. Many factors can be set such as test capacity per day, test compliance, contact compliance, early "breakout" from isolation, and duration time to receive test results.

##### Epidemiological Models
There are several different approaches to epidemiological models that have been developed for a long time and various experiences applying models to the current COVID-19 epidemic have been reported. This paper summarizes the various model approaches applied to COVID-19 *non-judgmentally*[1]. Time series forecasting of the most rigorous kind applied correctly to reported infection and death data through as late as May 1, 2020 seems challenged by incompleteness of data as both infections and deaths have been seriously under-reported[2]. SEIR simulations have different challenges because their input parameters, which  represent social behaviors and clinical factors,  are difficult to define given unknowns about the disease. Attempting to correlate the two kinds of models is difficult: time series forecasts are subject to data quality challenges; SEIR models differ substantially from reported data. At this juncture, it is more important to understand the dynamics of the epidemic that SEIR models can provide especially modeling various public health interventions, while critically examining the models for plausibility.

[1] "Prediction and analysis of Coronavirus Disease 2019," Lin Jia, Kewen Li, Yu Jiang, Xin Guo and Ting Zhao, https://arxiv.org/abs/2003.05447

[2] "Correcting under-reported COVID-19 case numbers: estimating the true scale of the pandemic," Kathleen M. Jagodnik, Forest Ray, Federico M. Giorgi, and Alexander Lachmann, medRxiv pre-print, https://www.medrxiv.org/content/10.1101/2020.03.14.20036178v2

