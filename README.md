# CovidSim

This is a classic SEIR (Susceptible, Exposed, Infected, Removed) simulation of the COVID outbreak of 2019-2020 with some new twists, written in the Julia programming language. [Look at some preliminary results...](https://github.com/lewisl/CovidSim/blob/master/reports/report%201/report%201.ipynb)

There are now both a group model, in the repository Covidgroupsim, and this individual level model here in this repository Covidilmsim.

The group model tracks groups of people in 8 categories by day in a given locale (city or region). Each locale has its own data structure for tracking people.  Multiple locales can be simulated in a single run. The individual model tracks individuals in a locale. For each individual a status and condition are tracked, as well as other traits and outcomes. 

The individual level model ("ilm" from now on) tracks each individual in a locale with an individual's specific traits and outcomes. Because a simulation doesn't know the individual people we use the same age groups, disease status conditions, and sickdays as in the group model. But, the ilm enables more complex policy scenarios to be simulated with more understandable logic. But, the ilm runs slower than the group model because each individual must be queried and updated.

The rest of this readme describes the individual level model.

The statuses are:

- Unexposed
- Infectious (summary of the 4 disease conditions)
- Recovered
- Dead
- and nonsick--to correspond to conditions of degrees of illness.

The conditions are:

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

	The disease spreads from those who are infected to those who are not. Transmissibility varies with the number of days that someone spreading has had it, with asymptomatic transmission assumed. Susceptibility varies by age group of the recipient.

- Transition

	A person who is sick with the virus transitions through the stages from nil to either recovered or dead, based on user-defined decision trees that vary by age group.

Basic tracking includes cumulative data series for each group, new daily values for each group, and detailed daily progression of spreading.  Charts are defined for cumulative data, daily data, and spreading progression.

There are many input parameters that control the behavior of the simulation. Key parameters that affect spreading are:
- contactfactors 

    Determine the number of people that infectious "spreaders" contact, on average, per day. These vary by age group and disease condition of the spreader.

- touchfactors 

    Determine the probability that a contact is consequential--significant enough to *potentially* transmit the virus 

- sendrisk and recvrisk 

    Determine the probability of actually transmitting the virus from sender, which varies by number of days the person has had the disease, and the probability of infecting the recipient, which varies by age group.

- r0 simulation

    The model is more complicated than assuming one r0 applies to the entire population. R0 is *not* an input; it is an outcome.  The factors above provide different effective transmission rates for different age groups, disease conditions, and stage of infection. The r0 simulation provides a sanity check on transmission to see the resulting r0 for a single cohort that includes all age groups and sickdays. The model defaults provide for an early stage R0 of roughly 1.8. (Early stage assumes that the infected group is small relative to the population so that transmission is *not* affected by a large group of non-susceptible people, who may be dead, recovered, or already infected).  The r0 simulation can be run "mid-stream" during a simulation to see how case scenarios and epidemic dynamics change shortrun r0, which is as much socially determined as biologically.

Transition of infected individuals (in the disease cell groups above) is controlled by input decision trees:

- decision trees 

    1 per age group, determine when the condition of an infected person shifts from nil, to mild, to sick, to severe, to recovering or dying. The tree provides different paths through the stages in varying number of days and probabilities, including skipping disease states.

- sanity check on decision trees

   Each decision tree (for an age group) must insure that all infected individuals resolve to recovered or dead at the end of the maximum sickday period (25 days by default). Total probability across all outcomes must sum to 1.0.  The sanity check can be quickly run on a set of decision trees for all 5 age groups. In addition to verifying that probabilities sum to 1, this provides the expected (mean) % of recovered and dead by age group, which can be compared to reported clinical outcomes. [Read more...](https://github.com/lewisl/CovidSim/blob/master/documentation/decision%20tree%20concept.md)

A benefit of the model is comparative ease for running a variety of test cases to examine the response of disease progression to events and potential policy interventions:

- seeding
  
    Travel modeling is planned. In the meantime, seeding events can be defined to introduce infectious people to a locale.  This can occur on any day and introduce people of any condition or sickday.  This enables "manually" causing travel of the disease to new locales when multiple locales are simulated. Multiple seeding events can easily be included in a single simulation run.

- isolation

    With a simple callback function approach, people can be isolated on a given day in a given locale and can be "un-isolated" later.

- social distancing

    The factors that drive spread of the virus can be changed with a complying and non-complying group. A subsequent "event" can change the degree of social distancing and the compliance to simulate varying degrees of "opening up."

- test, trace and isolate

    A group-based SEIR model cannot trace individual testing and outcomes but we can distribute tests for breadth, determine outcomes for the tested group, determine contacts, isolate those with positive test results and repeat through multiple generations of contacts. Many factors can be set such as test capacity per day, test compliance, contact compliance, early "breakout" from isolation, and sickday time to receive test results.

##### Epidemiological Models
There are several different approaches to epidemiological models that have been developed for a long time and various experiences applying models to the current COVID-19 epidemic have been reported. This paper summarizes the various model approaches applied to COVID-19 *non-judgmentally*[1]. Time series forecasting of the most rigorous kind applied correctly to reported infection and death data through as late as May 1, 2020 seems challenged by incompleteness of data as both infections and deaths may be seriously under-reported[2]. SEIR simulations have different challenges because their input parameters, which  represent social behaviors and clinical factors,  are difficult to define given unknowns about the disease. Attempting to correlate the two kinds of models is difficult: time series forecasts are subject to data quality challenges; SEIR models differ substantially from reported data. At this juncture, it may be more important to understand the dynamics of the epidemic that SEIR models provide, while critically examining such models for plausibility.

[1] "Prediction and analysis of Coronavirus Disease 2019," Lin Jia, Kewen Li， Yu Jiang, Xin Guo, and Ting zhao, https://arxiv.org/abs/2003.05447

[2] "Correcting under-reported COVID-19 case numbers: estimating the true scale of the pandemic," Kathleen M. Jagodnik, Forest Ray, Federico M. Giorgi, and Alexander Lachmann, medRxiv pre-print, https://www.medrxiv.org/content/10.1101/2020.03.14.20036178v2

