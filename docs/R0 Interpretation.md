### R0 Interpretation

R0 in the model is *not* a direct input as it is in some differential equation models. Instead, it is the result of the ```spread!``` function, which models interaction of people by group (agegroup, disease condition) and the probability of virus transmission.  So, it falls out from the simulated interactions, much as R0 in the real world can't be observed directly, but is the result of the transmissibility of the virus and the social patterns of people's interactions.

In the function ```how_many_contacts!```, the contactfactors determine how many people are contacted by infected people we refer to as "spreaders." In the function ```how_many_touched!``` touchfactors determine whether a contact is "consequential” and also distributes the contacts across the composition of the population.  Those who might be touched include not only the susceptible, but also people who have recovered and other infected people. They are "out there," too. This reduces the transmission of the virus because the spreader has a lower probability of contacting someone who is susceptible. This is the social dynamic of R0: as the body of people who acquire immunity by recovery or by vaccination (currently not modeled) grows larger, the effective R0 goes down because fewer people get the disease.

The combination of spreaders making contacts and recipients receptive to a touch is based on pseudo-random sampling from 2 distributions:

- The contactfactors input provides the scale input to a gamma distribution sample which produces a distribution of number of contacts. A larger contact factor results in wider dispersion of the sample output--e.g., a longer tail to the right while the mode is close to the left--typically between 0 and 2.

- The touchfactors input provides the "probability of success" input to a pseudo-random sample from the binomial distribution that determines the success of the touch.

The biological transmissibility of the disease, whether by respiratory droplets or aerosols or acquisition from physical surfaces, is still not fully understood. The model multiplies 2 probabilities: 

- the ```sendrisk``` is based on the number of days the spreader has had the disease including days when the infected person is asymptomatic (“nil”) 
- the ```recvrisk``` is based on the age group of the recipient

The inputs above are judgment inputs that are "sanity checked" to produce R0 values in the ranges that epidemiologists have calcuated (or assumed?) based on observations of the early stages of the spread of the Coronavirus in several countries. The tendency is that younger people, who believe that they are not infected, even though in the simulation--and in reality--they may be infectious asymptomatic or with very mild symptoms, will make more contacts. Younger recipients without observable symptoms will, likewise, be more accessible.  Older people and people in the disease categories of sick and severe will make and accept a lower rate of contact.

While R0 is very hard to "observe" in the wild, it is a useful diagnostic to interpret how aggressive a model is in spreading the disease. The model provides an R0 simulator, which calculates the R0 from a stage 1 cohort of spreaders; traces the spreaders through all 25 days of their disease; and simply adds up how many people get the virus from this stage 1 cohort.  Those infected "disappear" out of this r0 mini-simulation so that their spreading is not added to the spreading of the initial cohort. Why a cohort?  Why not just 1 person? We want to have a dispersion of people across age groups and disease conditions because we don't know what the characteristics of a single "patient 0" would be. Thus, we calculate an average across the entire range of characteristics in roughly the distribution that they occur.

This shows some of the difficulties of both calculating and interpreting R0. By moving this cohort through all days of the progression of their disease we have a range of diverse spreaders. Some make many contacts; some make only a few; some recover and stop spreading before 25 days are up; some die. This corresponds to the progress of transmission during the epidemic.

How, then, should we interpret R0?  If we had younger healthier, people with many social interations as the only spreaders, the disease would spread more rapidly with a higher R0.  If we had only older people, who rapidly became sick, they would restrain their contacts and would spread the disease to fewer people. The mix of people who progress through the disease differently for different durations results in an R0 calculation that is effectively an average as occurs in actual epidemic.

The other major factor is the context in which the virus spreads. Social distancing, quarantining, use of masks (of ambivalent benefit), ending mass events could reduce the extent of social contacts and the spread of the infection. This reduces effective R0. As the percentage of the population that has recovered--thereby acquiring at least partial immunity of some duration--increases, the spreaders encounter more people who are less likely to become infected. This is the strongest impact on R0 in the absence of an effective, broadly available vaccine. A large enough increase in recovered people resuls in a decline in the growth rate of active cases of the infection. This is sometimes called herd immunity.

Herd immunity also is amenable to varying interpretations. First, it does not result in immunity. People will still get sick--though progressively fewer at a slower rate. The disease continues for a time. But, the term does mean that the spread slows and eventually becomes a trickle. This provides simple math for imputing what R0 was. When active infected cases begin to decline, the observed percentage of the population that has acquired immunity--if you could observe it--can be used to calculate the average cumulative R0 up to that point:

- call ```ip``` the observed immunity percentage when infected cases begin to decline
- then ```R0 = 1 / (1 - ip)```
- if you actually believe you know R0, then you can calculate ```ip = (R0 - 1) / R0```

So, there is no absolutely correct, static value for R0. It is still useful to have a realistic approach to estimating R0 multiple times. We can judge if we think a simulation is running the infection too hot, not hot enough, or in a range that seems to correspond to observed disease progression. That's important because this simulation has lots of "knobs and dials" that influence the rate of disease transmission. The `r0_sim` function calculates R0, using the single cohort sample described above. We can calculate R0 using various input parameters to pick suitable parameters.  We can also calculate R0 every 10 days of the simulation to see how it changes with case scenarios and the simulation context as it changes over the days.

The table below shows a range of R0 values.

-  columns are scale values for contacts by spreaders
- rows are Binomial probabilities that a contact becomes a “consequential” touch 

This table slightly over-estimates R0 because it applies only at the start of disease progression before there are many recovered people.

|          | 1.1  | 1.2  | 1.3  | 1.4  | 1.5  | 1.6  | 1.7  | 1.8  | 1.9  | 2.0  |
| -------- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| **0.18** | 0.23 | 0.47 | 0.47 | 0.49 | 0.55 | 0.64 | 0.65 | 0.68 | 0.68 | 0.73 |
| **0.23** | 0.28 | 0.53 | 0.61 | 0.62 | 0.65 | 0.69 | 0.73 | 0.79 | 0.82 | 0.83 |
| **0.28** | 0.33 | 0.61 | 0.66 | 0.7  | 0.79 | 0.8  | 0.83 | 0.9  | 0.95 | 0.99 |
| **0.33** | 0.38 | 0.7  | 0.74 | 0.85 | 0.84 | 0.94 | 0.98 | 1.04 | 1.08 | 1.11 |
| **0.38** | 0.43 | 0.8  | 0.85 | 0.89 | 0.93 | 1.03 | 1.11 | 1.16 | 1.2  | 1.27 |
| **0.43** | 0.48 | 0.88 | 0.91 | 0.99 | 1.03 | 1.16 | 1.23 | 1.26 | 1.32 | 1.42 |
| **0.48** | 0.53 | 0.97 | 1.06 | 1.08 | 1.18 | 1.26 | 1.27 | 1.42 | 1.47 | 1.52 |
| **0.53** | 0.58 | 1.01 | 1.09 | 1.17 | 1.25 | 1.33 | 1.43 | 1.52 | 1.52 | 1.68 |
| **0.58** | 0.63 | 1.11 | 1.2  | 1.25 | 1.38 | 1.42 | 1.5  | 1.65 | 1.75 | 1.78 |
| **0.63** | 1.11 | 1.2  | 1.25 | 1.38 | 1.42 | 1.5  | 1.65 | 1.75 | 1.78 | 1.95 |

