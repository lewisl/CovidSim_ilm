#### Immunity effectiveness from vaccination:
    ##### rises for a short time after each shot
        - effect of multiple shots: based on clinical trials
        - number of shots recorded in vaxstatus as one of {:none, :first, :full, :booster}
        - data provides max effectiveness from trial based on vaxstatus and inbound variant
    ##### decays with time per vaccine
        - varies with each kind of effectiveness?
        - varies per inbound variant?
    ##### effectiveness per inbound variant
        - against contagion  (in spread)
        - against severe sickness  (in transition)
        - against death  (in transition)
            

#### Immunity from recovery
    ##### rises for a short time immediately after recovery
    ##### decays with time per acquired variant
       - varies per inbound variant (or capture that with effectiveness variation by variant)?
    ##### effectiveness per inbound variant
        - against contagion (in spread)
        - against severe sickness (in transition)
        - against death (in transition)
            

#### Combination of immunity from recovery and immunity from vaccination

#### Approaches & Questions
1. Combine multiplicatively and then squash between 0.05 and 0.95
    1) combining immunity factors by multiplying makes them get smaller, which doesn't make sense
    2) decay can sensibly be multiplicative to reduce effectiveness
    3) effectiveness per inbound variant can be multiplicative to reduce effectiveness for more contagious or lethal variant
2. factors have ranges between limits and combine with minimum/maximum
3. How much weight to agegroup for immunity effectiveness? for decay?
4. assume decay is the same reduction factor for each kind of effectiveness
    1) the more rapid and severe drop in contagion effectiveness to be captured by range of contagion effectiveness?