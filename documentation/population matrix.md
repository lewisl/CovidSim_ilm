# Population Matrix

A population matrix is the data structure used to track groups of people in the simulation. All population matrices have the same indices and sizes. These names are used to indicate the purposes of the matrices:

- opendatmx: holds people who are circulating within a locale
- isolatedmx: holds people who are in quarantine in a locale

History matrices have a slightly different structure:

- cumhistmx: holds the cumulative history by day of the status of the people in the simulation
- newhistmx: holds the new values at each day of the people in the simulation

opendatmx and isolatedmx are dicts:

- keys are locales;
- values are the actual population matrix with structure below

**Population Matrix Structure**

Dimensions are:

- rows: duration, which indicated how many days members of a cohort have had the virus
- columns: conditions => unexposed, infectious, recovered, dead, nil, mild, sick, and severe. nil, mild, sick, and severe are the only columns that use all duration rows as they are the infectious conditions. Infectious is a convenience total of all 4 disease conditions and is never updated directly.
- planes (third dimension value): agegrps.  Currently, we use 5 agegroups:
    - 0-20
    - 20-40
    - 40-60
    - 60-80
    - 80+

You may ask why not use DataFrames or one of the Julia packages that provide named arrays. DataFrames are slightly awkward, a bit slower, and must be 2 dimensional. All of the data in a population matrix is integer so there is no need for mixed types, which is the greatest benefit of DataFrames. The various packages implementing "named arrays" are at different degrees of maturity; add even more complexity than this rough-and-ready solution; and have slightly different APIs. The approach here is in no way general:  it only serves this one application. The package maintainer is probably the only person writing code to manipulate these data structures. One day I'll switch to one of the named array packages.

**Accessing Population Matrices**

While the structure is simple enough and Julia provides convenient syntax for accessing and modifying multi-dimensional arrays, the population matrices should not be accessed directly in case the structure changes (indeed, I changed from locales as a 4th dimension, to locales as a key to a dictionary of 3d arrays--and it all worked!).

Convenience functions enable accessing values, updating in place, adding in place, subtracting in place, and totalling across dimensions. The functions are very simple--you could do all of what they do in short, obvious code--but they isolate your code from the underlying physical layout.

Here are the functions:

- grab: ```grab(condition, agegrp, duration, locale, openmx)```
- input!: ```input!(val, condition, agegrp, duration, locale, dat)```
- plus!: ```plus!(val, condition, agegrp, duration, locale, dat)```
- minus!: ```minus!(val, condition, agegrp, duration, locale, dat)```

**Alternative Mappings**

Some data structures refer to the population conditions either in a different order or leaving out some conditions. To avoid explicit reference by integers, we provide maps that convert the name labels (unexposed through severe) to the appropriate values.

Here is an example of an alternative mapping:

```julia
"""
Map a condition index from rows in the data matrix to
indices for the transition probabilities:


               unexposed  infectious  recovered  dead   nil  mild  sick  severe
data rows         1          2            3        4     5     6     7     8
transition pr    -1         -1            1        6     2     3     4     5


Transition probability indices that return -1 are not used and will raise an error.

- Use with text literal in code as map2pr.nil => 2
- Use with variables that stand for the data rows as map2pr[nil] => 2 # if nil == 5
"""
const map2pr = (unexposed=-1, infectious=-1, recovered=1, dead=6, nil=2, mild=3, sick=4, severe=5)
```