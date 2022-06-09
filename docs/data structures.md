# Data Structures
This document provides overview descriptions of the major data structures used by the simulation. For input parameters both the parameter input files and the program data structures are described.

### Population Data
The population data is a table containing a row for each person with columns for the traits of each person. In Julia, the actual datastructure is a TypedTable that allows grouping columns of different types. It is created by function pop_data() in file setup.jl. It contains that traits of each person in a given locale for the current day of the simulation. Only  columns that contain a vector record the history of a person.

##### Trait Columns
|   Column |   Type   |  Description                            |
| ----------  | --------- | ----------------------------------- |
| **status** | Enum{Int} | one of unexposed, infectious[1], recovered, dead |
| **agegrp** | Enum{Int} | one of age0_19, age20_39, age40_59, age60_79, age80_up |
| **cond** | Enum{Int} | one of uninfected, nil, mild, sick, severe |
| **duration** | Integer | no. of days a person has been sick, from 1 to 25 |
| **variant** |  Symbol | provdied by parameter inputs, currently :base, :alpha, :delta, :omicron_ba1, :omicron_ba2 |
| **sickday** | Vector{Int} | array of each day a person became sick, to account for reinfection |
| **recovday** | Vector{Int} | array of days a person has recovered |
| **deadday** | Integer | day of death |
| **ring** | Symbol | currently not used, but to encode groups of people more likely to interact with each other than people in another ring |
| **sdcomply** | Symbol | social distancing case applicable to a person |
| **vaxstatus** | Symbol | one of :none, :first, :multiple, :full, :booster |
| **vaxrcvd** | Vector{Symbol} | vector of vaccines received. Currently one of :JnJ, :Moderna, :Pfizer |
| **vaxday** | Vector{Int} | vector of days a vaccine shot was received |
| **tested** | Bool[2] | true or false if a person has been tested in "test and trace" |
| **testday** | Integer | day of test |
| **quar** | Bool | true or false if a person is in quarantine |
| **quarday** | Integer | day a person starts quarantining |

Here is the output of a fragment of a popdat table

```
Table with 17 columns and 95626 rows
      status     agegrp   cond        duration  variant   sickday  recovday  ⋯
    ┌──────────────────────────────────────────
 1  │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 2  │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 3  │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 4  │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 5  │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 6  │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 7  │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 8  │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 9  │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 10 │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 11 │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 12 │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 13 │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 14 │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 15 │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 16 │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 17 │ unexposed  age0_19  uninfected  0         Symbol[]  [0]      [0]       ⋯
 ⋮  │     ⋮         ⋮         ⋮          ⋮         ⋮         ⋮        ⋮      ⋱
```



### History Series
The history of outcomes during a simulation is recorded in a history series with a row for each day and a column for each statistic being collected. The collection of columns is a TypedTable. The tables are collected in a Dict[3] with a key of locale. For each locale in the dict there are two tables accessed as a NamedTuple (a kind of "on-the-fly" struct) with fields for cum and new. The cum table represents the current value of a statistic on a specific day of the simulation. The new table represents the change in a statistic on a specific day of the simulation. There are lots and lots of columns. This is not that expensive in time or storage because the number of rows equals the number of days of the simulation. Updates are done at the end of each day across the entire population. Group totals for related columns are done only once per locale at the end of the simulation. Each statistic is collected by age group to enable comparisons by age group. Age groups are totaled for each statistic. All vaccine columns by age group and total are summed into totvaccinated. The columns for each infectious condition by age group and total are summed into totinfected, which is also the total for all variants.

The calday column has a date for the day. The simulation internally uses days from 1 to the last day. We allow the simulated first day to be aligned with a specific calendar day to allow the simulation to be compared to reported data. The patterns of the output should be plausible but don't expect the actual values to match closely at all.

##### History Series Column Groups

| Column Group | Values                                              | no.  |
| ------------ | --------------------------------------------------- | ---- |
| Status       | unexposed, infectious, recovered, dead by agegrp and total |  24    |
| Cond         | nil, mild, sick, severe, totinfected by agegrp and total | 30     |
| Vaccine      | Pfizer, Moderna, JnJ, totvaccinated by agegrp and total |   24   |
| Variant      |  base, alpha, delta, omicron_ba1, omicron_ba2 by agegrp and total  |  30    |

Here are examples of fully qualified column names unexposed_age80_up, totvaccinated_age60_79, totvaccinated_total. Unfortunately, there is no hierarchical grouping of columns in a TypedTable. Separate tables could be used but it is easier to programmatically reference the columns for updating as one table.

Here is the path for accessing history tables

```
-dict- locale
series[38015].cum.unexposed_age_0_19   
                 |                  |  
---NamedTuple----|                  |  
---- Typed Table columns------------|     

```

Here is the output of a fragment of a history series table

```
Table with 109 columns and 180 rows
      calday       unexposed_age0_19  unexposed_age20_39  unexposed_age40_59  ⋯
    ┌──────────────────────────────────────────
 1  │ 2020-01-01  0                  0                   0                   ⋯
 2  │ 2020-01-02  0                  0                   0                   ⋯
 3  │ 2020-01-03  0                  0                   0                   ⋯
 4  │ 2020-01-04  0                  0                   0                   ⋯
 5  │ 2020-01-05  0                  0                   0                   ⋯
 6  │ 2020-01-06  0                  0                   0                   ⋯
 7  │ 2020-01-07  0                  0                   0                   ⋯
 8  │ 2020-01-08  0                  0                   0                   ⋯
 9  │ 2020-01-09  0                  0                   0                   ⋯
 10 │ 2020-01-10  0                  0                   0                   ⋯
 11 │ 2020-01-11  0                  0                   0                   ⋯
 12 │ 2020-01-12  0                  0                   0                   ⋯
 13 │ 2020-01-13  0                  0                   0                   ⋯
 14 │ 2020-01-14  0                  0                   0                   ⋯
 15 │ 2020-01-15  0                  0                   0                   ⋯
 16 │ 2020-01-16  0                  0                   0                   ⋯
 17 │ 2020-01-17  0                  0                   0                   ⋯
 ⋮  │     ⋮               ⋮                  ⋮                   ⋮           ⋱
```

## Input parameters and data

### Transition Trees

Accessing a  transition array from a transition tree uses this path in Julia

```
 
--dict------   variant
transitionset[omicron_ba1].tree.age0_19[14]    
                           |    |
---Transitionparams--------|    |      |  |
------field tree to Agetree-----|      |  |
---------field age0_19 to a dict-------|  |
----------------- dict key duration ---|--| -> pointing
                         to a value that is 4 x 6 array

```

### Notes
[1] Infectious status represents currently exposed and not necessarily infectious. Actually being infectious is based on a probability for each day of the duration of a person's illness

[2] Bool is the type boolean in the Julia language