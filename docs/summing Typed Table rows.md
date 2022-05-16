This is not an issue. I just thought it was a cool approach that took me a while to figure out. 

It's a way to add a selection of columns from a wide-ish table of 108 columns of about 800 rows. The rows are days of a simulation and the columns are outcomes from each day. Most of the outcomes are collected for 5 age group cohorts. The age groups need to be totaled for each outcome each day of the simulation. We can wait until the simulation is done to do the totaling. There are 2 sets of statistics: cumulative and new.

Here is an expression that does it with no allocations! This is one line from a slightly longer function so I'll explain the variables.

```julia
getproperty(newhist, Symbol(item, "_", "total"))[:] .= .+(columns(getproperties(newhist, seriesbyage[item]))...)      
```

```newhist``` is a Table that contains the columns that track daily changes of the simulation. New results are updated into the Table each 'day' of the simulation.
```item``` is the base of column names that look like ```:unexposed_age20_39```.

- The left hand side uses getproperty to reference the total column for one of the outcomes being tracked. 
- ```[:]``` refers to all of the rows and enables the vector to be updated in place with no allocation of a new result--but we also need to make sure that no temporary arrays are allocated by the calcuation on the right hand side.  
- ```.=``` is Julia's notation for fusing the calculation on the right hand side into a loop.
- ```.+``` does element-wise addition of the columns.
- ```columns``` function creates a tuple of the columns of a table. Note that we splat the tuple so that we get the arguments to ```.+```.
- ```getproperties``` is a Table function to create a new table comprised of multiple columns from a wider table. The column names are a tuple in a previously created Dict of the age group columns for each outcome. That tuple looks like this:
```(:unexposed_age0_19, :unexposed_age20_39, :unexposed_age40_59, :unexposed_age60_79, :unexposed_age80_up)``` It is a bit unwieldy and surprisingly expensive in time to create on the fly within the loop. Since it never changes, it is created ahead of time as a Dict because there are 18 such tuples.

So, what this one-liner does is sum up all the rows of the 5 columns with no allocations. Since I need to do this for 18 outcomes, I can do it programmatically in a loop. This is a little bit cumbersome because the syntax for Tables favors literals--```newhist.unexposed_age20_39```--more than programmatic referencing.

Here is the complete function:

```julia
@inline function update_total_agegrps!(newhist, cumhist)
        
    for item in seriesgroups
        getproperty(newhist, Symbol(item, "_", "total"))[:] .= .+(columns(getproperties(newhist, seriesbyage[item]))...)      
        getproperty(cumhist, Symbol(item, "_", "total"))[:] .= .+(columns(getproperties(cumhist, seriesbyage[item]))...)     
    end
    
end
```