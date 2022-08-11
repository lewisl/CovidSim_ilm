#### Decision Tree Concept

In this SIR simulation of the COVID-19 outbreak (compartments for Susceptible, Infected, Removed (Recovered or Dead) further divided by age group), I originally used a decision tree to progression people who have gotten infected from nil (asymptomatic), to mild, to sick, to severe, to recovered, or (sadly) to dead. The simulated folks don't go through all conditions.  After a number of days of being infectious they reach progression stages every few days that probabilistically send them to a different condition or terminal state over the course of 25 days. The terminal states are outcomes of either dead or recovered.  It's essential that all the probabilities for all sequences leading to each terminal state must sum to 1.0.  Each age group has its own progression process so I sanity check all of them to make sure the probabilities of recovering or dying total 1.0.

Originally, I implemented this concept as actual trees in an acyclic top-down graph. However, in the simulation logic the nodes had no real meaning and weren't used. Instead, the model checked each infected person on specific days to determine if the person's condition should change based on the probabilities of each outcome, given the day of the progression and the person's current condition. Eventually, I changed from a tree approach to progression stages, current conditions, and probability of progression to a new condition (still sick) or terminal status (recovered or dead). 

The branches and leaves were gone. However, the probabilities of all progression sequences from nil (or asymptomatic illness) to a terminal state must sum to 1.0 and each person must reach a terminal state within the allotted number of simulation days. I changed the sanity check approach from navigating branches to leaf nodes to finding all sequences leading to terminal states. With some changes to the implementation, this is conceptually the same problem.

I thought I *should* be able to do this recursively, but my brain is more loopy than recursive.  I had done the text book recursion cases of factorial and following the branches of a binary tree with the classic paths down left and right branches. But, the progression sequences are more ragged.  There are a different number of outcomes at each progression point (though always 1 to 3) and each sequence can have a different number of steps (longest is 7).  For 2 days I beat myself up, got impatient, and gave up.  I built the paths by hand (the human brain can see the paths almost immediately) and then fed the paths as input to the sanity check. I put the problem down to get more important stuff done.

Later, I took it on again.  In half an hour in < 20 lines of code I solved it. It's not recursive, but it is conceptually recursive (I did learn something...) or we might call it pseudo-recursive.  The trickiest thing wasn't the recursive-lite approach. I realized I need to append to the current sequence without modifying it. In ML-like languages we'd do this with "cons", which creates a pair.  In ML, one abhors mutating a variable so you can only call with a new value (the argument) or return a new value. Recursion performs the loop with input parameters to the next recursive call being the changed value. Alternatively in Julia, we love to mutate data structures in place for efficiency.  There is append! to extend an array in place, but there is no append to create a new array. But, it's easy to use vcat create a new array that is an append--or "cons" of the old array (the head) and its new tail. 

Recursion requires an explicit test with a branch: one branch calls the recursive function again and the other stops recursion by (usually) returning a constant value to stop adding to the stack and start returning values for all of the pending calls on the stack. A while loop is a bit like recursion:  we stop when the conditions have been met that we are done, not by a fixed number of iterations. This while loop runs 17 times for 31 distinct sequences. Each time through we finish at least one sequence of a given length. 

```julia
while !isempty(todo)
   seq = popfirst!(todo)  
   lastnode = seq[end]
   breakday, fromcond, tocond = lastnode
   nxtidx = findfirst(isequal(breakday), breakdays) + 1
   for brk in breakdays[nxtidx:end]
      if tocond in keys(dt_by_age[brk])   # keys are the fromcond at the next break day so previous tocond == current fromcond
         for i in eachindex(dt_by_age[brk][tocond][:outcomes])
            outcome = dt_by_age[brk][tocond][:outcomes][i]
            prob = dt_by_age[brk][tocond][:probs][i]
            newseq = vcat(seq, (duration=brk, fromcond=tocond, tocond=outcome, prob=prob))
            if (outcome == dead) | (outcome == recovered)  # terminal node reached--no more nodes to add
               push!(done, newseq)
            else  # not at a terminal outcome: still more nodes to add
               push!(todo, newseq)
            end
         end
         break # we found the tocond as a matching fromcond
      end
   end
end
```



This is conceptually like recursion. Each addition to a path is tested: Instead of pushing things onto the stack, I push them onto the todo list. I pop an item from the todo list and there are 2 outcomes: either the next outcome at a stage is a terminal condition and that sequence is done; or the possible outcome will lead to another and the sequence is one more stage longer and gets pushed onto the todo list.

It's sort of inverted recursion because I extend the sequence by a stage before deciding what to do with it. In recursion, the completion of the paths comes up as the recursive function returns static results instead of making another call. But, it's short, simple and relatively obvious because at each step there are only 2 choices.  I am sure it's not very efficient and there are lots of allocations, but recursion wouldn't be that efficient and there would be lots of allocations pushed onto the stack. But, the progression steps are few--see the one below.

It's a bit tricky because a tree is a bit of a gnarly nested data structure: 

- a tree is a dict of progression stages--a day on which a new outcome must be chosen;
- each stage is a dict of current conditions for an infected person;
- for each current condition is a dict with a discrete distribution of next possible outcomes and their probabilities (it's a distribution because the probabilities sum to 1.0)
- terminal outcomes either recovered or dead;

Sequences proved to be a bit more complicated than a tree of nodes. In the tree approach, each node contained a pointer to the next node so it was easy to walk the tree. With the sequence approach we extend to multiple new sequences by adding each outcome to the current sequence. This is messier than walking paths of a tree. But, if you look at the innermost if-else, you'll see that we either have reached the end of a sequence, and add it to the "done" list; or we have found another intermediate sequence and we add it to the "memo" list, todo. 

Another complication is that decision stages to build a sequence might not be strictly sequential. A person with a given condition might skip a stage until "her time comes up" and she faces probabilistic distribution to her next condition. This means we might have to search more than one subsequent progression stage. This memo-izes nicely because we need only search subsequent progression stages--a declining number of searches.

Now, I can quickly tell for all the trees if the branch endpoints all add up to 1 and get the expected value proportions for recovering or dying.  This was really satisfying, if slightly morbid.

Here is what a tree looks like for a single age group:

```julia
#=
age0_19:                                # one age group
  5:                                    # progression day--number of days a person has been sick
    nil:                                # the person's current condition
      probs: [0.4, 0.5, 0.1]            # probabilities of each potential outcome
      outcomes: [nil, mild, sick]       # outcomes can be an illness condition or a terminal status
  9:                                    # another progression day...
    nil:                                 
      probs: [0.9, 0.1]
      outcomes: [recovered, sick]
    mild:
      probs: [1.0]
      outcomes: [mild]
    sick:
      probs: [0.95, 0.05]
      outcomes: [sick, severe]
  14:
    mild:
      probs: [1.0]
      outcomes: [recovered]
    sick:
      probs: [0.85, 0.12, 0.03]
      outcomes: [recovered, sick, severe]
    severe:
      probs: [0.692, 0.302, 0.006]
      outcomes: [recovered, severe, dead]
  19:
    severe:
      probs: [0.891, 0.106, 0.003]
      outcomes: [recovered, severe, dead]
  25:
    sick:
      probs: [0.976, 0.024]
      outcomes: [recovered, dead]
    severe:
      probs: [0.91, 0.09]
      outcomes: [recovered, dead]
=#

```


Here is the output for a valid collection of 5 trees, one for each age group. The leaf node probabilities all sum to one (3rd column is total probs for "recover" and 4th column is total probs for "dead").

```
5×4 Array{Float64,2}:
 1.0  1.0  0.999511  0.000488522
 2.0  1.0  0.999332  0.000668336
 3.0  1.0  0.99616   0.0038404
 4.0  1.0  0.980863  0.0191366
 5.0  1.0  0.850242  0.149758
```