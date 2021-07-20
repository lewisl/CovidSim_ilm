##### For spreading as simplified for contact tracing



        # how many people are contacted by each spreader?  Think of this as reaching out...
        # contacts is the potential number of people contacted by a spreader in each
        # cell by sickday (sickdaylim), infectious cond (4), and agegrp(5)

    #=  This originally ignores the conditions of the touched--assumes they are 
        all equally likely to be touched
        how_many_touched corrects this.
        We assume spreaders is small compared to all_accessible. At some point this might not be true:
        how_many_touched also handles this.
    =#


    # now to business: who gets touched in all conditions and agegroups?
    #=
        - folks in each cell of poscontacts touch a sample of the accessible reduced by the touch factor
        - draw a categorical sample for each cell to distribute them across the contact categories
        - we should care about not touching more than the number of accessible
        - touchfactors is (5,6): agegrps by unexposed, recovered, nil, sick, mild, severe
           therefore, we only use rows 1:4
    =#

