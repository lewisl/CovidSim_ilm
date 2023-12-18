
A model is a named tuple of all of the data structures that are created and initialized before a
simulation is run. You'll never use these types explicitly. The data structures are created by 
function setup_model and functions it calls. Data for parameters are loaded from yaml files. The
content of the parameter files are converted to appropriate Juia data structures.

A model contains these elements of the type shown:
    :ndays              => Int64
    :day1               => Date
    :locales            => Vector{Int64} 
    :dat                => NamedTuple{(:popdat, :agegrp_idx)}
        :popdat         => Dict{Int64, LazyTable} # the int is a locale identifier
        :LazyTable      => contains columns:
                           :status      => Symbol
                           :agegrp      => Symbol
                           :cond        => Symbol
                           :duration    => Int64 
                           :variant     => Vector{Symbol}
                           :sickday     => Vector{Symbol}
                           :recovday    => Vector{Int64}
                           :deadday     => Vector{Int64} 
                           :ring        => Int64
                           :sdcase      => Int64 ???
                           :vaxstatus   => Symbol 
                           :vaxrcvd     => Vector{Symbol} 
                           :vaxday      => Vector{Int64}
                           :tested      => Bool
                           :testday     => Int64
                           :quar        => Bool
                           :quarday     => Int64
    :series, 
    :geo, 
    :progressionset, 
    :vaxset, 
    :vaxschedset, 
    :infectset, 
    :social, 
    :trvec, 
    :variantlist, 
    :vaxlist, 
    :indoor_seq, 
    :seriescolnames
