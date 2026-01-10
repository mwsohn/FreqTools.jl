module FreqTools

################################################################################
##
## Dependencies
##
################################################################################

using DataFrames, NamedArrays, Printf, PrettyTables,
    Reexport, Statistics, HypothesisTests, StatsAPI, Tables

@reexport using FreqTables

##############################################################################
##
## Exported methods and types (in addition to everything reexported above)
##
##############################################################################

export tab

##############################################################################
##
## Load files
##
##############################################################################
include("tab.jl")

end