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
export TAB1OUT, TAB1OUT2, TAB2OUT, TAB3OUT

##############################################################################
##
## Load files
##
##############################################################################
include("helper.jl")
include("tab.jl")


end