module BilevelBenchmarkTools

#Personal dependencies
using BOLIB

include("compute_lambda.jl")
include("referee_strategies.jl")
include("agregate_efforts.jl")

end
