module ThermoCycleGlides
# Write your package code here.
using Clapeyron, Polynomials, Interpolations, FiniteDifferences
using ForwardDiff, Optim, Metaheuristics, CommonSolve,LinearAlgebra
using Plots


#NonlinearSolver - NR
include("NonlinearSolver/newton-raphson.jl")

# Thermo-fixes
include("thermoextensions/utils.jl")
include("thermoextensions/fix_instabilites.jl")

# Cycle Structs
include("CycleStructs/ORC.jl")
include("CycleStructs/HeatPump.jl")

#solve
include("Solve/solve.jl")

#opt
include("Optimizations/Optimize.jl")


#Plotting
include("Plots/cycleplots.jl")

end
