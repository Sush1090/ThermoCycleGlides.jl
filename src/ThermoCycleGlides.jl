module ThermoCycleGlides
# Write your package code here.
using Clapeyron, Polynomials, Interpolations, FiniteDifferences
using ForwardDiff, CommonSolve,LinearAlgebra
using Plots


abstract type ThermoCycleProblem end

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
