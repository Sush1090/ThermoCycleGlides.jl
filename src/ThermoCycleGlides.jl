module ThermoCycleGlides
# Write your package code here.
using Clapeyron, Polynomials, Interpolations, FiniteDifferences
using ForwardDiff, CommonSolve, LinearAlgebra
using Plots, StaticArrays,Metaheuristics

import LinearAlgebra: norm, rank

abstract type ThermoCycleProblem end

abstract type ThermoCycleSolution end

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
