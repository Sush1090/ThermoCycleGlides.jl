module ThermoCycleGlides

using Clapeyron, Polynomials, Interpolations, FiniteDifferences
using ForwardDiff, CommonSolve, LinearAlgebra
using Plots, StaticArrays,Metaheuristics

import LinearAlgebra: norm, rank
import Base: show, length, copy, similar, promote_type



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

#Optimization
include("Optimizations/Optimize.jl")


#Plotting
include("Plots/cycleplots.jl")


# PrecompileTools
include("precompile.jl")

    function show(prob::ThermoCycleProblem)
        show_parameters(prob)
    end
    function show(sol::SolutionState)
        show_parameters(sol)
    end

    export show
end
