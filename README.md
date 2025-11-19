# ThermoCycleGlides.jl

[![Build Status](https://github.com/Sush1090/ThermoCycleGlides.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Sush1090/ThermoCycleGlides.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sush1090.github.io/ThermoCycleGlides.jl/dev/)

This package aims to solve Heat Pump and ORC systems for given known temperature glides. For now it is robust for sub-critical cycles. 

The thermodynamic computations use Clapeyron.jl. 

# Installation
For the latest release, first type `]` and then:

```julia
pkg> add ThermoCycleGlides
```

For the developer version type:

```julia
pkg> add https://github.com/Sush1090/ThermoCycleGlides.jl
```

# Usage
Usage Heat Pump Example :

```julia
julia> using Clapeyron, ThermoCycleGlides

julia> fluid = cPR(["cyclopentane"],idealmodel = ReidIdeal);

julia> η_comp = 0.75; pp_cond = 2; pp_evap = 2;

julia> T_evap_in = 273.15 + 10; T_evap_out = 273.15 + 0; T_cond_in = 273.15 + 50;  T_cond_out = 273.15+60;

julia> ΔT_sc = 3; ΔT_sh = 10;

julia> hp = HeatPump(fluid=fluid,z=[1.0],T_evap_in=T_evap_in,T_evap_out = T_evap_out,T_cond_in = T_cond_in,T_cond_out=T_cond_out,η_comp=η_comp,pp_evap=pp_evap,pp_cond=pp_cond,ΔT_sc = ΔT_sc,ΔT_sh = ΔT_sh);

julia> sol_hp = solve(hp)
SolutionState{Float64, Int64}([0.12829257763187535, 1.4551588056895837], 16, 7, [0.0, 5.684341886080802e-14], [0.07660159441435545, 0.07660159441435545], [1.6566058479359296, 1.6566058479359296], true, 0, 2.5468671054250572e-15, 8.038873388460929e-14)
 
julia> COP(hp,sol_hp)
-3.735868783526992
```

To plot do the following;

```julia
plot_cycle(hp,sol_hp,N=300)
```

![HP_cyclopentane](docs/src/Images/hp_cyclopentane.png)


# Limitation

1. Fluid models are limited to the ones provided by default in Clapeyron.jl
2. For now the solver is stable for sub-critical parameters. So if incase the solver does converge please check if the parameters provided allow the solution to be subcritical. 
3. For mixtures, it is recommended to use parameters sufficently below the critical point.  
4. If for solving `autodiff = true` then for the first run there will be some compile time. Subsequent runs will be faster. 