"""
Maximize COP for the cycle. 
    Select Var on prob to optimize:
    from heat pump continous variables are : [z, T_evap_in, T_evap_out, T_cond_in, T_cond_out, η_comp, pp_evap, pp_cond, ΔT_sh, ΔT_sc]
    For these [η_comp, pp_evap, pp_cond] are fixed for the cycle.
    For a Heat pump -> target temperature glide for the condensor is a necessity. Hence [T_cond_in, T_cond_out] are fixed. And inlet to evaporator of secondary fluid is fixed.
    Hence we can only optimize [z,T_evap_out,ΔT_sc,ΔT_sh]. z can only be optimized if fluid is a mixture.
"""
function generate_opt_fun(hp::HeatPump, opt_fields::Vector{String})
    all_syms = fieldnames(typeof(hp))
    opt_syms = Symbol.(opt_fields)

    # Validate optimisation field names
    for s in opt_syms
        if !(s in all_syms)
            error("Invalid field name: $s")
        end
    end

    # Return a function that updates via NamedTuple splatting (no Dict)
    function update_fun(x::AbstractVector{T}) where {T<:Real}
        vals = (; zip(all_syms, getfield.(Ref(hp), all_syms))...)
        updated_vals = merge(vals, NamedTuple{opt_syms}(x))
        return HeatPump(; updated_vals...)
    end

    return update_fun
end


function obj(obj_fun_update::Function,x::AbstractVector{T};autodiff::Bool=true,N::Int=20) where {T<:Real}
    prob = obj_fun_update(x)
    sol,res = solve(prob, autodiff = autodiff,N = N)
    if norm(res) > 1e-2
        @warn "The solution did not converge. Residues: $res, norm: $(norm(res))"
    end
    return COP(prob,sol)
end


# function optimize_hp(prob::HeatPump, opt_fields::Vector{String};
#     iterations::Int=100, xtol::Real=1e-8, ftol::Real=1e-8, autodiff::Bool=true, N::Int=20)

#     # Generate the update function for the HeatPump
#     obj_fun_update = generate_opt_fun(prob, opt_fields)

#     # Define the objective function
#     obj_fun(x) = obj(obj_fun_update, x; autodiff=autodiff, N=N)

#     # Initial guess for optimization
#     x0 = zeros(eltype(prob.z), length(opt_fields))

#     # Perform optimization
#     result = optimize(obj_fun, x0, BFGS(); maxiters=iterations, xtol=xtol, ftol=ftol)

#     return result
    
# end