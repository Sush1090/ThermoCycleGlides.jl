abstract type ThermoCycleGlidesData end

struct CompressorData <:ThermoCycleGlidesData 
    fluid::EoSModel
    z::AbstractVector
    p_ratio::Vector{Float64}
    eff::Vector{Float64}
end

struct ExpanderData <:ThermoCycleGlidesData 
    p_ratio::Vector{Float64}
    eff::Vector{Float64}
end


function generate_map(data::ThermoCycleGlidesData,order = 3)
    poly = fit(data.p_ratio, data.eff, order)
    η_fun(p) = clamp(poly(p), 0.0, 1.0)
    return η_fun
end


export CompressorData, ExpanderData, generate_map


function parameter_fitter(data::CompressorData,x0::AbstractVector; T_ref = 300.0, p_ref = 101325.0)
    @assert length(x0) == 2
    π_data = data.p_ratio
    η_data = data.eff

    function objective(x)
        η_design = x[1]
        r_v      = x[2]

        # Penalise non-physical values
        if η_design <= 0 || η_design > 1 || r_v <= 1
            return 1e6
        end

        f = off_design_compressor_relation(
            data.fluid, data.z, η_design, r_v;
            p_ref = p_ref, T_ref = T_ref
        )

        err = 0.0
        @inbounds for i in eachindex(π_data)
            η_model = f(π_data[i])
            err += (η_model - η_data[i])^2
        end

        return err
    end

    # Initial guess
    x0 = x0   # [η_design, built_in_volume_ratio]

    result = Optim.optimize(objective, x0, NelderMead())

    return Optim.minimizer(result), result
end

export parameter_fitter