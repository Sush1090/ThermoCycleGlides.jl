abstract type ThermoCycleGlidesData end

struct CompressorData <:ThermoCycleGlidesData 
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