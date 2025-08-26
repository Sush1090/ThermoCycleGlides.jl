

function fix_nan!(x::AbstractVector)
    isnan_mask = isnan.(x)
    if all(!isnan, x)
        return x  # no NaNs to fix
    end

    # Indices of non-NaN values
    idxs = findall(!isnan, x)
    values = x[idxs]

    # Interpolation over the full index range
    interp = LinearInterpolation(idxs, values, extrapolation_bc=Line())

    for i in eachindex(x)
        if isnan(x[i])
            x[i] = interp(i)
        end
    end

    return x
end

function fix_nan(x::AbstractVector)
    y = copy(x)
    fix_nan!(y)
    return y
end

function fix_nan_poly!(x::AbstractVector; degree::Int=3)
    n = length(x)
    isnan_mask = isnan.(x)

    if all(!isnan, x)
        return x  # nothing to fix
    end

    # Known values
    known_idxs = findall(!isnan, x)
    known_vals = x[known_idxs]

    if length(known_idxs) <= degree
        error("Not enough non-NaN points for polynomial interpolation of degree $degree")
    end

    # Fit polynomial to known data
    p = fit(known_idxs, known_vals, degree)

    # Replace NaNs
    for i in eachindex(x)
        if isnan(x[i])
            x[i] = p(i)
        end
    end

    return x
end

function fix_nan_poly(x::AbstractVector; degree::Int=3)
    y = copy(x)
    fix_nan_poly!(y; degree=degree)
    return y
end

