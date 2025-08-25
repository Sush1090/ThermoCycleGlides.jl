function box_projection(x::Array{T,1},lb::Array{TT,1},ub::Array{TTT,1}) where {T <: Real, TT <: Real, TTT <: Real}
    y = copy(x)
    for i in eachindex(x)
        if y[i] < lb[i]
            y[i] = lb[i]
        elseif y[i] > ub[i]
            y[i] = ub[i]
        end
    end
    return y
end

function box_projection!(x::Array{T,1},lb::Array{TT,1},ub::Array{TTT,1}) where {T <: Real, TT <: Real, TTT <: Real}
    for i in eachindex(x)
        if x[i] < lb[i]
            x[i] = lb[i]
        elseif x[i] > ub[i]
            x[i] = ub[i]
        end
    end
end


function constrained_newton_fd(f::Function,x::Array{T,1},
    lb::Array{TT,1},ub::Array{TTT,1};xtol::TOL=1e-8,ftol::TOL=1e-8,iterations::S=100,
    fd_method::M = central_fdm(2, 1)) where  {T <: Real, S <: Integer, TT <: Real, TTT <: Real, M <: FiniteDifferences.AdaptedFiniteDifferenceMethod,TOL<:Real}

    n = length(x)
    xk = copy(x)
    xn = similar(x)
    jk = Array{T,2}(undef,n,n)

    lenx = zero(T)
    lenf = zero(T)

    iter = 0
    while true
        jk .= FiniteDifferences.jacobian(fd_method,f,xk)[1]
   
        if !all(isfinite,jk)
            error("The jacobian has non-finite elements")
        end
           
        if rank(jk) == n # Compute a Newton step
            xn .= xk - jk\f(xk)
        else # Compute a Levenberg-Marquardt step
            λ = 0.001*norm(f(xk))^2
         
            xn .= xk - (jk'jk + λ*I)\(jk'f(xk))
        end

        box_projection!(xn,lb,ub)
        # @show xn

        lenx = maximum(abs,norm((xn.-xk))/norm(xk))
        lenf = maximum(abs,norm((f(xn).-f(xk)))/norm(f(xk)))
        # @show lenx,lenf
    
        xk .= xn

        iter += 1

        if iter >= iterations || (lenx <= xtol || lenf <= ftol)
            break
        end
  
    end
  
    return xk

end

function constrained_newton_ad(f::Function,x::Array{T,1},
    lb::Array{TT,1},ub::Array{TTT,1};xtol::TOL=1e-8,ftol::TOL=1e-8,iterations::S=100) where 
    {T <: Real, S <: Integer, TT <: Real, TTT <: Real,TOL<:Real}

    type_promoted = promote_type(eltype(x), eltype(lb), eltype(ub))
    xk = copy(x)
    n = length(x)
    xn = similar(x)
    jk = similar(x, type_promoted, n, n)
    # @show x, jn, 
    for iter in 1:iterations
        jk = ForwardDiff.jacobian(f, xk)

        fx = f(xk)
        # @show jk
        # if cond(jk) < 1e10
            xn = xk - jk \ fx
        # else
            # λ = 1e-3
            # xn .= xk - (jk' * jk + λ * I) \ (jk' * fx)
        # end

        # simple box projection (not differentiable, but works if outside AD path)
        for i in eachindex(xn)
            xn[i] = clamp(xn[i], lb[i], ub[i])
        end

        if norm(xn - xk) < xtol || norm(f(xn)) < ftol
            break
        end
        lenx = maximum(abs,norm((xn.-xk))/norm(xk))
        lenf = maximum(abs,norm((f(xn).-f(xk)))/norm(f(xk)))
        if iter >= iterations || (lenx <= xtol || lenf <= ftol)
            break
        end

        xk = xn
    end

    return xk
end

