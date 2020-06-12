abstract type AbstractChebyshevBasis{P} <: AbstractGegenbauerBasis{P} end

polynomial_type(::Type{<:AbstractChebyshevBasis}, V::Type) = MP.polynomialtype(V, Float64)

reccurence_first_coef(::Type{<:AbstractChebyshevBasis}, degree) = 2
reccurence_third_coef(::Type{<:AbstractChebyshevBasis}, degree) = -1
reccurence_deno_coef(::Type{<:AbstractChebyshevBasis}, degree) = 1

"""
    struct ChebyshevBasisFirstKind{P} <: AbstractChebyshevBasis{P}
        elements::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\frac{1}{\\sqrt{1 - x^2}}`` over the interval ``[-1, 1]``.
"""
struct ChebyshevBasisFirstKind{P} <: AbstractChebyshevBasis{P}
    elements::Vector{P}
end

const ChebyshevBasis{P} = ChebyshevBasisFirstKind{P}

degree_one_univariate_polynomial(::Type{<:ChebyshevBasisFirstKind}, variable::MP.AbstractVariable) = MA.@rewrite(variable + 0)

function scalar_product_function(::Type{<:ChebyshevBasisFirstKind})
    function sp(i::Int)
        if i == 0
            return π
        elseif mod(i, 2) == 1
            return 0
        else
            n = Int(i/2)
            return π*factorial(i)/(2^i * factorial(n)^2)
        end
    end
    return sp
end

"""
    struct ChebyshevBasisSecondKind{P} <: AbstractChebyshevBasis{P}
        elements::Vector{P}
    end

Orthogonal polynomial with respect to the univariate weight function ``w(x) = \\sqrt{1 - x^2}`` over the interval ``[-1, 1]``.
"""
struct ChebyshevBasisSecondKind{P} <: AbstractChebyshevBasis{P}
    elements::Vector{P}
end

degree_one_univariate_polynomial(::Type{<:ChebyshevBasisSecondKind}, variable::MP.AbstractVariable) = MA.@rewrite(2variable + 0)


function scalar_product_function(::Type{<:ChebyshevBasisSecondKind})
    function sp(i::Int)
        if i == 0
            return π/2
        elseif mod(i, 2) == 1
            return 0
        else
            n = Int(i/2)
            return π*factorial(i)/(2^(i+1) * factorial(n) * factorial(n +1))
        end
    end
    return sp
end

