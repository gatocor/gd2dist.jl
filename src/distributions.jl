
norm(x,μ,σ2) = exp(-(x-μ)^2/(2*σ2))/sqrt(2*pi*σ2)

lognorm(x,μ,σ2) = -(x-μ)^2/(2*σ2)-log(2*pi)/2-log(σ2)/2

logσ2(σ2,x2,n) = if σ2 > 0; -x2/(2*σ2)-n*log(2*pi)/2-n*log(σ2)/2 else -Inf end

logσ2(σ2,x2,n,σ2₂) = if σ2 > 0; -x2/(2*(σ2+σ2₂))-n*log(2*pi)/2-n*log(σ2+σ2₂)/2 else -Inf end

function logσ2sum(σ2,x2,n,σ2₂)
    
    t = 0
    if length(σ2₂) != length(n)
        for i in 1:length(σ2₂)
            t += logσ2(σ2,x2[i],n[i],σ2₂[i])
        end
        t += logσ2(σ2,x2[end],n[end])
    else
        for i in 1:length(σ2₂)
            t += logσ2(σ2,x2[i],n[i],σ2₂[i])
        end    
    end
    
    return t
end

gamma(x,θ,κ) = x^(κ-1)*exp(-x/θ)/(SpecialFunctions.gamma(κ)*θ^κ)

function gamma(x,θ₁,κ₁,θ₂,κ₂)

    m = κ₁*θ₁+κ₂*θ₂
    v = κ₁*θ₁^2+κ₂*θ₂^2
    θ = v/m
    κ = m/θ

    return gamma(x,θ,κ)
end

loggamma(x,θ,κ) = (κ-1)*log(x)-x/θ-SpecialFunctions.loggamma(κ)-κ*log(θ)

function loggamma(x,θ₁,κ₁,θ₂,κ₂)

    m = κ₁*θ₁+κ₂*θ₂
    v = κ₁*θ₁^2+κ₂*θ₂^2
    θ = v/m
    κ = m/θ

    return loggamma(x,θ,κ)
end

function loggammaθ(θ₁,xsum,xlog,n,κ₁,θ₂,κ₂)

    m = κ₁*θ₁+κ₂*θ₂
    v = κ₁*θ₁^2+κ₂*θ₂^2
    θ = v/m
    κ = m/θ

    return (κ-1)*xlog-xsum/θ-n*SpecialFunctions.loggamma(κ)-n*κ*log(θ)

end

function loggammaκ(κ₁,xsum,xlog,n,θ₁,θ₂,κ₂)

    m = κ₁*θ₁+κ₂*θ₂
    v = κ₁*θ₁^2+κ₂*θ₂^2
    θ = v/m
    κ = m/θ

    return (κ-1)*xlog-xsum/θ-n*SpecialFunctions.loggamma(κ)-n*κ*log(θ)

end

function loggammaθsum(θ₁,xsum,xlog,n,κ₁,θ₂,κ₂)

    t = 0
    if length(xsum) != length(κ₁)
        for i in 1:length(κ₂)
            t += loggammaθ(θ₁,xsum[i],xlog[i],n,κ₁,θ₂[i],κ₂[i])
        end
        t += loggammaθ(θ₁,xsum[end],xlog[end],n,κ₁,0,0)
    else
        for i in 1:length(κ₂)
            t += loggammaθ(θ₁,xsum[i],xlog[i],n,κ₁,θ₂[i],κ₂[i])
        end
    end

    return t

end

function loggammaκsum(κ₁,xsum,xlog,θ₁,θ₂,κ₂)

    t = 0
    if length(xsum) != length(κ₁)
        for i in 1:length(κ₂)
            t += loggammaθ(κ₁,xsum[i],xlog[i],n,θ₁,θ₂[i],κ₂[i])
        end
        t += loggammaθ(κ₁,xsum[end],xlog[end],n,θ₁,0,0)
    else
        for i in 1:length(κ₂)
            t += loggammaθ(κ₁,xsum[i],xlog[i],n,θ₁,θ₂[i],κ₂[i])
        end
    end

    return t

end