
using  ApproxOperator, LinearAlgebra, Printf

include("input.jl")


err = 0.0
for ap in elements["Ω"]
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        uʰ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            uʰ += N[i]*sin(xᵢ.x)
        end
        u = sin(ξ.x)
        err += (u-uʰ)^2*𝑤
    end
end
err = err^0.5