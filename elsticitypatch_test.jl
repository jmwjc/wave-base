
using Revise, ApproxOperator, BenchmarkTools, YAML

elements,nodes = ApproxOperator.importmsh_fem("./msh/test.msh")
# elements,nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")
# nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")

nₚ = length(nodes)

set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γ"])

E = 3e6
ν=0.3
u(x,y) = x+y
v(x,y) = x+y
∂u∂x(x,y) = 1.0
∂u∂y(x,y) = 1.0
∂v∂x(x,y) = 1.0
∂v∂y(x,y) = 1.0
ApproxOperator.prescribe!(elements["Γ"],:g₁=>(x,y,z)->u(x,y))
ApproxOperator.prescribe!(elements["Γ"],:g₂=>(x,y,z)->v(x,y))
ApproxOperator.prescribe!(elements["Γ"],:n₁₁=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Γ"],:n₁₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γ"],:n₂₂=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
ApproxOperator.prescribe!(elements["Ω"],:v=>(x,y,z)->v(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->∂u∂x(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->∂u∂y(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->∂v∂x(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->∂v∂y(x,y))

ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢgᵢds}(:α=>1e13*E),
    Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
    Operator{:∫wVdΓ}
]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
ops[1].(elements["Ω"];k=k)
ops[2].(elements["Γ"];k=k,f=f)

d = k\f
push!(getfield(nodes[1],:data),:d₁=>(1,d[1:2:2*nₚ-1]))
push!(getfield(nodes[1],:data),:d₂=>(1,d[2:2:2*nₚ]))
Hₑ_PlaneStress = ops[3](elements["Ω"])

# Θ = π
# β = 0.25
# γ = 0.5
# Δt = 0.01
# 𝑓 = 100
# total_time = 10.0
# times = 0.0:Δt:total_time
# d = zeros(nₚ)
# x = zeros(length(times))
# deflection = zeros(length(times))
# dexact = zeros(length(times))
# v = zeros(nₚ)
# aₙ = zeros(nₚ)
# for (n,t) in enumerate(times)

#     prescribe!(elements["Γ"],:V=>(x,y,z)->F₀*sin(2Θ*𝑓*t))   
                       
#     fₙ = zeros(nₚ)
#     ops[4](elements["Γ"],fₙ)

#     # predictor phase
#     d .+= Δt*v + Δt^2/2.0*(1.0-2.0*β)*aₙ
#     v .+= Δt*(1.0-γ)*aₙ
#     a = (m + β*Δt^2*(k+kα))\(fₙ+fα-(k+kα)*d)
#     # Corrector phase
#     d .+= β*Δt^2*a
#     v .+= γ*Δt*a
#     aₙ .= a

#     # cal deflection
#     ξ = elements["Γ"][1].𝓖[1]
#     N = ξ[:𝝭]
#     for (i,xᵢ) in enumerate(elements["Γ"][1].𝓒)
#         I = xᵢ.𝐼
#         deflection[n] += N[i]*d[I]
#     end 

#     # cal exact solution
#     dexact[n] = w(5.0,t)

# end

# f = Figure()
# ax = Axis(f[1,1])

# scatterlines!(times,deflection)
# lines!(times,dexact)

# f