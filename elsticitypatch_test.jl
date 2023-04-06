
using Revise, ApproxOperator, BenchmarkTools

elements,nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")
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
    Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν)
]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
ops[1].(elements["Ω"];k=k)
ops[2].(elements["Γ"];k=k,f=f)

d = k\f
push!(getfield(nodes[1],:data),:d₁=>(1,d[1:2:2*nₚ-1]))
push!(getfield(nodes[1],:data),:d₂=>(1,d[2:2:2*nₚ]))
Hₑ_PlaneStress = ops[3](elements["Ω"])