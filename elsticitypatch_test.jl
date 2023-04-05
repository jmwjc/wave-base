
using Revise, ApproxOperator, BenchmarkTools

elements,nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")
# nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")

nₚ = length(nodes)

set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γ"])

E = 3e6
v=0.3
D=12
I=D^3/12
EI=E*I
u(x,y) = x+y
ApproxOperator.prescribe!(elements["Γ"],:g=>(x,y,z)->u(x,y))
ApproxOperator.prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))

ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:k=>1.0),
    Operator{:∫vᵢgᵢds}(:α=>1e7*E),
    Operator{:Hₑ_PlaneStress}()
]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
ops[1].(elements["Ω"];k=k)
ops[2].(elements["Γ"];k=k,f=f)

d = k\f
push!(getfield(nodes[1],:data),:d=>(1,d))
Hₑ_PlaneStress = ops[3](elements["Ω"])