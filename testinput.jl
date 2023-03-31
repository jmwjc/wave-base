
using Revise, ApproxOperator, BenchmarkTools

elements,nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")
# nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")

nₚ = length(nodes)

set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γ"])

u(x,y) = x+y
ApproxOperator.prescribe!(elements["Γ"],:g=>(x,y,z)->u(x,y))
ApproxOperator.prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))

ops = [
    Operator{:∫∫∇v∇udxdy}(:k=>1.0),
    Operator{:∫vgdΓ}(:α=>1e13),
    Operator{:L₂}()
]

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
ops[1].(elements["Ω"];k=k)
ops[2].(elements["Γ"];k=k,f=f)

d = k\f
push!(getfield(nodes[1],:data),:d=>(1,d))
L₂ = ops[3](elements["Ω"])