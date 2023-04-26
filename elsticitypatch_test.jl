
using Revise, ApproxOperator, BenchmarkTools
include("D:/wave-base/importmsh.jl")
elements,nodes = importmsh_fem("./msh/test.msh")
# elements,nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")
# nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")

nₚ = length(nodes)

set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γ"])
set𝝭!.(elements["Γᵗ"])
E = 3e6
ν=0.3
u(x,y) = 0.0
v(x,y) = 0.0

ApproxOperator.prescribe!(elements["Γ"],:g₁=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γ"],:g₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γ"],:n₁₁=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Γ"],:n₁₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γ"],:n₂₂=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Ω"],:u=>(x,y,z)->u(x,y))
ApproxOperator.prescribe!(elements["Ω"],:v=>(x,y,z)->v(x,y))
ApproxOperator.prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Ω"],:∂v∂x=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Ω"],:∂v∂y=>(x,y,z)->0.0)

ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢgᵢds}(:α=>1e13*E),
    Operator{:Hₑ_PlaneStress}(:E=>E,:ν=>ν),
    Operator{:∫vᵢtᵢds}()
]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
ops[1].(elements["Ω"];k=k)
ops[2].(elements["Γ"];k=k,f=f)

d = k\f
push!(getfield(nodes[1],:data),:d₁=>(1,d[1:2:2*nₚ-1]))
push!(getfield(nodes[1],:data),:d₂=>(1,d[2:2:2*nₚ]))
Hₑ_PlaneStress = ops[3](elements["Ω"])

Θ = π
β = 0.25
γ = 0.5
Δt = 0.01
𝑓 = 100
total_time = 10.0
times = 0.0:Δt:total_time
d = zeros(2nₚ)
x = zeros(length(times))
deflection = zeros(length(times))
dexact = zeros(length(times))
v = zeros(2nₚ)
aₙ = zeros(2nₚ)
for (n,t) in enumerate(times)

    prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->F₀*sin(2Θ*𝑓*t))                 
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->F₀*sin(2Θ*𝑓*t))
    fₙ = zeros(2nₚ)
    ops[4](elements["Γᵗ"],fₙ)

    # predictor phase
    d .+= Δt*v + Δt^2/2.0*(1.0-2.0*β)*aₙ
    v .+= Δt*(1.0-γ)*aₙ
    a = (m + β*Δt^2*(k+kα))\(fₙ+fα-(k+kα)*d)
    # Corrector phase
    d .+= β*Δt^2*a 
    v .+= γ*Δt*a
    aₙ .= a

end

f = Figure()
ax = Axis(f[1,1])

scatterlines!(times,deflection)
lines!(times,dexact)

f