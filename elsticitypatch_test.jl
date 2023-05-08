
using Revise, ApproxOperator, BenchmarkTools, Printf
include("importmsh.jl")
elements,nodes = import_fem("./msh/test_4.msh")
# elements,nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")
# nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")

nₚ = length(nodes)
nₑ = length(elements["Ω"])

set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γ"])
set𝝭!.(elements["Γᵗ"])
E = 3e6
ν=0.3
Cᵢᵢᵢᵢ = E/(1-ν^2)
Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
Cᵢⱼᵢⱼ = E/2/(1+ν)

ρ = 1.0

prescribe!(elements["Γ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:n₁₁=>(x,y,z,n₁,n₂)->n₁*n₁)
prescribe!(elements["Γ"],:n₁₂=>(x,y,z,n₁,n₂)->n₁*n₂)
prescribe!(elements["Γ"],:n₂₂=>(x,y,z,n₁,n₂)->n₂*n₂)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)                 

ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢgᵢds}(:α=>1e13*E),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫∫ρvᵢuᵢdxdy}(:ρ=>ρ)
]

k = zeros(2*nₚ,2*nₚ)
m = zeros(2*nₚ,2*nₚ)
kα = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
fα = zeros(2*nₚ)
ops[1](elements["Ω"],k)
ops[4](elements["Ω"],m)
ops[2](elements["Γ"],kα,fα)

d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
push!(nodes,:d₁=>d₁,:d₂=>d₂)

F₀ = 1e-7
Θ = π
β = 0.25
γ = 0.5
𝑓 = 100000
force_time = 1/𝑓
Δt = force_time/10
total_time = 100*Δt
times = 0.0:Δt:total_time
d = zeros(2nₚ)
v = zeros(2nₚ)
aₙ = zeros(2nₚ)
for (n,t) in enumerate(times)

    if t ≤ force_time
        prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->F₀*sin(2Θ*𝑓*t))
    else
        prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->0.0)
    end
    fill!(f,0.0)
    ops[3](elements["Γᵗ"],f)

    # predictor phase
    global d .+= Δt*v + Δt^2/2.0*(1.0-2.0*β)*aₙ
    global v .+= Δt*(1.0-γ)*aₙ
    a = (m + β*Δt^2*(k+kα))\(f+fα-(k+kα)*d)
    # Corrector phase
    global d .+= β*Δt^2*a 
    global v .+= γ*Δt*a
    global aₙ .= a


    d₁ .= d[1:2:2*nₚ]
    d₂ .= d[2:2:2*nₚ]

    fo = open("./vtk/figure"*string(n,pad=4)*".vtk","w")
    @printf fo "# vtk DataFile Version 2.0\n"
    @printf fo "Test\n"
    @printf fo "ASCII\n"
    @printf fo "DATASET POLYDATA\n"
    @printf fo "POINTS %i float\n" nₚ
    for p in nodes
        @printf fo "%f %f %f\n" p.x p.y p.z
    end
    @printf fo "POLYGONS %i %i\n" nₑ 4*nₑ
    for ap in elements["Ω"]
        𝓒 = ap.𝓒
        @printf fo "%i %i %i %i\n" 3 (x.𝐼 for x in 𝓒)...
    end
    @printf fo "POINT_DATA %i\n" nₚ
    @printf fo "SCALARS UX float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for p in nodes
        @printf fo "%f\n" p.d₁
    end
    @printf fo "SCALARS UY float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for p in nodes
        @printf fo "%f\n" p.d₂
    end
    @printf fo "CELL_DATA %i\n" nₑ
    @printf fo "TENSORS STRESS float\n"
    for ap in elements["Ω"]
        𝓒 = ap.𝓒
        𝓖 = ap.𝓖
        for (i,ξ) in enumerate(𝓖)
            if i == 1
                B₁ = ξ[:∂𝝭∂x]
                B₂ = ξ[:∂𝝭∂y]
                ε₁₁ = 0.0
                ε₂₂ = 0.0
                ε₁₂ = 0.0
                for (j,xⱼ) in enumerate(𝓒)
                    ε₁₁ += B₁[j]*xⱼ.d₁
                    ε₂₂ += B₂[j]*xⱼ.d₂
                    ε₁₂ += B₁[j]*xⱼ.d₂ + B₂[j]*xⱼ.d₁
                end
                σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
                σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁+Cᵢᵢᵢᵢ*ε₂₂
                σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
                @printf fo "%f %f %f\n" σ₁₁ σ₁₂ 0.0
                @printf fo "%f %f %f\n" σ₁₂ σ₂₂ 0.0
                @printf fo "%f %f %f\n" 0.0 0.0 0.0
                break
            end
        end
    end
    close(fo)
end