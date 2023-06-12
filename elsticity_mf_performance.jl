
using Revise, ApproxOperator, BenchmarkTools, Printf, SparseArrays, Pardiso, TimerOutputs
include("input.jl")
elements,nodes,elms,nds = import_gauss_quadratic("./msh/test_50.msh","./msh/test_50.msh",:TriGI3)
# elements,nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")
# nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")

const to = TimerOutput()
ps = MKLPardisoSolver()
set_matrixtype!(ps,2)

nₚ = length(nodes)
nₑ = length(elements["Ω"])
nₒ = length(nds)
nₑₒ = length(elms["Ω"])

s = 1.5*410/50*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

@timeit to "shape function" begin
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ"])
set𝝭!(elements["Γᵗ"])
set∇𝝭!(elements["Ωₒ"])
end

E = 3e6
ν = 0.3
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
    Operator{:∫vᵢgᵢds}(:α=>1e7*E),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫∫ρvᵢuᵢdxdy}(:ρ=>ρ)
]

@timeit to "assembly matrix" begin
# k = zeros(2*nₚ,2*nₚ)
# m = zeros(2*nₚ,2*nₚ)
# kα = zeros(2*nₚ,2*nₚ)
k = spzeros(2*nₚ,2*nₚ)
m = spzeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
fα = zeros(2*nₚ)

ops[1](elements["Ω"],k)
ops[4](elements["Ω"],m)
ops[2](elements["Γ"],k,fα)
end

d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
push!(nodes,:d₁=>d₁,:d₂=>d₂)

F₀ = 1
Θ = π
β = 0.25
γ = 0.5
𝑓 = 100
force_time = 1/𝑓
Δt = π*force_time/80
total_time = 250*Δt
times = 0.0:Δt:total_time
d = zeros(2nₚ)
v = zeros(2nₚ)
a = zeros(2nₚ)
aₙ = zeros(2nₚ)

u₁ = zeros(nₒ)
u₂ = zeros(nₒ)
σ₁₁ = zeros(nₒ)
σ₂₂ = zeros(nₒ)
σ₁₂ = zeros(nₒ)
for (n,t) in enumerate(times)

    if t ≤ force_time
        prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->F₀*sin(2Θ*𝑓*t))
    else
        prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->0.0)
    end
    fill!(f,0.0)
    ops[3](elements["Γᵗ"],f)

    @timeit to "solve" begin
    # predictor phase
    global d .+= Δt*v + Δt^2/2.0*(1.0-2.0*β)*aₙ
    global v .+= Δt*(1.0-γ)*aₙ
    # a = (m + β*Δt^2*k)\(f+fα-k*d)
    solve!(ps,a,m + β*Δt^2*k,f+fα-k*d)

    # Corrector phase
    global d .+= β*Δt^2*a 
    global v .+= γ*Δt*a
    global aₙ .= a
    end


    @timeit to "output" begin
    d₁ .= d[1:2:2*nₚ]
    d₂ .= d[2:2:2*nₚ]

    fill!(u₁,0.0)
    fill!(u₂,0.0)
    fill!(σ₁₁,0.0)
    fill!(σ₂₂,0.0)
    fill!(σ₁₂,0.0)
    for (j,p) in enumerate(elements["Ωₒ"])
        ξ, = p.𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        for (i,xᵢ) in enumerate(p.𝓒)
            u₁[j] += N[i]*xᵢ.d₁
            u₂[j] += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₁[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₁
        end
        σ₁₁[j] = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂[j] = Cᵢᵢⱼⱼ*ε₁₁+Cᵢᵢᵢᵢ*ε₂₂
        σ₁₂[j] = Cᵢⱼᵢⱼ*ε₁₂
    end

    fo = open("./vtk/50/figure"*string(n,pad=4)*".vtk","w")
    @printf fo "# vtk DataFile Version 2.0\n"
    @printf fo "Test\n"
    @printf fo "ASCII\n"
    @printf fo "DATASET POLYDATA\n"
    @printf fo "POINTS %i float\n" nₒ
    for p in nds
        @printf fo "%f %f %f\n" p.x p.y p.z
    end
    @printf fo "POLYGONS %i %i\n" nₑₒ 4*nₑₒ
    for ap in elms["Ω"]
        𝓒 = ap.vertices
        @printf fo "%i %i %i %i\n" 3 (x.i-1 for x in 𝓒)...
    end
    @printf fo "POINT_DATA %i\n" nₒ
    @printf fo "SCALARS UX float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for i in 1:nₒ
        @printf fo "%f\n" u₁[i]
    end
    @printf fo "SCALARS UY float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for i in 1:nₒ
        @printf fo "%f\n" u₂[i]
    end
    @printf fo "TENSORS STRESS float\n"
    for i in 1:nₒ
        @printf fo "%f %f %f\n" σ₁₁[i] σ₁₂[i] 0.0
        @printf fo "%f %f %f\n" σ₁₂[i] σ₂₂[i] 0.0
        @printf fo "%f %f %f\n" 0.0 0.0 0.0
    end
    close(fo)
    end
end
show(to)