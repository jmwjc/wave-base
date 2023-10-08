using Revise, ApproxOperator, BenchmarkTools, GLMakie
# CairoMakie 
include("input_CG.jl")

ndiv = 5
elements,nodes = import_meshfree("./msh/mf_5.msh",ndiv)

nₚ = length(nodes)
s = 2.9/5*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

f = Figure()
ax = Axis(f[1, 1])
set𝝭!(elements["Ω"])
I = 26
a = elements["Ω"][1]
𝓖 = collect(a.𝓖)
x = 0.0:1.0/ndiv:1.0
y = 0.0:1.0/ndiv:1.0
𝝭 = zeros(ndiv+1,ndiv+1)
for i in 1:ndiv+1
   for j in 1:ndiv+1
      n = (i-1)*(ndiv+1) + j
      ξ = 𝓖[n]
      N = ξ[:𝝭]
      𝝭[i,j] = N[I]
   end
end
surface(x,y,𝝭)
# lines!(ax, y, 𝝭[3,:])
#  f

