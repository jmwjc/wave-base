using Revise, ApproxOperator, BenchmarkTools, GLMakie
include("input_CG.jl")

ndiv = 100
elements,nodes = import_meshfree("./msh/test_30.msh",ndiv)

nₚ = length(nodes)

s = 1.5/30*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])

I = 1
a = elements["Ω"][1]
𝓖 = collect(a.𝓖)
x = 0.0:1.0/ndiv:1.0
y = 0.0:1.0/ndiv:1.0
𝝭 = zeros(ndiv+1,ndiv+1)
f = Figure()
for i in 1:ndiv+1
   for j in 1:ndiv+1
      n = (i-1)*(ndiv+1) + j
      ξ = 𝓖[n]
      N = ξ[:𝝭]
      𝝭[i,j] = N[n]
   end
end
sufrace(x,y,𝝭)

f