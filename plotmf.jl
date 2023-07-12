using Revise, ApproxOperator, BenchmarkTools, GLMakie
include("input_CG.jl")

elements,nodes = import_meshfree("./msh/test_30.msh",30)

nₚ = length(nodes)

s = 1.5*410/60*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!.(elements["Ω"])

f = Figure()

for elm in elements["Ω"]
   xs = [x.x for x in elm.𝓖]
   ys = [y.y for y in elm.𝓖]
   zs = []
   print(zs)
   surface(xs, ys, zs, axis=(type=Axis3,))
end
f