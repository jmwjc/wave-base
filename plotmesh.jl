using Revise, ApproxOperator, CairoMakie

elements,points,entities = ApproxOperator.importcomsol("圆形骨料.mphtxt")

f = Figure()
ax = Axis(f[1,1])
ax.aspect = AxisAspect(1)

for a in elements["Ω"]
    xs = [p.x for p in a.vertices[[1,2,3,1]]]
    ys = [p.y for p in a.vertices[[1,2,3,1]]]
    lines!(xs,ys,linewidth = 0.5, color = :black)
end

index = findall(x->x∈(0,1,2,6,8),entities["Γ"])
# index = findall(x->x∈(14,),entities["Γ"])
for a in elements["Γ"][index]
    xs = [p.x for p in a.vertices]
    ys = [p.y for p in a.vertices]
    lines!(xs,ys,linewidth = 1.0, color = :blue)
end

f