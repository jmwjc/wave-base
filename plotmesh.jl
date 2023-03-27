using Revise, ApproxOperator, CairoMakie

elements,points = ApproxOperator.importcomsol("圆形骨料.mphtxt")

f = Figure()
ax = Axis(f[1,1])
ax.aspect = AxisAspect(1)

for a in elements["Ω"]
    xs = [p.x for p in a.vertices[[1,2,3,1]]]
    ys = [p.y for p in a.vertices[[1,2,3,1]]]
    lines!(xs,ys,linewidth = 0.5, color = :black)
end

for a in elements["Γ"]
    xs = [p.x for p in a.vertices]
    ys = [p.y for p in a.vertices]
    lines!(xs,ys,linewidth = 1.0, color = :blue)
end

f