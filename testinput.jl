
using Revise, ApproxOperator

elements,nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")
# nodes = ApproxOperator.importcomsol_fem("圆形骨料.mphtxt")

set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γ"])