function import_meshfree(filename::String,ndiv::Int)
    ~,nds = ApproxOperator.importmsh(filename)
    nₚ = length(nds)
    nodes = Node{(:𝐼,),1}[]
    x = zeros(nₚ)
    y = zeros(nₚ)
    z = zeros(nₚ)
    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    for (i,p) in enumerate(nds)
        node = Node{(:𝐼,),1}((i,),data)
        node.x = p.x
        node.y = p.y
        node.z = p.z
        push!(nodes,node)
    end
    # parameters = (:Wave2D,:□,:CubicSpline)
    parameters = (:Wave2D,:□,:Gaussian)
    n𝒑 = 15

    elements = Dict([
        "Ω"=>ReproducingKernel{parameters...,:Tri3}[],
    ])

    𝓒 = nodes
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]

    ng = (ndiv+1)^2
    ns = length(nds)
    data = Dict([
        :x=>(2,zeros(ng)),
        :y=>(2,zeros(ng)),
        :z=>(2,zeros(ng)),
        :𝝭=>(4,zeros(ng*ns)),
        :∂𝝭∂x=>(4,zeros(ng*ns)),
        :∂𝝭∂y=>(4,zeros(ng*ns)),
        :𝗠=>(0,zeros(n𝒑)),
        :∂𝗠∂x=>(0,zeros(n𝒑)),
        :∂𝗠∂y=>(0,zeros(n𝒑)),
    ])
    G = 0
    s = 0
    for x in 0.0:1.0/ndiv:1.0
        for y in 0.0:1.0/ndiv:1.0
            G += 1
            ξ = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((1,G,1,s),data)
            ξ.x = x
            ξ.y = y
            ξ.z = 0.0
            push!(𝓖,ξ)
            s += nₚ
        end
    end

    element = ReproducingKernel{parameters...,:Tri3}((0,nₚ,𝓒),(0,ng,𝓖))
    push!(elements["Ω"],element)

    return elements, nodes
end