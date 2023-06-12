function import_gauss_quadratic(filename::String,fileoname::String,s::Symbol)
    elms,nds = ApproxOperator.importmsh(filename)
    elms_output,nds_output = ApproxOperator.importmsh(fileoname)
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
    sp = ApproxOperator.RegularGrid(x,y,z,n=1,γ=2)

    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21
    # parameters = (:Wave2D,:□,:CubicSpline)
    # n𝒑 = 15
    scheme = ApproxOperator.quadraturerule(s)

    elements = Dict([
        "Ω"=>ReproducingKernel{parameters...,:Tri3}[],
        "Ωₒ"=>ReproducingKernel{parameters...,:Poi1}[],
        "Γᵗ"=>ReproducingKernel{parameters...,:Poi1}[],
        "Γ"=>ReproducingKernel{parameters...,:Seg2}[],
    ])

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    ng = length(scheme[:w])
    ns = 0
    nₑ = length(elms["Ω"])

    for (C,a) in enumerate(elms["Ω"])
        indices = Set{Int}()
        for i in 1:ng
            ξ = scheme[:ξ][i]
            η = scheme[:η][i]
            x,y,z = a(ξ,η)
            union!(indices,sp(x,y,z))
        end
        nc = length(indices)
        for i in indices
            push!(𝓒,nodes[i])
        end
        element = ReproducingKernel{parameters...,:Tri3}((c,nc,𝓒),(g,ng,𝓖))
        push!(elements["Ω"],element)

        c += nc
        g += ng
        ns += nc*ng
    end

    data = Dict([
        :ξ=>(1,scheme[:ξ]),
        :η=>(1,scheme[:η]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(g)),
        :y=>(2,zeros(g)),
        :z=>(2,zeros(g)),
        :𝑤=>(2,zeros(g)),
        :𝝭=>(4,zeros(ns)),
        :∂𝝭∂x=>(4,zeros(ns)),
        :∂𝝭∂y=>(4,zeros(ns)),
        :𝗠=>(0,zeros(n𝒑)),
        :∂𝗠∂x=>(0,zeros(n𝒑)),
        :∂𝗠∂y=>(0,zeros(n𝒑)),
    ])
    
    G = 0
    s = 0
    for (C,a) in enumerate(elms["Ω"])
        𝐴 = ApproxOperator.get𝐴(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data)
            ξ = x.ξ
            η = x.η
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐴*x.w
            push!(𝓖,x)
            s += getfield(elements["Ω"][C],:𝓒)[2]
        end
    end
 
    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    ns = 0
    ng = 1
    gauss_scheme = :PoiGI1
    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    nₑ = length(elms["Γᵗ"])
    for (C,a) in enumerate(elms["Γᵗ"])
        indices = Set{Int}()
        for i in 1:ng
            x,y,z = a.x,a.y,a.z
            union!(indices,sp(x,y,z))
        end
        nc = length(indices)
        for i in indices
            push!(𝓒,nodes[i])
        end
        element = ReproducingKernel{parameters...,:Poi1}((c,nc,𝓒),(g,ng,𝓖))
        push!(elements["Γᵗ"],element)
        c += nc
        g += ng
        ns += ng*nc
    end

    G = 0
    s = 0
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,[0.]),
        :y=>(2,[0.]),
        :z=>(2,[0.]),
        :𝑤=>(2,[1.]),
        :𝝭=>(4,zeros(ns)),
        :𝗠=>(0,zeros(n𝒑))
       ])
    for (C,a) in enumerate(elms["Γᵗ"])
        x_ = a.x
        y_ = a.y
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            x.x = x_
            x.y = y_
            x.z = 0.0
            x.𝑤 = 1.0
            push!(𝓖,x)
            s += getfield(elements["Γᵗ"][C],:𝓒)[2]
        end
    end

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    ng = 3
    ns = 0
    gauss_scheme = :SegGI3
    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    nₑ = length(elms["Γ"])

    for (C,a) in enumerate(elms["Γ"])
        indices = Set{Int}()
        for i in 1:ng
            ξ = scheme[:ξ][i]
            x,y,z = a(ξ)
            union!(indices,sp(x,y,z))
        end
        nc = length(indices)
        for i in indices
            push!(𝓒,nodes[i])
        end
        element = ReproducingKernel{parameters...,:Seg2}((c,nc,𝓒),(g,ng,𝓖))
        push!(elements["Γ"],element)
        c += nc
        g += ng
        ns += ng*nc
    end
       
    G = 0
    s = 0
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :n₁=>(3,zeros(nₑ)),
        :n₂=>(3,zeros(nₑ)),
        :𝗠=>(0,zeros(n𝒑)),
        :∂𝗠∂x=>(0,zeros(n𝒑)),
        :∂𝗠∂y=>(0,zeros(n𝒑)),
        :𝝭=>(4,zeros(ns)),
        :∂𝝭∂x=>(4,zeros(ns)),
        :∂𝝭∂y=>(4,zeros(ns))
    ])
    for (C,a) in enumerate(elms["Γ"])
        𝐿 = ApproxOperator.get𝐿(a)
        x₁ = a.vertices[1].x
        x₂ = a.vertices[2].x
        y₁ = a.vertices[1].y
        y₂ = a.vertices[2].y
        n₁ = (y₂-y₁)/𝐿
        n₂ = (x₁-x₂)/𝐿
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            x_,y_,z_ = a(ξ)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐿*x.w/2
            push!(𝓖,x)
            s += getfield(elements["Γ"][C],:𝓒)[2]
        end
        elements["Γ"][C].n₁ = n₁
        elements["Γ"][C].n₂ = n₂
    end

    scheme = ApproxOperator.quadraturerule(:PoiGI1)

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    ng = length(scheme[:w])
    ns = 0
    nₑ = length(nds_output)

    for (C,a) in enumerate(nds_output)
        indices = Set{Int}()
        for i in 1:ng
            x,y,z = a.x,a.y,a.z
            union!(indices,sp(x,y,z))
        end
        nc = length(indices)
        for i in indices
            push!(𝓒,nodes[i])
        end
        element = ReproducingKernel{parameters...,:Poi1}((c,nc,𝓒),(g,ng,𝓖))
        push!(elements["Ωₒ"],element)

        c += nc
        g += ng
        ns += nc*ng
    end

    data = Dict([
        :x=>(2,zeros(g)),
        :y=>(2,zeros(g)),
        :z=>(2,zeros(g)),
        :𝝭=>(4,zeros(ns)),
        :∂𝝭∂x=>(4,zeros(ns)),
        :∂𝝭∂y=>(4,zeros(ns)),
        :𝗠=>(0,zeros(n𝒑)),
        :∂𝗠∂x=>(0,zeros(n𝒑)),
        :∂𝗠∂y=>(0,zeros(n𝒑)),
    ])
    
    G = 0
    s = 0
    for (C,a) in enumerate(nds_output)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data)
            x.x = a.x
            x.y = a.y
            x.z = a.z
            push!(𝓖,x)
            s += getfield(elements["Ωₒ"][C],:𝓒)[2]
        end
    end

    return elements,nodes,elms_output,nds_output
end