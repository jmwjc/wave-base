
"""
importmsh
"""
function importmsh(filename::String)
    fid = open(filename,"r")
    readline(fid)
    line = readline(fid)
    v_,f_,d_ = split(line," ")
    version = parse(Float64,v_)
    filetype = parse(Int,f_)
    datasize = parse(Int,d_)
    readline(fid)
    if version == 4.1
        elements,nodes,entities = import_msh_4(fid)
    elseif version == 2.2
        elements,nodes,entities = import_msh_2(fid)
    else
        println("Version does not match!")
    end
    return elements, nodes, entities
end

function import_msh_4(fid::IO) end

function import_msh_2(fid::IO)
    etype = Dict(1=>:Seg2,2=>:Tri3,3=>:Quad4,8=>:Seg3,9=>:Tri6,15=>:Point)
    points = Point[]
    elements = Dict{String,Any}()
    entities = Dict{String,Any}()
    physicalnames = Dict{Int,String}()
    for line in eachline(fid)
        if line == "\$PhysicalNames"
            numPhysicalNames = parse(Int,readline(fid))
            physicalnames=>Dict{Int,String}()
            for i in 1:numPhysicalNames
                line = readline(fid)
                d_,p_,n_ = split(line," ")
                dimension = parse(Int,d_)
                physicalTag = parse(Int,p_)
                name = strip(n_,'\"')
                physicalnames[physicalTag] = name
            end
            readline(fid)
        elseif line == "\$Nodes"
            line = readline(fid)
            nₚ = parse(Int,line)
            for i in 1:nₚ
                line = readline(fid)
                i,x,y,z = split(line," ")
                i = parse(Int,i)
                x = parse(Float64,x)
                y = parse(Float64,y)
                z = parse(Float64,z)
                push!(points,Point(i,x,y,z))
            end
            readline(fid)
        elseif line == "\$Elements"
            line = readline(fid)
            nₑ = parse(Int,line)
            for i in 1:nₑ
                line = readline(fid)
                entries = split(line," ")
                elmN_ = entries[1]
                elmT_ = entries[2]
                numT_ = entries[3]
                phyT_ = entries[4]
                elmE_ = entries[5]
                l_ = entries[6:end]
                elmNumber = parse(Int,elmN_)
                elmType = parse(Int,elmT_)
                numTag = parse(Int,numT_)
                phyTag = parse(Int,phyT_)
                elmEntary = parse(Int,elmE_)
                nodeList = parse.(Int,l_)
                name = physicalnames[phyTag]
                type = eval(etype[elmType])
                if ~haskey(elements,name)
                    elements[name] = type[]
                    entities[name] = Int[]
                end
                if type == Point
                   push!(elements[name],points[nodeList...])
                else
                   push!(elements[name],type(Tuple(points[i] for i in nodeList)))
                end
                push!(entities[name],elmEntary)
            end
        end
    end
    return elements, points, entities
end

function importmsh_fem(filename::String)
    elms,nds = importmsh(filename)
    nₚ = length(nds)
    nodes = Node{(:𝐼,),1}[]
    data = Dict([:x=>(1,zeros(nₚ)),:y=>(1,zeros(nₚ)),:z=>(1,zeros(nₚ))])
    for (i,p) in enumerate(nds)
        node = Node{(:𝐼,),1}((i,),data)
        node.x = p.x
        node.y = p.y
        node.z = p.z
        push!(nodes,node)
    end

    elements = Dict(["Ω"=>Element{:Tri3}[],"Γ"=>Element{:Seg2}[],"Γᵗ"=>Element{:Poi1}[]])

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 3
    gauss_scheme = :TriGI3
    nₑ = length(elms["Ω"])

    scheme = quadraturerule(gauss_scheme)
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :η=>(1,scheme[:η]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*3)),
        :∂𝝭∂x=>(4,zeros(ng*nₑ*3)),
        :∂𝝭∂y=>(4,zeros(ng*nₑ*3)),
    ])
    for (C,a) in enumerate(elms["Ω"])
        element = Element{:Tri3}((c,3,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 3

        𝐴 = get𝐴(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            η = x.η
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐴*x.w
            push!(𝓖,x)
            s += 3
        end
        g += ng
        push!(elements["Ω"],element)
    end
    
    
    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 2 
    gauss_scheme = :SegGI2
    scheme = quadraturerule(gauss_scheme)
    nₑ = length(elms["Γ"])

    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*2)),
    ])
    for (C,a) in enumerate(elms["Γ"])
        element = Element{:Seg2}((c,2,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 2
       
        𝐿 = get𝐿(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            x_,y_,z_ = a(ξ)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐿*x.w
            push!(𝓖,x)
            s += 2
        end
        g += ng
        push!(elements["Γ"],element)
    end


    data = Dict([:x=>(1,[0]),:y=>(1,[-205]),:z=>(1,[0])])
    𝓒 = [nodes[5]]
    c = 0
    g = 0
    G = 0
    s = 0
    gauss_scheme = :PoiGI1
    nₑ = length(elms["Γᵗ"])
    data_𝓖 = Dict([
         :ξ=>(1,scheme[:ξ]),
         :w=>(1,scheme[:w]),
         :x=>(2,[0.]),
         :y=>(2,[-205]),
         :z=>(2,[0]),
         :𝑤=>(2,[1]),
         :𝝭=>(4,[1]),
     ])
     element = Element{:Poi1}((c,1,𝓒),(g,1,𝓖))
     push!(elements["Γᵗ"],element)

     return elements,nodes
end
    
