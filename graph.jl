struct Graph
    n :: Int # |V|
    m :: Int # |E|
    u :: Array{Int, 1}
    v :: Array{Int, 1} # uv is an edge
    w :: Array{Int, 1} # edge weight
    nbr :: Array{Array{Int, 1}, 1}
    name :: String
    degree :: Array{Int, 1}
end


include("core.jl")

function get_graph(ffname)
    n = 0
    Label = Dict{Int32, Int32}()
    Origin = Dict{Int32, Int32}()

    getID(x :: Int) = haskey(Label, x) ? Label[x] : Label[x] = n += 1
    fname = string("./graphs/",ffname,".txt")
    fin = open(fname, "r")

    str = readline(fin)
    u = Int[]
    v = Int[]
    str = strip(str)
    tot = 0
    i=0
    while str != ""
        if occursin("#", str)
            str = strip(readline(fin))
            continue
        end
        str = split(str)
        if length(str)!=2
            str = strip(readline(fin))
            continue
        end
        i = i + 1
        if length(str)<2
            println("error", i)
        end

        x   = parse(Int, str[1])
        y   = parse(Int, str[2])
        u1 = getID(x)
        v1 = getID(y)
        Origin[u1] = x
        Origin[v1] = y
        push!(u, u1)
        push!(v, v1)
        tot += 1
        str = strip(readline(fin))
    end
    nbr=[ [ ] for i in 1:n ]
    for i=1:tot
        u1=u[i];
        v1=v[i];
        push!(nbr[u1],v1);
        push!(nbr[v1],u1);
    end

    degree = zeros(n)
    for i=1:n
        degree[i] = length(nbr[i])
    end
    w = zeros(Float64, tot)
    for i=1:tot
        w[i] = 1
    end

    close(fin)
    return Graph(n, tot, u, v, w, nbr, ffname, degree)
end
