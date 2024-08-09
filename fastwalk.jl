using Random
using Dates
using DataStructures
using Laplacians
using Base.Threads
include("logw.jl")
include("graph.jl")

seed_value = Dates.millisecond(now())
Random.seed!(seed_value)

function loop_erased_walk(g::Graph, v, T)
    n = g.n
    intree = falses(n)
    next = fill(-1, n)
    root = fill(-1, n)
    re = zeros(n)
    for i=1:T
        intree .= false
        next .= -1
        root .= -1    
        intree[v] = true
        for j in 1:n
            u = j
            re[u] = re[u] + 1/T
            while !intree[u]
                next[u] = rand(g.nbr[u])
                u = next[u]
                re[u] = re[u] + 1/T
            end
            re[u] = re[u] - 1/T
            u = j
            while !intree[u] 
                intree[u] = true
                u = next[u]
            end
        end
    end
    return re
end


function loop_erased_walk_parra(g::Graph, v, T)
    n = g.n
    re = zeros(n)
    local_re = [zeros(n) for _ in 1:nthreads()]
    @threads for i=1:T
        t_id = threadid()
        intree = falses(n)
        next = fill(-1, n)
        root = fill(-1, n)    
        intree[v] = true
        for j in 1:n
            u = j
            local_re[t_id][u] = local_re[t_id][u] + 1/T
            while !intree[u]
                next[u] = rand(g.nbr[u])
                u = next[u]
                local_re[t_id][u] = local_re[t_id][u] + 1/T
            end
            local_re[t_id][u] = local_re[t_id][u] - 1/T
            u = j
            while !intree[u] 
                intree[u] = true
                u = next[u]
            end
        end
    end
    for t_id in 1:nthreads()
        re .+= local_re[t_id]
    end
    return re
end


function calcolumn(A, deg, deg_sqrt, v, n, m)
    sol = approxchol_lap(A; tol=1e-6)
    y = zeros(n)
    yy = zeros(n)
    yy[v] = deg_sqrt[v]
    for i=1:n
        y[i] = yy[i] - 1 / (2 * m) * deg[i] * deg_sqrt[v]
    end
    z = zeros(n)
    z .= sol(y)
    yyy = dot(deg, z)
    r = zeros(n)
    for i=1:n
        r[i] = deg_sqrt[i] * z[i] - 1 / (2 * m) * deg_sqrt[i] * yyy
    end
    return r
end

function fast(g, A, v, T, deg, deg_sqrt)
    n = g.n
    m = g.m
    re = loop_erased_walk_parra(g, v, T)
    col_v = calcolumn(A, deg, deg_sqrt, v, n, m)
    diags = zeros(n)
    for u=1:n
        diags[u] = re[u] - deg[u] / deg[v] * col_v[v] + 2 * deg_sqrt[u] /  deg_sqrt[v] * col_v[u]
    end
    return diags
end

str = "facebook"
if '#' in str
    continue
end
G = get_graph(str)
n = G.n
m = G.m
L = lapsp(G);
A = adjsp(G);
deg = zeros(n)
deg_sqrt = zeros(n)
for i=1:n
    deg[i] = L[i,i]
    deg_sqrt[i] = sqrt(deg[i])
end
