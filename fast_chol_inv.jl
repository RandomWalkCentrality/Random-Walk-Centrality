using Random
using Dates
using DataStructures
using SparseArrays
using LinearAlgebra
using Base.Threads
using ProgressMeter
using MAT

seed_value = Dates.millisecond(now())
Random.seed!(seed_value)

function fast_cholesky(graph_name::SubString{String}, ϵ::Float64, base_window_size::Int)
    file = matopen("ichol_factor_$(graph_name).mat")
    R = read(file, "L")
    close(file)
    
    R = SparseMatrixCSC(R)
    n = size(R, 1)
    
    # Compute the number of non-zero elements in each column as an estimate of the degree
     degrees = [nnz(R[:, i]) for i in 1:n]
     max_degree = maximum(degrees)
     
     tau_v = zeros(n)
     local_tau_v = [zeros(n) for _ in 1:nthreads()]
 
     chunk_size = max(1, div(n, nthreads()))
     
     Base.Threads.@threads for t in 1:nthreads()
         start_idx = (t-1) * chunk_size + 1
         end_idx = t == nthreads() ? n : t * chunk_size
         local_z_star = spzeros(n)
         local_Z_tilde = Dict{Int, SparseVector{Float64, Int64}}()
         
        for j in end_idx:-1:start_idx
            # Dynamically adjust window size
            adaptive_window_size = base_window_size * (1 + degrees[j] / max_degree)
             
            fill!(local_z_star, 0)
            compute_z_star!(local_z_star, R, local_Z_tilde, j)
             
            if nnz(local_z_star) <= 2 * log(n)
                 local_Z_tilde[j] = copy(local_z_star)
            else
                 k = find_truncation_k(local_z_star, ϵ)
                 local_Z_tilde[j] = spzeros(n)
                 trunc_k!(view(local_Z_tilde[j], :), local_z_star, k)
            end

            local_tau_v[t][j] = norm(local_Z_tilde[j])^2
             
            # Maintain adaptive sliding window
            while length(local_Z_tilde) > adaptive_window_size
                oldest_key = minimum(keys(local_Z_tilde))
                delete!(local_Z_tilde, oldest_key)
            end
        end
    end
 
     for t in 1:nthreads()
         tau_v .+= local_tau_v[t]
     end
     
     return tau_v
end

function compute_z_star!(z_star::SparseVector, R::SparseMatrixCSC, Z_tilde::Dict{Int, SparseVector{Float64, Int64}}, j::Int)
    fill!(z_star, 0)
    z_star[j] = 1 / R[j,j]
    
    for i in reverse(nzrange(R, j))
        row = R.rowval[i]
        if row > j && haskey(Z_tilde, row)
            val = R.nzval[i] / R[j,j]
            axpy!(-val, Z_tilde[row], z_star)
        end
    end
    dropzeros!(z_star)
end

function find_truncation_k(z_star::SparseVector, ϵ::Float64)
    total = sum(abs, z_star.nzval)
    if total < 1e-10
        @warn "z_star total is very small: $total"
        return nnz(z_star)
    end
    cumsum = 0.0
    for (i, val) in enumerate(sort!(abs.(z_star.nzval), rev=true))
        cumsum += val
        if cumsum / total > 1 - ϵ
            return i
        end
    end
    return nnz(z_star)
end

function trunc_k!(z_trunc::SubArray, z::SparseVector, k::Int)
    fill!(z_trunc, 0)
    if k >= nnz(z)
        copyto!(z_trunc, z)
    else
        pq = PriorityQueue{Int, Float64}(Base.Order.Reverse)
        for (idx, val) in zip(z.nzind, z.nzval)
            abs_val = abs(val)
            if length(pq) < k
                enqueue!(pq, idx => abs_val)
            elseif abs_val > peek(pq)[2]
                dequeue!(pq)
                enqueue!(pq, idx => abs_val)
            end
        end
        for (idx, _) in pq
            z_trunc[idx] = z[idx]
        end
    end
end

function run_fast_cholesky(graph_name, v, ϵ, strategy_name, num_thread)
    @info "running fast_cholesky..."
    diag_logf = open(string("./output_ic/", graph_name, "_fc_app_diags_extra_", ϵ, "_", strategy_name, "_v_", v, "_nth_", num_thread, ".txt"), "w")
    time_logf = open(string("./output_ic/", graph_name, "_time_extra_", strategy_name, "_v_", v, "_nth_", num_thread, ".txt"), "a")

    file = matopen("ichol_factor_$(graph_name).mat")
    R = read(file, "L")
    close(file)
    n = size(R, 1) + 1

    result, time, bytes, gctime, memallocs = @timed fast_cholesky(graph_name, ϵ, Int(20*round(log(n))))

    for e in result
        println(diag_logf, e)
    end
    println(time_logf, "$(ϵ) $time")

    close(diag_logf)
    close(time_logf)

    return result, time
end

function process_graph(graph_name)
    # Read v (vertex with maximum degree) from the MATLAB-generated file
    file = matopen("max_degree_vertex_$(graph_name).mat")
    v = Int(round(read(file, "v")))
    close(file)
    file = matopen("max_degree_vertex_$(graph_name).mat")
    v = Int(round(read(file, "v")))
    close(file)
    
    num_thread = 1
    time_logf = open(string("./output_ic/", graph_name, "_time_degree_v_extra_", v, "_nth_", num_thread, ".txt"), "w")
    close(time_logf)
    
    base_result = nothing
    for (i, ϵ) in enumerate([1e-3,1e-4,1e-5,1e-6,1e-7])
        app_diags, elapsed_time = run_fast_cholesky(graph_name, v, ϵ, "degree", num_thread)
        
        if i == 1
            base_result = app_diags
            @info "ϵ=$ϵ: time=$(round(elapsed_time, digits=3))s"
        else
            max_diff = maximum(abs.(app_diags .- base_result))
            @info "ϵ=$ϵ: max_diff=$(round(max_diff, digits=6)), time=$(round(elapsed_time, digits=3))s"
        end
    end
end

function main()
    fname = open("./filename_ic.txt", "r")
    str = readline(fname)
    nn = parse(Int, str)

    for _ in 1:nn
        graph_name = readline(fname)
        graph_name = strip(graph_name)
        '#' in graph_name && continue
        
        println("Processing graph: $graph_name")
        
        process_graph(graph_name)
        
        println() # Add a blank line between graphs
    end
end

main()