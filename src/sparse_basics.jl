using Graphs

function sparsity_term_2D(v_vals, n_v)
    res = 0.0

    # cut down computations by a factor of 2
    
    for i1 in 1:n_v-1
        for i2 in i1+1:n_v
            @inbounds res += abs(v_vals[i1] - v_vals[i2]) # iterate over x velocities
            @inbounds res += abs(v_vals[i1+n_v] - v_vals[i2+n_v]) # iterate over y velocities
        end
    end

    
    # for i1 in 1:n_v-1
    #     for i2 in i1+1:n_v
    #         res += abs(v_vals[i1+n_v] - v_vals[i2+n_v])
    #     end
    # end

    return 2*res/(n_v^2)
end

function sparsity_term_2D_g(v_vals, n_v)
    res = 0.0
    
    for i1 in 1:n_v-1
        for i2 in i1+1:n_v
            @inbounds res += sqrt((v_vals[i1] - v_vals[i2])^2+(v_vals[i1+n_v] - v_vals[i2+n_v])^2)
        end
    end

    return 2*res/(n_v^2)
end


function get_bin(v, vmin, dv)
    return trunc(Int16, (v-vmin) / dv) + 1
end

function entropy_2D(v_vals, n_v, vx_min, vy_min, dvx, dvy, nvxint, nvyint, f_scaling)
    res = Dict{Int64,Float64}()

    for i in 1:n_v
        binx = get_bin(v_vals[i], vx_min, dvx)
        biny = get_bin(v_vals[n_v+i], vy_min, dvy)

        bin_n = biny*nvxint + binx

        if !haskey(res, bin_n)
            res[bin_n] = 1
        else
            res[bin_n] += 1
        end
    end

    S = 0.0
    unity = 0.0
    for (key, value) in res
        S += dvx * dvy * f_scaling * value * log(f_scaling * value)
        unity += dvx * dvy * f_scaling * value
    end

    return S / unity - log(unity)
end

function compute_ndens(v_vals, f_scaling)
    return length(v_vals) * f_scaling
end

function compute_vx(v_vals, f_scaling)
    return sum(v_vals)
end

function compute_vx_2D(v_vals, n_v, f_scaling)
    return sum(v_vals[1:n_v])
end

function compute_vy_2D(v_vals, n_v, f_scaling)
    return sum(v_vals[n_v+1:end])
end

function compute_moment(M, ndens, v_vals, f_scaling)
    vx = compute_vx(v_vals, f_scaling)
    return sum((v_vals .- vx).^M) * f_scaling / ndens
end

function compute_moment_2D(M, ndens, v_vals, n_v, f_scaling)
    vx = compute_vx_2D(v_vals, n_v, f_scaling)
    vy = compute_vy_2D(v_vals, n_v, f_scaling)
    return sum(((v_vals .- vx).^M[1]) .* ((v_vals .- vy).^M[2])) * f_scaling / ndens
end

function get_init_positions_top_moment(m_vec, m_orders)
    pos = abs(m_vec[end])^(1.0/m_orders[end])
    return [-pos, pos]
end

function pretty_print_moments(m_vec, m_orders)
    for (mo, mvv) in zip(m_orders, m_vec)
        println("M$mo = $mvv")
    end
end

function pretty_print_moments_2D(m_vec, m_orders)
    for (mo, mvv) in zip(m_orders, m_vec)
        println("M($(mo[1]),$(mo[2])) = $mvv")
    end
end

function pretty_print_moments_with_reference(m_vec, m_vec_ref, m_orders)
    for (mo, mvv, mvvr) in zip(m_orders, m_vec, m_vec_ref)
        err = abs(mvvr - mvv)
        println("M$mo = $mvv, MRef$mo = $mvvr, err=$err")
    end
end

function pretty_print_moments_with_reference_2D(m_vec, m_vec_ref, m_orders)
    for (mo, mvv, mvvr) in zip(m_orders, m_vec, m_vec_ref)
        err = abs(mvvr - mvv)
        println("M($(mo[1]),$(mo[2])) = $mvv, MRef$mo = $mvvr, err=$err")
    end
end

function get_distances(v_vals)
    distances = abs.(v_vals[2:end] - v_vals[1:end-1])
    return minimum(distances), sum(distances) / length(distances), maximum(distances)
end

function get_number_of_collapsed_points(v_vals, thresh=1e-5)
    distances = abs.(v_vals[2:end] - v_vals[1:end-1])
    return sum(distances .< thresh)
end


function get_number_of_collapsed_points_2D(v_vals, n_v, thresh=1e-5)
    n_close = 0
    thr2 = thresh^2
    for i in 1:n_v-1
        for j in i+1:n_v
            dist2 = (v_vals[i]-v_vals[j])^2 + (v_vals[n_v+i]-v_vals[n_v+j])^2
            if dist2 < thr2
                n_close += 1
            end
        end
    end
    return n_close
end

function create_adjacency_matrix(A, v_vals, n_v, thresh=1e-5)
    A[:,:] .= 0
    thr2 = thresh^2
    for i in 1:n_v-1
        for j in i+1:n_v
            dist2 = (v_vals[i]-v_vals[j])^2 + (v_vals[n_v+i]-v_vals[n_v+j])^2
            if dist2 < thr2
                A[i,j] = 1
                A[j,i] = 1
            end
        end
    end
end

function get_number_of_collapsed_points_graph(A, v_vals, n_v, thresh=1e-5)
    create_adjacency_matrix(A, v_vals, n_v, thresh)
    G = SimpleGraph(A)
    return length(connected_components(G))
end

function compute_moment_vx(M, ndens, v_vals, f_scaling, vx)
    return sum((v_vals .- vx).^M) * f_scaling / ndens
end

function moment_constraint(v_vals,
                           which_moments_to_conserve, reference_moments, ndens, f_scaling)
    result = zeros(length(which_moments_to_conserve))
    vx = compute_vx(v_vals, f_scaling)
    for (i, M) in enumerate(which_moments_to_conserve)
        result[i] = compute_moment_vx(M, ndens, v_vals, f_scaling, vx) - reference_moments[i]
    end
    return result
end