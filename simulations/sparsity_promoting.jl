using PRIMA
using CSV
using DataFrames
include("../src/basics.jl")
include("../src/sparse_basics.jl")


function mytarget(v_vals, n_v)
    return 1.0
end

function myconstraint_prealloc2(res, v_vals, n_v, moment_constraints, moment_generating_function)
    for i in 1:length(moment_constraints)
        @inbounds mc = moment_constraints[i]

        sum_mom = 0.0
        for j in 1:n_v
            @inbounds sum_mom += v_vals[j]^mc[1] * v_vals[n_v+j]^mc[2]
        end
        computed_mom = sum_mom / n_v
        @inbounds res[i] = computed_mom - moment_generating_function(mc)
    end
    return res
end


function check_next_moments(v_vals, n_v, moments_to_check, moment_generating_function)
    res = zeros(length(moments_to_check))
    computed = zeros(length(moments_to_check))
    analytical = zeros(length(moments_to_check))
    for i in 1:length(moments_to_check)
        mc = moments_to_check[i]

        # computed_mom = sum(v_vals[1:n_v].^mc[1])*sum(v_vals[n_v+1:end].^mc[2])/(n_v^2)
        computed_mom = sum((v_vals[1:n_v].^mc[1]) .* (v_vals[n_v+1:end].^mc[2]))/n_v
        res[i] = computed_mom - moment_generating_function(mc)
        computed[i] = computed_mom
        analytical[i] = moment_generating_function(mc)
        println("M($(mc[1]), $(mc[2]))_an = ", moment_generating_function(mc), " , comp: ", computed_mom, ", err = ", res[i])
    end

    return Dict("error"=>res, "computed"=>computed, "analytical"=>analytical)
end

function solve(n_v, λ1, n_subs, moms_conserved, nextmoments, vdf_func, metric)
    d_thresh1 = 1e-2
    d_thresh2 = 1e-3
    d_thresh3 = 1e-4

    lo_vel = -4.0
    hi_vel = 4.0
    dvx = abs(hi_vel - lo_vel) / n_subs
    dvy = dvx

    xlo = repeat([lo_vel], n_v*2)
    xhi = repeat([hi_vel], n_v*2)
    f_scale = 1.0 / n_v

    println("----------")
    println("λ1=$(λ1)")
    println("Conserving: ", moms_conserved)
    println("Predicting: ", nextmoments)

    res = zeros(length(moms_conserved))

    constraint_prima = (v_vals) -> myconstraint_prealloc2(res, v_vals, n_v, moms_conserved, vdf_func)

    sparsity_term_lambda = sparsity_term_2D
    if metric != "L1_comp"
        sparsity_term_lambda = sparsity_term_2D_g
    end

    target_prima = (v_vals) -> sparsity_term_lambda(v_vals, n_v)
    target_prima_me = (v_vals) -> sparsity_term_lambda(v_vals, n_v) + entropy_2D(v_vals, n_v, lo_vel, lo_vel, dvx, dvy, n_subs, n_subs, f_scale) / λ1

    target_prima_dummy = (v_vals) -> mytarget(v_vals, n_v)

    sol0 = vcat(Vector(LinRange(-1.0, 1.0, n_v)),Vector(LinRange(-1.0, 1.0, n_v)))
    println("Constraint error (sol0): ", maximum(abs.(constraint_prima(sol0))))
    E0 = entropy_2D(sol0, n_v, lo_vel, lo_vel, dvx, dvy, n_subs, n_subs, f_scale)
    println("S(sol0): ", E0)

    x_constraints, info = cobyla(target_prima_dummy, sol0; xl=xlo, xu=xhi, nonlinear_eq=constraint_prima)
    # x, info = cobyla(target_prima, x_constraints; xl=xlo, xu=xhi, nonlinear_eq=constraint_prima)
    x, info = cobyla(target_prima_me, x_constraints; xl=xlo, xu=xhi, nonlinear_eq=constraint_prima)

    println("N constraints: ", length(moms_conserved))
    constraint_error = maximum(abs.(constraint_prima(x)))
    println("Constraint error: ", constraint_error)

    A = zeros(Int16, n_v, n_v)
    nr1 = get_number_of_collapsed_points_graph(A, x, n_v, d_thresh1)
    nr2 = get_number_of_collapsed_points_graph(A, x, n_v, d_thresh2)
    nr3 = get_number_of_collapsed_points_graph(A, x, n_v, d_thresh3)
    println("Np reduced: ", nr1, ", ", nr2, ", ", nr3)
    errors = check_next_moments(x, n_v, nextmoments, vdf_func)
    Ex = entropy_2D(x, n_v, lo_vel, lo_vel, dvx, dvy, n_subs, n_subs, f_scale)
    println("S(x), L1(x): ", Ex, ", ", sparsity_term_2D(x, n_v))

    return Dict("E(constraints)"=>constraint_error,
                "Nreduced(1e-2)"=>nr1, "Nreduced(1e-3)"=>nr2, "Nreduced(1e-4)"=>nr3,
                "sol"=>x, "E0"=>E0, "Ex"=>Ex,
                "errors"=>errors)
end


function solve_constraints_only(n_v, moms_conserved, nextmoments, vdf_func, metric)

    # n_v = 30
    # λ1 = 1e9
    d_thresh1 = 1e-2
    d_thresh2 = 1e-3
    d_thresh3 = 1e-4

    lo_vel = -4.0
    hi_vel = 4.0
    # n_subs = 400

    xlo = repeat([lo_vel], n_v*2)
    xhi = repeat([hi_vel], n_v*2)
    f_scale = 1.0 / n_v

    # moms_conserved = generate_moments(1, 3)
    # nextmoments = generate_moments(4,4)
    println("----------")
    println("Conserving: ", moms_conserved)
    println("Predicting: ", nextmoments)

    res = zeros(length(moms_conserved))

    constraint_prima = (v_vals) -> myconstraint_prealloc2(res, v_vals, n_v, moms_conserved, vdf_func)

    sparsity_term_lambda = sparsity_term_2D
    if metric != "L1_comp"
        sparsity_term_lambda = sparsity_term_2D_g
    end

    target_prima_dummy = (v_vals) -> mytarget(v_vals, n_v)
    #3 (generic function with 1 method)

    sol0 = vcat(Vector(LinRange(-1.0, 1.0, n_v)),Vector(LinRange(-1.0, 1.0, n_v)))
    println("Constraint error (sol0): ", maximum(abs.(constraint_prima(sol0))))
    
    x, info = cobyla(target_prima_dummy, sol0; xl=xlo, xu=xhi, nonlinear_eq=constraint_prima)
    
    println("N constraints: ", length(moms_conserved))
    constraint_error = maximum(abs.(constraint_prima(x)))
    println("Constraint error: ", constraint_error)
    
    A = zeros(Int16, n_v, n_v)
    nr1 = get_number_of_collapsed_points_graph(A, x, n_v, d_thresh1)
    nr2 = get_number_of_collapsed_points_graph(A, x, n_v, d_thresh2)
    nr3 = get_number_of_collapsed_points_graph(A, x, n_v, d_thresh3)
    println("Np reduced: ", nr1, ", ", nr2, ", ", nr3)
    errors = check_next_moments(x, n_v, nextmoments, vdf_func)

    return Dict("E(constraints)"=>constraint_error,
                "Nreduced(1e-2)"=>nr1, "Nreduced(1e-3)"=>nr2, "Nreduced(1e-4)"=>nr3,
                "sol"=>x,
                "errors"=>errors)
end

function main(n_v_arr, λ1_arr, n_subs_arr, max_mom_c, max_next_mom, vdf_name, vdf_func; write_out=false, metric="L1_comp")

    path_to_output_prefix = "out/" * metric * "/" * vdf_name
    mkpath(path_to_output_prefix)

    moms_conserved = generate_moments(1, max_mom_c)
    nextmoments = generate_moments(max_mom_c+1,max_next_mom)

    for n_v in n_v_arr
        for n_subs in n_subs_arr
            results_arr = []
            for λ1 in λ1_arr
                res = solve(n_v, λ1, n_subs, moms_conserved, nextmoments, vdf_func, metric)
                push!(results_arr, res)
            end
            if write_out
                df_log = DataFrame(L1_lambda=λ1_arr,
                                   nreduced_1em2=[res["Nreduced(1e-2)"] for res in results_arr],
                                   nreduced_1em3=[res["Nreduced(1e-3)"] for res in results_arr],
                                   nreduced_1em4=[res["Nreduced(1e-4)"] for res in results_arr],
                                   error_constraints=[res["E(constraints)"] for res in results_arr],
                                   E0=[res["E0"] for res in results_arr],
                                   Ex=[res["Ex"] for res in results_arr]
                                   )
                CSV.write(path_to_output_prefix * "/" * "log_nv$(n_v)_nsub$(n_subs)_mc$(max_mom_c)", df_log)

                df_moments = DataFrame(L1_lambda=λ1_arr)
                for (i, nm) in enumerate(nextmoments)
                    df_moments[!, "error_$(nm[1])_$(nm[2])"] = [res["errors"]["error"][i] for res in results_arr]
                    df_moments[!, "computed_$(nm[1])_$(nm[2])"] = [res["errors"]["computed"][i] for res in results_arr]
                    df_moments[!, "analytical_$(nm[1])_$(nm[2])"] = [res["errors"]["analytical"][i] for res in results_arr]
                end
                CSV.write(path_to_output_prefix * "/" * "moments_prediction_nv$(n_v)_nsub$(n_subs)_mc$(max_mom_c)", df_moments)

                df_sol = DataFrame(L1_lambda=λ1_arr)
                for i in 1:n_v
                    df_sol[!, "vx_$i"] = [res["sol"][i] for res in results_arr]
                    df_sol[!, "vy_$i"] = [res["sol"][n_v+i] for res in results_arr]
                end
                CSV.write(path_to_output_prefix * "/" * "solution_nv$(n_v)_nsub$(n_subs)_mc$(max_mom_c)", df_sol)
            end
        end
    end
end


function constraint_only_main(n_v_arr, max_mom_c, max_next_mom, vdf_name, vdf_func; write_out=false, metric="L1_comp")

    path_to_output_prefix = "out/" * metric * "/" * vdf_name
    mkpath(path_to_output_prefix)

    moms_conserved = generate_moments(1, max_mom_c)
    nextmoments = generate_moments(max_mom_c+1,max_next_mom)

    for n_v in n_v_arr
        results_arr = []
        res = solve_constraints_only(n_v, moms_conserved, nextmoments, vdf_func, metric)
        push!(results_arr, res)
        if write_out
            df_log = DataFrame(nreduced_1em2=[res["Nreduced(1e-2)"] for res in results_arr],
                                nreduced_1em3=[res["Nreduced(1e-3)"] for res in results_arr],
                                nreduced_1em4=[res["Nreduced(1e-4)"] for res in results_arr],
                                error_constraints=[res["E(constraints)"] for res in results_arr]
                                )
            CSV.write(path_to_output_prefix * "/" * "constraints_only_log_nv$(n_v)_mc$(max_mom_c)", df_log)

            # -1.0 just a sanity check in the output later on, so we know that the solution is just a constraint-satisfying one
            df_moments = DataFrame(L1_lambda=[-1.0])
            for (i, nm) in enumerate(nextmoments)
                df_moments[!, "error_$(nm[1])_$(nm[2])"] = [res["errors"]["error"][i] for res in results_arr]
                df_moments[!, "computed_$(nm[1])_$(nm[2])"] = [res["errors"]["computed"][i] for res in results_arr]
                df_moments[!, "analytical_$(nm[1])_$(nm[2])"] = [res["errors"]["analytical"][i] for res in results_arr]
            end
            CSV.write(path_to_output_prefix * "/" * "constraints_only_moments_prediction_nv$(n_v)_mc$(max_mom_c)", df_moments)

            df_sol = DataFrame(L1_lambda=[-1.0])
            for i in 1:n_v
                df_sol[!, "vx_$i"] = [res["sol"][i] for res in results_arr]
                df_sol[!, "vy_$i"] = [res["sol"][n_v+i] for res in results_arr]
            end
            CSV.write(path_to_output_prefix * "/" * "constraints_only_solution_nv$(n_v)", df_sol)
        end
    end
end

const lambdas = [1e-3]
# const lambdas = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 2.5e2, 5e2, 7.5e2, 1e3, 5e3, 1e4, 5e4, 1e5, 1e6]
const entropy_subs = [10]

const mom_bimodal1 = (n) -> moment_Bimodal_2D(n, 0.25, 0.75, 1.2, 0.6, -0.4, -0.2, 0.5, 1.5)

for base_constraints in [2,3,4]
    println("BC=$(base_constraints)")
    
    # component-wise
    main([30], lambdas,
         entropy_subs, base_constraints, 6, "Bim1", mom_bimodal1, write_out=true, metric="L1_comp")
    # sqrt(g^2)
    main([30], lambdas,
         entropy_subs, base_constraints, 6, "Bim1", mom_bimodal1, write_out=true, metric="L1_g")
end