using SpecialFunctions
using HCubature

function generate_moments(m_start, m_end)
    res = Vector{Vector{Int}}()
    for m in m_start:m_end
        for i in 0:m
            push!(res, [i, m-i])
        end
    end
    return res
end

function moment_Maxwell_2D(n)
    if (n[1] % 2 == 0) && (n[2] % 2 == 0)
        return gamma(0.5 * (n[1] + 1)) * gamma(0.5 * (n[2] + 1)) / π
    else
        return 0.0
    end
end

function moment_Druyvesteyn(n)
    if (n % 2 == 0)
        return A_Druyvesteyn * gamma((n+1)/4.0) / (2 * alpha_Druyvesteyn^((n+1)/4.0))
    else
        return 0.0
    end
end

function integrand_Druyvesteyn_2D(n, vx, vy)
    # v_x^n[1] * v_y^n[2] * exp(-(v_x^2+v_y^2)^2)
    return (vx^n[1]) * (vy^n[2]) * exp(-((vx^2 + vy^2)^2)/π)
end

function integrand_Druyvesteyn_2D_transform(n, t)
    # t ∈ [0,1]
    return integrand_Druyvesteyn_2D(n, t[1] / (1.0-t[1]^2), t[2] / (1.0-t[2]^2)) * (1.0+t[1]^2)/((1.0-t[1]^2)^2) * (1.0+t[2]^2)/((1.0-t[2]^2)^2)
end

function moment_Druyvesteyn_2D(n)
    if (n[1] % 2 == 0) && (n[2] % 2 == 0)
        vdf_d = (t) -> integrand_Druyvesteyn_2D_transform(n, t)
        hc = hcubature(vdf_d, [-1.0, -1.0], [1.0, 1.0])
        # println(hc[2])
        return (2.0/(π^2)) * hc[1]
    else
        return 0.0
    end
end

function integrand_Bimodal_2D(n, vx, vy, rho1, rho2, vx0_1, vy0_1, vx0_2, vy0_2, T1, T2)
    # v_x^n[1] * v_y^n[2] * exp(-(v_x^2+v_y^2)^2)
    return (vx^n[1]) * (vy^n[2]) * ((rho1/T1) * exp(-((vx-vx0_1)^2 + (vy-vy0_1)^2)/T1) + (rho2/T2) * exp(-((vx-vx0_2)^2 + (vy-vy0_2)^2)/T2))
end

function integrand_Bimodal_2D_transform(n, t, rho1, rho2, vx0_1, vy0_1, vx0_2, vy0_2, T1, T2)
    # t ∈ [0,1]
    return integrand_Bimodal_2D(n, t[1] / (1.0-t[1]^2), t[2] / (1.0-t[2]^2), rho1, rho2, vx0_1, vy0_1, vx0_2, vy0_2, T1, T2) * (1.0+t[1]^2)/((1.0-t[1]^2)^2) * (1.0+t[2]^2)/((1.0-t[2]^2)^2)
end


function moment_Bimodal_2D(n, rho1, rho2, vx0_1, vy0_1, vx0_2, vy0_2, T1, T2)
    if ((n[1] == 1) && (n[2] == 0)) || ((n[1] == 0) && (n[2] == 1)) # Vx=0, Vy=0 due to choice of bimodal parameters, but quadrature takes forever to compute
        return 0.0
    else 
    # if (n[1] != 1) || (n[2] != 1)
        vdf_d = (t) -> integrand_Bimodal_2D_transform(n, t, rho1, rho2, vx0_1, vy0_1, vx0_2, vy0_2, T1, T2)
        hc = hcubature(vdf_d, [-1.0, -1.0], [1.0, 1.0])
        # println(hc[2])
        return (1.0/(π)) * hc[1]
    # else
    #     return 0.0
    end
end

function moment_triple_dirac_2D(n, vx1, vx2, vx3, vy1, vy2, vy3, w1, w2, w3)
    return w1 * (vx1^n[1]) * (vy1^n[2]) + w2 * (vx2^n[1]) * (vy2^n[2]) + w3 * (vx3^n[1]) * (vy3^n[2])
end

function moment_triple_dirac_2D_example(n)
    # 
    # return moment_triple_dirac_2D(n, -8.0, 2.0, 2.0, -2.0, 0.0, 4.0, 0.2, 0.7, 0.1)
    return moment_triple_dirac_2D(n, -1.0, 0.25, 0.25, -0.25, 0.0, 1.0, 0.2, 0.75, 0.05)
end


function moment_symmetric_double_dirac_2D_example(n)
    # 
    return moment_triple_dirac_2D(n, -2.0, 1.0, 0.0, -1.0, 0.5, 0.0, 1.0/3.0, 2.0/3.0, 0.0)
end
