export Recombination
export Xe_Saha_H, Xe_Saha_H_He, Xe_Peebles, Xe
export τ, τ′, τ′′
export g, g′, g′′
export s
export x_switch_Peebles, x_decoupling, x_recombination, x_reionization_H, x_reionization_He

struct Recombination
    bg::Background
    xswitch::Float64
    Xe_Peebles::ODESolution # free electron fraction
    τ::ODESolution # optical depth
    s::ODESolution # sound horizon

    function Recombination(bg::Background; x0=-20.0, xswitch=NaN)
        if isnan(xswitch)
            xswitch = x_switch_Peebles(bg.par)
        end

        Xe_Peebles = integrate_Xe_Peebles(bg.par, xswitch)
        τ = integrate_τ(bg.par, Xe_Peebles, xswitch; x0=x0)
        s = integrate_s(bg.par; x0=x0)

        new(bg, xswitch, Xe_Peebles, τ, s)
    end
end

# for convenience
ρcrit(par::Parameters, x) = 3*H(par,x)^2 / (8*π*G)
ρb(par::Parameters, x)    = Ωb(par, x) * ρcrit(par, x)
nb(par::Parameters, x)    = ρb(par, x) / mH
nH(par::Parameters, x)    = (1-par.Yp) * nb(par, x)
nHe(par::Parameters, x)   = par.Yp/4 * nb(par, x)
Tγ(par::Parameters, x)    = par.Tγ0 / a(x)
λdB(m, T) = h / √(2*π*m*kB*T) # de Broglie wavelength for mass m and temperature T

function Xe_Saha_H(par::Parameters, x)
    Tb = Tγ(par, x)
    a = 1
    b = 1 / λdB(me,Tb)^3 * exp(-EH1ion/(kB*Tb)) / nH(par,x)
    c = -b

    # when b >> 1, the quadratic equation solution is
    #   (-b + √(b^2+4*b)) / 2
    # = b/2 * (-1 + √(1+4/b))
    # ≈ b/2 * (-1 + 1 + 2/b - 2/b^2)  (b >> 1)
    # = 1 - 1/b
    return b < 1e10 ? quadroots(a, b, c)[2] : 1 - 1/b # choose Taylor expansion for large b
end

function Xe_Saha_H_He(par::Parameters, x; tol=1e-15)
    # Fixed-point iterate for Xe, starting with the fast H-only Saha equation as the initial guess
    # If Xe0 is small, the iterations converges very slowly, but then there is almost no He anyway, so we just take the H-only solution
    Xe0 = Xe_Saha_H(par, x)

    # To find Xe = XH+ + Yp/(4*(1-Yp)) * (XHe+ + 2*XHe++),
    # begin with an initial guess for Xe and iteratively solve the system of Saha equations
    # (1) ne * XHe+ / (1 - XHe+ - XHe++) = 2 / λe^3 * exp(-EHe1ion / (kB*Tb))
    # (2) ne * XHe++ / XHe+ = 4 / λe^3 * exp(-EHe2ion/(kB*Tb))
    # (3) ne * XH+ / (1-XH+) = 1 / λe^3 * exp(-EH1ion/(kB*Tb))
    # for {XH+, XHe+, XHe++} = {XH1, XHe1, XHe2}.
    function Xe_Saha_H_He_fixed_point(Xe)
        Tb = Tγ(par, x)
        ne = Xe * nH(par, x)
        R1 = 1 / λdB(me,Tb)^3 * exp(-EH1ion /(kB*Tb)) / ne
        R2 = 2 / λdB(me,Tb)^3 * exp(-EHe1ion/(kB*Tb)) / ne
        R3 = 4 / λdB(me,Tb)^3 * exp(-EHe2ion/(kB*Tb)) / ne
        XH1  = 1 / (1 + 1/R1)
        XHe1 = 1 / (1 + 1/R2 + R3)
        XHe2 = 1 / (1 + 1/R3 + 1/(R2*R3))
        return XH1 + par.fHe * (XHe1 + 2*XHe2)
    end

    return Xe0 < 0.5 ? Xe0 : fixed_point_iterate(Xe_Saha_H_He_fixed_point, Xe0; tol=tol)
end

function x_switch_Peebles(par::Parameters)
    return find_zero(x -> Xe_Saha_H_He(par, x) - 0.999, (-20.0, +20.0), rtol=1e-20, atol=1e-20)
end

function Xe_reionization(par::Parameters, x)
    if !par.reionization
        return 0.0
    end

    y(z) = (1+z)^(3/2)
    dy_dz(z) = 3/2 * (1+z)^(1/2)
    Δy(z, Δz) = dy_dz(z) * Δz
    smoothstep(y, Δy, h) = h/2 * (1 + tanh(y / Δy))
    Xe_reionization_total  = smoothstep(y(par.z_reion_H ) - y(z(x)), Δy(par.z_reion_H,  par.Δz_reion_H),  1+par.fHe)
    Xe_reionization_total += smoothstep(y(par.z_reion_He) - y(z(x)), Δy(par.z_reion_He, par.Δz_reion_He), 0+par.fHe)
    return Xe_reionization_total
end

function integrate_Xe_Peebles(par::Parameters, xswitch)
    return integrate_Xe_Peebles(par, xswitch, Xe_Saha_H_He(par, xswitch))
end

function integrate_Xe_Peebles(par::Parameters, x1, Xe1)
    function dXe_dx(x, Xe)
        Tb = Tγ(par,x) # K
        n_1s = (1-Xe) * nH(par,x) # 1/m^3
        Λ_2s_1s = 8.227 # 1/s
        Λ_α = H(par,x) * (3*EH1ion/(ħ*c))^3 / ((8*π)^2 * n_1s) # 1/s
        ϕ2 = 0.448 * log(EH1ion/(kB*Tb)) # dimensionless
        α2 = 64*π / √(27*π) * (α/me)^2 * √(EH1ion/(kB*Tb)) * ϕ2 * ħ^2/c # m^3/s
        β  = α2 / λdB(me,Tb)^3 * exp(-EH1ion/(kB*Tb))
        β2 = α2 / λdB(me,Tb)^3 * exp(-EH1ion/(4*kB*Tb)) # 1/s (compute this instead of β2 = β * exp(3*EH1ion/(4*kB*Tb)) to avoid exp overflow)
        C_r = (Λ_2s_1s + Λ_α) / (Λ_2s_1s + Λ_α + β2)
        return C_r / H(par,x) * (β*(1-Xe) - nH(par,x)*α2*Xe^2)
    end
    return solve(ODEProblem{false}((Xe,_,x) -> dXe_dx(x, Xe), Xe1, (x1, 0.0)), Tsit5(); abstol=1e-8, reltol=1e-8)
end

Xe(par::Parameters, Xe_Peebles::ODESolution, xswitch, x) = (x < xswitch ? Xe_Saha_H_He(par, x) : Xe_Peebles(x)) + Xe_reionization(par, x)
Xe(rec::Recombination, x) = Xe(rec.bg.par, rec.Xe_Peebles, rec.xswitch, x)

function integrate_τ(par::Parameters, Xe_Peebles::ODESolution, xswitch; x0=-20.0)
    ne(x) = nH(par,x) * Xe(par, Xe_Peebles, xswitch, x)
    dτ_dx(x) = -ne(x) * σT * c / H(par,x)
    return solve(ODEProblem((τ,_,x) -> dτ_dx(x), 0.0, (0.0, x0)), Tsit5(); abstol=1e-10, reltol=1e-10)
end

τ(rec::Recombination, x) = rec.τ(x)
τ′(rec::Recombination, x) = ForwardDiff.derivative(x ->  τ(rec, x), x)
τ′′(rec::Recombination, x) = ForwardDiff.derivative(x -> τ′(rec, x), x)

g(rec::Recombination, x) = ForwardDiff.derivative(x -> exp(-τ(rec, x)), x)
g′(rec::Recombination, x) = ForwardDiff.derivative(x ->  g(rec, x), x)
g′′(rec::Recombination, x) = ForwardDiff.derivative(x -> g′(rec, x), x)

function integrate_s(par::Parameters; x0=-20.0)
    R(x) = 4*par.Ωγ0 / (3*par.Ωb0*a(x))
    cs(x) = c * √(R(x) / (3*(1+R(x))))
    ds_dx(x, s) = cs(x) / aH(par, x)
    s0 = cs(x0) / aH(par, x0)
    return solve(ODEProblem((s,_,x) -> ds_dx(x,s), s0, (x0, 0.0)), Tsit5(); abstol=1e-10, reltol=1e-10)
end

s(rec::Recombination, x) = rec.s(x)

x_reionization_H(par::Parameters)  = par.reionization ? x(1/(1+par.z_reion_H))  : NaN
x_reionization_He(par::Parameters) = par.reionization ? x(1/(1+par.z_reion_He)) : NaN
x_decoupling(rec::Recombination) = find_zero(x -> τ′(rec,x)^2 - τ′′(rec,x), (-20.0, -3.0)) # equivalent to g′=0 without the exponential; exclude reionization peak for x > -3 # TODO: -3.0 -> dynamic limit from reionization redshift
x_recombination(rec::Recombination) = find_zero(x -> Xe(rec,x) - 0.1, (-20.0, -3.0)) # exclude reionization for x > -3 # TODO: -3.0 -> dynamic limit from reionization redshift
