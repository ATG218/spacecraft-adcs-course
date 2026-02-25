#I wdot + w x (Iw) = tau

# I         : 3x3 inertia matrix about CoM
# w         : [w1, w2. w3]^T 
# H = Iw(t) : Iw angular mom body framerate
# w x (Iw)  : gyroscopic coupling term????

using LinearAlgebra
using StaticArrays
using DifferentialEquations
using GLMakie
using GeometryBasics


# Theres an error on 94 I marked, not sure how to fix it
#perhaps useful for question
rpm_to_rads(rpm) = rpm * 2π / 60
rads_to_rpm(w)   = w * 60 / (2π)

I_body = @SMatrix [
    0.020  0.000  0.000;
    0.000  0.030  0.000;            # kg*m^2
    0.000  0.000  0.050
]

E = eigen(Matrix(I_body))           # eigenvectors + eigenvalues
Q = SMatrix{3,3}(E.vectors)         # columns are the principal axes in body frame
Iprin = SVector{3}(E.values)          # principal moments [I1,I2,I3]
#Iprin = sort(Iprin)                #idk if we needa sort
I1, I2, I3 = Iprin

println("Principal moments (kg m^2): ", Iprin)
println("Minor = ", I1, "  Intermediate = ", I2, "  Major = ", I3)

function euler_free!(dw, w, p, t) #from class
    I1, I2, I3 = p
    w1, w2, w3 = w

    dw[1] = (I2 - I3)/I1 * w2*w3 
    dw[2] = (I3 - I1)/I2 * w3*w1
    dw[3] = (I1 - I2)/I3 * w1*w2
end

#q2 stability
wmag = rpm_to_rads(10.0)   # 10 RPM -> rad/s
H0 = I3 * wmag             # momentum magnitude for major-axis spin

err = rpm_to_rads(0.05)

w_major = @SVector [0.0, 0.0, H0/I3]  
w_inter = @SVector [0.0, H0/I2, 0.0]
w_minor = @SVector [H0/I1, 0.0, 0.0]

w_major_pert = w_major + @SVector [err, 0.0, 0.0]  
w_inter_pert = w_inter + @SVector [0.0, 0.0, err]
w_minor_pert = w_minor + @SVector [0.0, err, 0.0]

println("Major-axis RPM: ", rads_to_rpm(w_major[3]))
println("Inter-axis RPM: ", rads_to_rpm(w_inter[2]))
println("Minor-axis RPM: ", rads_to_rpm(w_minor[1]))


#sim it
p = (I1, I2, I3)
tspan = (0.0, 8000.0)

function simulate(w0; saveat=0.5)
    prob = ODEProblem(euler_free!, Vector(w0), tspan, p)
    solve(prob, Vern9(); reltol=1e-10, abstol=1e-12, saveat=saveat) #arb tolerances
end

sol_major = simulate(w_major_pert)
sol_inter = simulate(w_inter_pert)
sol_minor = simulate(w_minor_pert)

#momentum sphere!!
function H_traj(sol)
    H = Vector{SVector{3,Float64}}(undef, length(sol.u))
    for (k, wk) in enumerate(sol.u)
        w = @SVector [wk[1], wk[2], wk[3]]
        H[k] = @SVector [I1*w[1], I2*w[2], I3*w[3]]
    end
    H
end

H_major = H_traj(sol_major)
H_inter = H_traj(sol_inter)
H_minor = H_traj(sol_minor)

Hnorm = norm(H_major[1])
println("Hnorm:", Hnorm)


#theres some issue with this i dont know how to fix, everything else makes sense
eq = [
    @SVector [Hnorm, 0.0, 0.0], @SVector [-Hnorm, 0.0, 0.0],
    @SVector [0.0, Hnorm, 0.0], @SVector [0.0, -Hnorm, 0.0],
    @SVector [0.0, 0.0, Hnorm], @SVector [0.0, 0.0, -Hnorm]
]

#animation stuff i dont get
function to_arrays(H)
    hx = Float32.([h[1] for h in H])
    hy = Float32.([h[2] for h in H])
    hz = Float32.([h[3] for h in H])
    hx, hy, hz
end

hx1, hy1, hz1 = to_arrays(H_major)
hx2, hy2, hz2 = to_arrays(H_inter)
hx3, hy3, hz3 = to_arrays(H_minor)

fig = Figure(size=(1200, 800))
ax = Axis3(fig[1,1], aspect=:data)
hidedecorations!(ax)

# momentum sphere
sphere = Sphere(Point3f(0f0,0f0,0f0), Float32(Hnorm))
mesh!(ax, sphere; transparency=true, alpha=0.15)

# equilibrium points
eqx = Float32.([e[1] for e in eq]); eqy = Float32.([e[2] for e in eq]); eqz = Float32.([e[3] for e in eq])
scatter!(ax, eqx, eqy, eqz; markersize=12)

# trajectories
lines!(ax, hx1, hy1, hz1)  # major
lines!(ax, hx2, hy2, hz2)  # intermediate
lines!(ax, hx3, hy3, hz3)  # minor

# slider scrubs a moving marker on one trajectory (pick which)
sl = Slider(fig[2,1], range=1:length(hx1), startvalue=1)
pos = Observable(Point3f(hx1[1], hy1[1], hz1[1]))
scatter!(ax, pos; markersize=18)

on(sl.value) do k
    kk = Int(k)
    pos[] = Point3f(hx1[kk], hy1[kk], hz1[kk])
end

display(fig)