"""
    cold_dark_matter(g; name = :c, kwargs...)

Create a particle species for cold dark matter in the spacetime with metric `g`.
"""
function cold_dark_matter(g; name = :c, kwargs...)
    description = "Cold dark matter"
    return matter(g; adiabatic = true, name, description, kwargs...)
end
