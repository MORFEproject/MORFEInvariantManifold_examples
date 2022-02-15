using MORFEInvariantManifold

mesh_file = "cantilever_hexa.mphtxt"

### DOMAINS INFO
# list of the domain in arrays. Each entry contains the domain tags
domains_list = [                 
[1]
]
# materials : array of strings. One material for each domain
materials = [
"Titanium"
]

### BOUNDARIES INFO
# list of boundaries. Each entry contains the surfaces associated to a boundary.
boundaries_list = [
[1]
]

# constrained degrees of freedom of the differents surfaces defined above
constrained_dof = [
[1, 1, 1]
]

# values of the Dirichlet boundary condition
bc_vals = [
[0.0, 0.0, 0.0]
]

# damping
α = 0.0
β = 0.0

# parametrisation infos
Φₗᵢₛₜ = [1]
style = 'c'
max_order = 5


Cp, rdyn = MORFE_mech_autonomous(mesh_file,domains_list,materials,
                                 boundaries_list,constrained_dof,bc_vals,
                                 α,β,
                                 Φₗᵢₛₜ,style,max_order);