# Import the module. Using "using" keyword to make all export functions accessible.
using MORFEInvariantManifold

# Name of the mesh file. The one of this example is a COMSOL mesh format.
mesh_file = "arch_4.mphtxt"

# add material
MORFE_add_material("polysilicon",2.33e-3,160e3,0.22)

### DOMAINS INFO
# domain_list is a vector that stores vectors of integers. 
# Each subvector is a domain. Each integer is a COMSOL volume.
domains_list = [                 
[1,2]
]
# materials is an array of strings. 
# materials[i] embeds the material associated to the domain defined in domains_list[j]
materials = [
"polysilicon"
]

### BOUNDARIES INFO
# boundaries_list is a vector that stores vectors of integers.
# Each subvector is a boundary. Each integer is a COMSOL face
boundaries_list = [
[1,11]
]

# constrained degrees of freedom of the surfaces defined above.
# It is an array n_{surface}*dim, with n_{surface} number of boundaries 
# defined above and dim dimension of the problem (3D = 3). The triplet
# identifies which direction is constrained in the model.
constrained_dof = [
[1, 1, 1]
]

# It prescribes the Dirichlet boundary conditions of the mechanical
# problem. Currently it handles only homogeneous boundary conditions
# so leave it to zero.
bc_vals = [
[0.0, 0.0, 0.0]
]

# Rayleigh damping coefficients.
# The damping matrix is defined as C = α M + β K
α = 0.0
β = 0.0

# PARAMETRISATION INFO
# Φₗᵢₛₜ: vector of integer. Each integer correspond to which mode is included in the reduced model.
Φₗᵢₛₜ = [1]
# style: parametrisation style. 'g' graph, 'r' real normal form, 'c' complex normal form
style = 'c'
# max_order: maximum order of the asymptotic expansion of the autonomous problem
max_order = 9

odir, Cp, rdyn = MORFE_mech_autonomous(mesh_file,domains_list,materials,
                                       boundaries_list,constrained_dof,bc_vals,
                                       α,β,
                                       Φₗᵢₛₜ,style,max_order)

x0 = [0.0,0.3]
integration_length = 100.0
forward=true
MaxNumPoints=100
minstep=1e-8
maxstep=20.0
ncol=4.0
ntst=40.0
MORFE_integrate_rdyn_backbone(odir,x0,integration_length,forward,MaxNumPoints,minstep,maxstep,ncol,ntst)

frf = MORFE_compute_backbone_modal(odir)