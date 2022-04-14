# Import the module. Using "using" keyword to make all export functions accessible.
using MORFEInvariantManifold
using MATLAB

# Name of the mesh file. The one of this example is a COMSOL mesh format.
mesh_file = "arch_ref.mphtxt"

# add material
MORFE_add_material("polysilicon",2.32e-3,160e3,0.22)

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
α = 0.8418240182000999/500.0
β = 0.0

# FORCING INFO
# list of excitation frequencies used to excite the system. 
# A parametrisation is performed for each frequency so 
# avoid using more thatn 1 or 2 frequencies
Ω_list = [1.0]
# For the moment we have modal-proportional forcing. 
# k_modes is a vector with length equal to the number of 
# excitation frequencies. For each frequency, a list of modes
# used to define the excitation forcing is reported.
# For instance, for the first excitation
# frequency Ω_list[i] the system is excited with a forcing vector
# defined by the combination of κ_modes[i] modes, wighted by the 
# load multipliers defined in κ_list[i], i.e.:
# F = sum_{j=1}^{n_frequencies} M κ_modes[i][j]*κ_list[i][j]
# the phase of the input force is given by κ_phase. 'c' stands for cosine,
# while 's' stands for sine.
κ_modes = [[1]]
κ_list = [[0.8418240182000999^2*0.06]]
κ_phase = [['c']]

# PARAMETRISATION INFO
# Φₗᵢₛₜ: vector of integer. Each integer correspond to which mode is included in the reduced model.
Φₗᵢₛₜ = [1]
# style: parametrisation style. 'g' graph, 'r' real normal form, 'c' complex normal form
style = 'c'
# max_order_a: maximum order of the asymptotic expansion of the autonomous problem
max_order_a = 9
# max_order_na: maximum order of the asymptotic expansion of the nonautonomous problem
max_order_na = 0

odir, Cp, Cp_na, rdyn = MORFE_mech_nonautonomous(mesh_file,domains_list,materials,
                                                 boundaries_list,constrained_dof,bc_vals,
                                                 Ω_list,κ_modes,κ_list,κ_phase,
                                                 α,β,
                                                 Φₗᵢₛₜ,style,max_order_a,max_order_na);


#
zero_amplitude = [0.0,0.0]
harmonics_init = [1.0,0.0]
time_integration_length = 15000.0
param_init = [0.0,1.0,0.8]
forward=true
MaxNumPoints=100.0
minstep=1e-8
maxstep=20.0
ncol=4.0
ntst=40.0
cont_param=3.0
analysis_number=1
MORFE_integrate_rdyn_frc(odir,zero_amplitude,harmonics_init,param_init,cont_param,
                                  time_integration_length,forward,MaxNumPoints,minstep,maxstep,ncol,ntst,analysis_number)
#
frf = MORFE_compute_frc_modal(odir,Ω_list)

# post proc
ω₀ = imag(Cp[2].f[1,1]) # extract eigenfrequency 
#
L = 6.4                 # characteristic length of the vibration amplitude
N=size(Cp[2].W)[1];
ϕ = maximum(abs.(Cp[2].W[Int(N/2)+1:N,1])) # maximum value of the mode
#

x = frf[:,1]./ω₀
y = frf[:,2].*(ϕ/L)
# plot results
eval_string("plot($x,$y)")
eval_string("xlim([1.000,1.05])")
eval_string("ylim([0.,0.5132])")