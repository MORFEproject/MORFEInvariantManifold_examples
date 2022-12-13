# Import the module. Using "using" keyword to make all export functions accessible.
using MORFEInvariantManifold
using MATLAB

# Name of the mesh file. The one of this example is a COMSOL mesh format.
# Meshes associated to the blade are used to benchmark the scalability of the code.
mesh_file = "plate_1.mphtxt"
# blade_1 is the smallest mesh, blade_9 is the biggest.
# unzip blade_9 before using it

### DOMAINS INFO
# domain_list is a vector that stores vectors of integers. 
# Each subvector is a domain. Each integer is a COMSOL volume.
domains_list = [                 
[1]
]

### MATERIAL 
# define material properties
material_name = "Titanium";
density = 4400.0;
young_modulus = 104e9;
poisson_ratio = 0.3;
# add material to the library. 
# Once added is saved and can be used without redefining it
MORFE_add_material(material_name,density,young_modulus,poisson_ratio)
# assign materials to the domain:
# materials is an array of strings. 
# materials[i] embeds the material associated to the domain defined in domains_list[j]
materials = [
"Titanium"
]

### BOUNDARIES INFO
# boundaries_list is a vector that stores vectors of integers.
# Each subvector is a boundary. Each integer is a COMSOL face
boundaries_list = [[1]]

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
α = 0.3
β = 0.0

# Nonlinear options
nls = 1
matROT = [0.0, 0.0, 100.0] * 2 * pi / 60

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
# F = sum_{j = 1}^{n_frequencies} M κ_modes[i][j]*κ_list[i][j]
# the phase of the input force is given by κ_phase. 'c' stands for cosine,
# while 's' stands for sine.
κ_modes = [[1]]
κ_list = [[1.0]]
κ_phase = [['c']]


# PARAMETRISATION INFO
# Φₗᵢₛₜ: vector of integer. Each integer correspond to which mode is included in the reduced model.
Φₗᵢₛₜ = [1]
neig = maximum(Φₗᵢₛₜ)
# style: parametrisation style. 'g' graph, 'r' real normal form, 'c' complex normal form
style = 'r'
# max_order_a: maximum order of the asymptotic expansion of the autonomous problem
max_order_a = 7
# max_order_na: maximum order of the asymptotic expansion of the nonautonomous problem
max_order_na = 0


@time odir, Cp, Cp_na, rdyn = MORFE_mech_nonautonomous(mesh_file,domains_list,materials,
                                       boundaries_list,constrained_dof,bc_vals,
                                       Ω_list,κ_modes,κ_list,κ_phase,
                                       α,β,
                                       Φₗᵢₛₜ,style,max_order_a,
                                       max_order_na,neig,nls,matROT
                                       );

# post proc
ω₀ = imag(Cp[2].f[1,1]) # extract eigenfrequency
T₀ = 2*pi /ω₀  # extract period
#
L = 1.02  # characteristic length of the vibration amplitude
N = size(Cp[2].W)[1]; # number of physical DOFs 
ϕ = maximum(abs.(Cp[2].W[Int(N/2)+1:N,1])) # maximum value of the mode
#
#
# param_init:  [artificial damping parameter = 0,  κ multiplier = 1, initial Ω = 0.8*ω₀]
param_init = [0.0,1.0,0.95*ω₀] 
cont_param = 3.0 # continuation parameter is Ω
#
zero_amplitude = [0.0,0.0]
harmonics_init = [1.0,0.0]
time_integration_length = 500*T₀
forward = true
MaxNumPoints = 110.0
minstep = 1e-8
maxstep = 0.15
initstep = 0.1
ncol = 4.0
ntst = 40.0

ms = MORFE_integrate_rdyn_frc(odir,zero_amplitude,harmonics_init,param_init,cont_param,
                                       time_integration_length,forward,MaxNumPoints,
                                       minstep,maxstep,initstep,ncol,ntst)
#
frf = MORFE_compute_frc_modal(odir,Ω_list)



x = frf[:,1]./ω₀
y = frf[:,2].*(ϕ/L)
# plot results
put_variable(ms,:x,x)
put_variable(ms,:y,y)
put_variable(ms,:max_order_a,max_order_a)

show_msession(ms) # do not close the pop up matlab windows until done with the analyses
eval_string(ms,"
figure(1);hold on
plot(x,y,'DisplayName',strcat(\"Order \",num2str(max_order_a)))
xlabel('\$\\Omega/\\omega_1\$','Interpreter','latex');
ylabel('max[\$u_1 \\phi_1\$]/\$L\$','Interpreter','latex');
xlim([0.96,1.02])
ylim([0.,0.14])
legend();
")
close(ms)
