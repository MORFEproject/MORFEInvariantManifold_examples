# Import the module. Using "using" keyword to make all export functions accessible.
using MORFEInvariantManifold
using MATLAB

# Name of the mesh file. The one of this example is a COMSOL mesh format.
mesh_file = "arch_I.mphtxt"


### DOMAINS INFO
# domain_list is a vector that stores vectors of integers. 
# Each subvector is a domain. Each integer is a COMSOL volume.
domains_list = [                 
[1,2]
]

### MATERIAL 
# define material properties
material_name = "polysilicon";
density = 2.32e-3;
young_modulus = 160e3;
poisson_ratio = 0.22;
# add material to the library. 
# Once added is saved and can be used without redefining it
MORFE_add_material(material_name,density,young_modulus,poisson_ratio)
# assign materials to the domain:
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
style = 'r'
# max_order: maximum order of the asymptotic expansion of the autonomous problem
max_order = 5

@time odir, Cp, rdyn = MORFE_mech_autonomous(mesh_file,domains_list,materials,
                                       boundaries_list,constrained_dof,bc_vals,
                                       α,β,
                                       Φₗᵢₛₜ,style,max_order)
# post proc
ω₀ = imag(Cp[2].f[1,1]) # extract eigenfrequency
T₀ = 2*pi /ω₀  # extract period
#
L = 6.4   # characteristic length of the vibration amplitude
N = size(Cp[2].W)[1]; # number of physical DOFs 
ϕ = maximum(abs.(Cp[2].W[Int(N/2)+1:N,1])) # maximum value of the mode

small_amplitude = [0.0,0.3]
time_integration_length = 2*T₀
forward = true
MaxNumPoints = 90.0
minstep = 1e-8
maxstep = 10.0
initstep = 0.1
ncol = 4.0
ntst = 40.0

ms = MORFE_integrate_rdyn_backbone(odir,small_amplitude,
                                    time_integration_length,forward,MaxNumPoints,
                                    minstep,maxstep,initstep,ncol,ntst)
#
frf = MORFE_compute_backbone_modal(odir)


x = frf[:,1]./ω₀
y = frf[:,2].*(ϕ/L)
# plot results
put_variable(ms,:x,x)
put_variable(ms,:y,y)
put_variable(ms,:max_order,max_order)

show_msession(ms) # do not close the pop up matlab windows until done with the analyses
eval_string(ms,"
figure(1);hold on
plot(x,y,'DisplayName',strcat(\"Order \",num2str(max_order)))
xlim([0.98,1.04])
ylim([0.,0.7631])
xlabel('\$\\omega/\\omega_1\$','Interpreter','latex');
ylabel('max[\$u_1 \\phi_1\$]/\$L\$','Interpreter','latex');
legend()
")
close(ms)
