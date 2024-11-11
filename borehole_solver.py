#
# Imports
#

import dolfin
import gmsh
import meshio
from fenics import *
import matplotlib.pyplot as plt
import time

#
# Define constants
#

#Ground: Granite
rho_ground = 2915 # density:2915  kg/m**3     #<-------------------------change
c_p_ground = 2240*10**3 # heat capacity: 861 J/m*K
k_ground = 3.15 # Thermal conductivity: 3.27 W/m*k
ground_t = 8.4 # Average ground temperature (BC)

#inner_pipe: MDPE
rho_inner_pipe = 938.5 # density: 938.5 Kg/m**3   #<-------------------------change
c_p_inner_pipe = 2300 # heat capacity: 2300 J/kg*K
k_inner_pipe = 0.40 # Thermal conductivity: 0.14 W/m*k

#Fluid (outer pipe): Water
rho_fluid = 999 # density: 999 kg/m^3      #<-------------------------change
c_p_fluid = 4.19*10**(6) # heat capacity: 4.19*10**(6) J/kg*K            #<-------------------------change
k_fluid = 0.590 # Thermal conductivity: 0.59 W/m*K
mu_fluid = 1.138*10**(-3) # dynmic viscosity: 1.138*10**(-3) Pa*s

fluid_t = 16. # Water temperature inlet (BC)
inflow_profile = ('0.', '-0.000596/(0.02*pi)')    # Volumetric flow 0.596 L/s----# m/s    #<-------------------------change
fluid_pressure = 1.0133e5 # Atmospheric pressure: 1.0133e5 Pa


#
# Define External force
#

alpha = 207*10**(-6) # Thermal expasion coeffcient
g = 9.81 # Graviational acceleration: 9.81 m/s^2
heat_source = Constant ( 0.0 ) # The heat input rate: 0 W
f_f  = Constant((0, 0)) # External force field applied to the fluid

#
# Define Load mesh
#

import meshio
from fenics import *
import matplotlib.pyplot as plt

# Step 1: Load mesh #####################################################################################

# READS MSH FILE AND CREATE XDMF FILES
mesh_from_file = meshio.read("Coaxial_Borehole_Heat_Exchanger.msh")
def create_mesh(mesh, cell_type, prune_z=False):
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:, :2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(
        points=points, cells={cell_type: cells}, cell_data={"name_to_read": [cell_data]}
    )
    return out_mesh

line_mesh = create_mesh(mesh_from_file, "line", prune_z=True)
meshio.write("facet_mesh.xdmf", line_mesh)

triangle_mesh = create_mesh(mesh_from_file, "triangle", prune_z=True)
meshio.write("mesh.xdmf", triangle_mesh)

# READS FILES AND CREATE MESH
line_mesh = create_mesh(mesh_from_file, "line", prune_z=True)
meshio.write("facet_mesh.xdmf", line_mesh)

triangle_mesh = create_mesh(mesh_from_file, "triangle", prune_z=True)
meshio.write("mesh.xdmf", triangle_mesh)

mesh = Mesh()
with XDMFFile("mesh.xdmf") as infile:
    infile.read(mesh)
#
mvc_2d = MeshValueCollection("size_t", mesh, 2)

with XDMFFile("mesh.xdmf") as infile:
    infile.read(mvc_2d, "name_to_read")
mf2 = cpp.mesh.MeshFunctionSizet(mesh, mvc_2d)
#
mvc_1d = MeshValueCollection("size_t", mesh, 1)


with XDMFFile("facet_mesh.xdmf") as infile:
    infile.read(mvc_1d, "name_to_read")
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc_1d)

print('num_cells   :', mesh.num_cells())
print('num_vertices:', mesh.num_vertices())
print('cell_type   :', mesh.ufl_cell())

subdomain_markers = mf2
boundary_markers = mf

## Step 1.5: Create submeshes and ds/dx ###############################################################

# Extract SubMesh for physical group tag 1,2,3
ground = SubMesh(mesh, subdomain_markers, 1)
fluid = SubMesh(mesh, subdomain_markers, 2)
int_pipe = SubMesh(mesh, subdomain_markers, 3)

# All mesh
dx=Measure("dx", domain=mesh, subdomain_data=subdomain_markers)
ds=Measure("ds", domain=mesh, subdomain_data=boundary_markers)

# air
boundary_fluid = MeshFunction("size_t",fluid,fluid.topology().dim()-1,0)
dx_fluid=Measure("dx", domain=fluid, subdomain_data=boundary_fluid)
ds_fluid=Measure("ds", domain=fluid, subdomain_data=boundary_fluid)

#
# Define boundaries
#

# FLUID ##########################################################################3

V = VectorElement("CG", triangle, 2)
Q = FiniteElement("CG", triangle, 1)
#
W = FunctionSpace(fluid, V * Q)

# Define boundaries ------------------------
#
class inflow_fluid(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], ground_len) and x[0]<=outer_pipe_radius and on_boundary
inflow = inflow_fluid()

class outflow_fluid(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], ground_len) and x[0] > outer_pipe_radius and on_boundary
outflow  = outflow_fluid()

class walls_boundary(SubDomain):
    def inside(self, x, on_boundary):
      top0 = 0
      top1 = outer_pipe_radius - inner_pipe_radius
      top2 = outer_pipe_radius + inner_pipe_radius
      top3 = 2*outer_pipe_radius
      #corners = (near(x[0], top0) or near(x[0], top1) or near(x[0], top2) or near(x[0], top3))
      corners = ( near(x[0], top1) or near(x[0], top2) or near(x[0], top3))
      bottom_ext = (near(x[1], ground_len-outer_pipe_len) )
      bottom_inner_pipe = near(x[1], ground_len-inner_pipe_len)
      return (corners or bottom_inner_pipe or bottom_ext) and on_boundary
walls = walls_boundary()

class vertical_simmetry_fluid(SubDomain):
    def inside(self, x, on_boundary):
      top0 = 0
      return (near(x[0], top0) and on_boundary)
simmetry  = vertical_simmetry_fluid()
#
#mark boundaries
#
#
inflow.mark(boundary_fluid, 2)
outflow.mark(boundary_fluid, 3)
walls.mark(boundary_fluid, 4)
simmetry.mark(boundary_fluid, 5)

# symmetry
# velocity in x-axys should be 0,and y-axis is open


bcu_inflow = DirichletBC(W.sub(0), Expression(inflow_profile, degree=2), boundary_fluid,2)
bcp_outflow = DirichletBC(W.sub(1), Constant(fluid_pressure), boundary_fluid,3)
bcp_outflow1 = DirichletBC(W.sub(1), Constant(0), boundary_fluid,2)
bcu_walls = DirichletBC(W.sub(0), Constant((0, 0)), boundary_fluid,4)
bcu_simmetry = DirichletBC(W.sub(0).sub(0), Constant(0), boundary_fluid,5)
bc_fluid = [bcu_inflow,bcu_walls,bcu_simmetry,bcp_outflow]

#plot(fluid)

# SOLID ##########################################################################3

# spaces for the thermal problem
Q1 = FunctionSpace(mesh, 'P', 1)

#
# Define boundaries --------------for the thermal problem----------
#
ground_bc = DirichletBC(Q1, Constant(ground_t), boundary_markers, 101)
inlet_bc = DirichletBC(Q1, Constant(fluid_t), boundary_markers, 102)

bcs = [ground_bc,inlet_bc]

# All other boundaries are isolated

class Tini(UserExpression):
    def __init__(self, markers, **kwargs):
        self.markers = markers
        super().__init__(**kwargs)
    def eval_cell(self, values, x, cell):
        if self.markers[cell.index] == 2:
            values[0] = ground_t #fluid_t
        else:
            values[0] = ground_t


#
# Variational form
#

T = TrialFunction(Q1) # trial and test functions for heat eq
ww = TestFunction(Q1)

T_init = Tini(subdomain_markers, degree=1)
T_old = interpolate( T_init, Q1 )
p_init = plot(T_old,title='Initialization', mode='color')

# Navier–Stokes #####################################################

# Define expressions used in variational forms
k  = Constant(dt)
mu = Constant(mu_fluid)

nu = mu/rho_fluid

# Navier stokes equation ###############################################################################

# Define body force with Boussinesq approximation
Boussinesq = (alpha/1000.)  * (T_old - T_init) * g
#
# Define functions for the current and previous time step
w = Function(W)  # current solution
u, p = split(w)
v, q = TestFunctions(W)

w_n = Function(W)  # Previous solution
u_n, p_n = split(w_n)#w_n.split()
#
nu = mu/rho_fluid     #    kinematic viscosity
#
U      = 0.5*(u_n + u)
F_momentum = (inner(u - u_n, v)/k + nu*inner(grad(u), grad(v)) + inner((-1./rho_fluid)*grad(p) + grad(u)*u, v) - div(u)*q)*dx_fluid
F_momentum += - Boussinesq * v[1] * dx_fluid
#F_momentum += Constant(DOLFIN_EPS)*p*q*dx_fluid
#
J=derivative(F_momentum, w)

# HEAT equation ####################################################################3

# Define variational problem for heat equation

L = rho_ground * c_p_ground * T_old * ww * dx(1) # Ground L
L += rho_inner_pipe * c_p_inner_pipe * T_old * ww * dx(3) # Internal pipe L
L += rho_fluid * c_p_fluid * T_old * ww * dx(2) # fluid L
L += dt * (heat_source/45.) * ww * dx(2) # Add load in IC only

A = rho_ground *c_p_ground *T*ww*dx(1) + dt*k_ground*dot(grad(T),grad(ww))*dx(1) # Ground A
A += rho_inner_pipe *c_p_inner_pipe *T*ww*dx(3) + dt*k_inner_pipe*dot(grad(T),grad(ww))*dx(3) # Internal Pipe A
A += rho_fluid * c_p_fluid * T * ww * dx(2) + dt * k_fluid * dot(grad(T), grad(ww)) * dx(2) + rho_fluid * c_p_fluid * dot(u_n, grad(T)) * ww * dx(2) # Fluid A
A += Constant(DOLFIN_EPS) * T * ww * dx

#
T_ = Function ( Q1 )


#
# Time loop #
#
t_init = 0.0
t_num = 10
dt = 0.01

error_u = 0.0
error_p =0.0
t = t_init
i = 0

for j in range ( 0, 7000 ):
  t = t + dt # update time
  print('j = ', j)
  #
  solve(F_momentum==0, w, bcs=bc_fluid, J=J
         ,solver_parameters={"newton_solver":{"linear_solver":'umfpack', "relative_tolerance":1E-3, "absolute_tolerance":1E-3, "maximum_iterations":100}}
         ,form_compiler_parameters={"cpp_optimize": True})
#
  v1, p1= w.split(True)
  v1n, p1n= w_n.split(True)
  if j>1:
    error_u = (v1.vector().max() - v1n.vector().max())/v1n.vector().max()
    error_p = (p1.vector().max() - p1n.vector().max())/p1n.vector().max()
#
  # Step 4: Heat equation
  solve ( A == L, T_, bcs )
  error_T = (T_.vector().max() - T_old.vector().max())/T_old.vector().max()
  print("[error u, error p , error T ] = ", (error_u,error_p,error_T))

  print('Re = ',rho_fluid*v1.vector().max()*outer_pipe_radius/mu_fluid)
  print("Max T: ", T_.vector().max())
#
  w_n.assign(w)
  T_old.assign ( T_ )

end_time = time.time()