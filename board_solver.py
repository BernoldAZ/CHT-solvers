
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

#Board
rho_board = 1.25 * 10**(-9) # density: 1250kg/m^3 ->  1.25 × 10^-9 ton/mm^3
c_p_board = 1300 # heat capacity: 1300J/kg*K -> 1300  mJ/ton*K
k_board = 0.35 # Thermal conductivity: 0.35W/m*k -> 0.35 mW/mm*K

#IC
rho_ic =  2.3 * 10**(-9) # density: 2300kg/m^3 -> 2.3 × 10^-9 ton/mm^3
c_p_ic = 710 # heat capacity: 710J/kg*K -> 710  mJ/ton*K
k_ic = 150 # Thermal conductivity: 150W/m*K -> 150 mW/mm*K

#Air
rho_air =  1.16*10**(-12) # density: 1.16kg/m^3 -> 1.16 x 10^-12 ton/mm^3
c_p_air = 1.0082 # heat capacity: 1008.2J/kg*K -> 1.0082 mJ/ton*K
k_air = 0.0282 # Thermal conductivity: 0.0282W/m*K -> 0.0282 mW/mm*K
mu_air =  1.78 * 10**(-11) # dynmic viscosity: 1.78*10^-5kg/m*s ->   1.78 × 10^-11 ton/mm*s

#  Define time
t_init = 0.0
t_final = 10.0
t_num = 404
dt = ( t_final - t_init ) / t_num

# Define External force

alpha = 0.0034 # Thermal expasion coeffcient: 0.0034/K
g = 9810 # Graviational acceleration: 9.81 m/s^2 -> 9810 mm/s^2
heat_source = Constant ( 0.83333 ) # The heat source term: 8.3333*10^5 W/m^3 -> 0.83333 mW/mm^3
f_f  = Constant((0, 0)) # External force field applied to the air

#inflow_profile = ('0.', '50.')    # Velocity entering the model Use it for forced convection
v0 = Constant( (0.0, 0.0) )
air_pressure = 0.10133  # Atmospheric pressure: 1.0133e5  Pa -> 0.10133 MPa
inlet_temperature = Constant(300.0) # Temperature of the entering fluid

# Define expressions used in variational forms
k  = Constant(dt)
mu = Constant(mu_air)
nu = mu/rho_air     #    kinematic viscosity

#
# Define Load mesh
#

# READS MSH FILE AND CREATE XDMF FILES
mesh_from_file = meshio.read("Convection_Cooling_Circuit_Boards.msh")
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

subdomain_markers = mf2
boundary_markers = mf

## Step 1.5: Create submeshes and ds/dx 

# Extract SubMesh for physical group tag
air = SubMesh(mesh, subdomain_markers, 2)

# All mesh
dx=Measure("dx", domain=mesh, subdomain_data=subdomain_markers)
ds=Measure("ds", domain=mesh, subdomain_data=boundary_markers)

# air
boundary_air = MeshFunction("size_t",air,air.topology().dim()-1,0)
dx_air=Measure("dx", domain=air, subdomain_data=boundary_air)
ds_air=Measure("ds", domain=air, subdomain_data=boundary_air)

#
# Define boundaries
#

# FLUID

# spaces for the navier stokes problem
V = VectorElement("CG", triangle, 2)
Q = FiniteElement("CG", triangle, 1)
W = FunctionSpace(air, V * Q)

# Define boundaries 
class inflow_fluid(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.) and on_boundary
inflow = inflow_fluid()

class outflow_fluid(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 130.) and on_boundary
outflow  = outflow_fluid()

class walls_air(SubDomain):
    def inside(self, x, on_boundary):
      left_right = ( x[1] > 0 and x[1] < 130 and on_boundary)
      corners = (near(x[0], 2.) or near(x[0], 12.)) and on_boundary
      return  left_right or corners
walls = walls_air()

#mark boundaries
inflow.mark(boundary_air, 2)
outflow.mark(boundary_air, 3)
walls.mark(boundary_air, 4)

#bcu_inflow = DirichletBC(W.sub(0), Constant(inflow_profile), boundary_air,2) # Use it for forced convection
bcp_outflow = DirichletBC(W.sub(1), Constant(air_pressure), boundary_air,3)  
bcu_walls = DirichletBC(W.sub(0), v0, boundary_air,4)
#
#bc_air = [bcu_inflow,bcu_walls,bcp_outflow] # Use it for forced convection
bc_air = [bcp_outflow,bcu_walls]
#

# SOLID 

class PeriodicBoundary(SubDomain):
    def __init__(self):
        SubDomain.__init__(self)

    # Map the left boundary to the right boundary
    def inside(self, x, on_boundary):
        # Adjust the condition for "inside" to match left boundary
        return near(x[0], 0) and on_boundary

    def map(self, x, y):
        # Map boundary
        y[0] = x[0] + (board_thickness +air_space_thickness) # Constants from gmsh
        y[1] = x[1]

pbc = PeriodicBoundary()

# spaces for the thermal problem
Q1 = FunctionSpace(mesh, 'P', 1, constrained_domain=pbc)

#mark boundaries
bcs = DirichletBC(Q1, inlet_temperature, boundary_markers,201)
# All other boundaries are isolated

# Initial state
class Tini(UserExpression):
    def __init__(self, markers, **kwargs):
        self.markers = markers
        super().__init__(**kwargs)
    def eval_cell(self, values, x, cell):
        if self.markers[cell.index] == 3:
            values[0] = 300.0 
        else:
            values[0] = 300.0 

#
# Variational form
#

# trial and test functions for heat eq
T = TrialFunction(Q1)
ww = TestFunction(Q1)

# Temperature vector in previous step
T_init = Tini(subdomain_markers, degree=1)
T_old = interpolate( T_init, Q1 )

# Navier–Stokes #

# Define body force with Boussinesq approximation
Boussinesq = (alpha/2090)  * (T_old - T_init) * g   # Alpha reduced to easy convergence

# Define functions for the current and previous time step
# current solution
w = Function(W)
u, p = split(w)
v, q = TestFunctions(W)

# Previous solution
w_n = Function(W)
u_n, p_n = split(w_n)

# Variational form of navier stokes
F_momentum = (inner(u - u_n, v)/k + nu*inner(grad(u), grad(v)) + inner((-1./rho_air)*grad(p) + grad(u)*u, v) + div(u)*q)*dx_air
F_momentum +=  Boussinesq * v[1] * dx_air

J=derivative(F_momentum, w)

# Heat equation #

# Define variational problem for heat equation
# linear form
L = rho_board * c_p_board * T_old * ww * dx(1) # Board L
L += rho_ic * c_p_ic * T_old * ww * dx(3) # IC L
L += rho_air * c_p_air * T_old * ww * dx(2) # AIR L
L += dt * (heat_source/45.0) * ww * dx(3) # Add load in IC only

# Bilinear form
A = rho_board *c_p_board *T*ww*dx(1) + dt*k_board*dot(grad(T),grad(ww))*dx(1) # Board A
A += rho_ic *c_p_ic *T*ww*dx(3) + dt*k_ic*dot(grad(T),grad(ww))*dx(3) # IC A
A += rho_air * c_p_air * T * ww * dx(2) + dt * k_air * dot(grad(T), grad(ww)) * dx(2) + dt * rho_air * c_p_air * inner(dot(u_n, grad(T)) , ww) * dx(2) # AIR A

# esults
T_ = Function ( Q1 )

error_u = 0.0
error_p =0.0
t = t_init
i = 0
tol = 1E-5

start_time = time.time()


#
# Time loop #
#

for j in range ( 0, t_num + 1 ):
  
  # update time
  t = t + dt

  # solve navier stokes
  solve(F_momentum==0, w, bcs=bc_air, J=J
         ,solver_parameters={"newton_solver":{"linear_solver":'umfpack', "relative_tolerance":1E-3, "absolute_tolerance":1E-3, "maximum_iterations":50}}
         ,form_compiler_parameters={"cpp_optimize": True})
  
  # get solutions at current step
  v1, p1= w.split(True)
  v1n, p1n= w_n.split(True)
  if j>1:
    error_u = (v1.vector().max() - v1n.vector().max())/v1n.vector().max()
    error_p = (p1.vector().max() - p1n.vector().max())/p1n.vector().max()

  # solve Heat equation
  solve ( A == L, T_, bcs )
  error_T = (T_.vector().max() - T_old.vector().max())/T_old.vector().max()

  # Update previous solution
  w_n.assign(w)
  T_old.assign ( T_ )

  # Early stopping condition
  if error_u <= tol and error_p <= tol and error_T <= tol:
    break # Stops the program if there is no change in u,p,T

end_time = time.time()

