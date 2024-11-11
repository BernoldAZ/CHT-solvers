import gmsh

# Parameters based on domain scheme
inner_pipe_radius = 0.020*5      # 40 mm diameter for inner pipe
inner_pipe_len = 165/5 # Effective length of the internal pipe (in meters)
outer_pipe_radius = 0.0575*5  # 0.115 borehole mm diameter
outer_pipe_len = inner_pipe_len+2*inner_pipe_radius  #
ground_len = 170/5 # Effective length of the external pipe (in meters)
ground_width = 5

lc_finer = 0.08;  # Finer mesh around certain regions
lc_less_finer = 0.08;  # Finer mesh around certain regions
lc_coarser = 1.0;   # Coarser mesh elsewhere

# Initialize GMSH
gmsh.initialize()

# Add a new model
gmsh.model.add("Coaxial_Borehole_Heat_Exchanger")

# Define mesh size (number of elements)

# Create the points
# Ground
ground1_LD = gmsh.model.geo.addPoint(0, 0, 0)
ground1_LU = gmsh.model.geo.addPoint(0, ground_len, 0)
ground1_RU = gmsh.model.geo.addPoint(ground_width, ground_len, 0)
ground1_RD = gmsh.model.geo.addPoint(ground_width, 0, 0)

# Outer pipe
out_pipe_LD = gmsh.model.geo.addPoint(0, ground_len-outer_pipe_len, 0)
out_pipe_LU = gmsh.model.geo.addPoint(0, ground_len, 0)
out_pipe_RU = gmsh.model.geo.addPoint(2*outer_pipe_radius, ground_len, 0)
out_pipe_RD = gmsh.model.geo.addPoint(2*outer_pipe_radius, ground_len-outer_pipe_len, 0)

# inner pipe
inner_pipe_LD = gmsh.model.geo.addPoint(outer_pipe_radius-inner_pipe_radius, ground_len-inner_pipe_len, 0)
inner_pipe_LU = gmsh.model.geo.addPoint(outer_pipe_radius-inner_pipe_radius, ground_len, 0)
inner_pipe_RU = gmsh.model.geo.addPoint(outer_pipe_radius+inner_pipe_radius, ground_len, 0)
inner_pipe_RD = gmsh.model.geo.addPoint(outer_pipe_radius+inner_pipe_radius, ground_len-inner_pipe_len, 0)


# Create the lines

#Ground
left_ground = gmsh.model.geo.addLine(1, 5)  #
left_ext_pipe = gmsh.model.geo.addLine(5, 2)  #
top_ext_pipe1 = gmsh.model.geo.addLine(2, 10)  #
top_int_pipe = gmsh.model.geo.addLine(10, 11)  #
top_ext_pipe2 = gmsh.model.geo.addLine(11, 7)  #
top_ground = gmsh.model.geo.addLine(7, 3)  #
right_ground = gmsh.model.geo.addLine(3, 4)  #
bottom_ground = gmsh.model.geo.addLine(4, 1)  #

#External pipe
# left_ext_pipe
# top_ext_pipe1
# top_int_pipe
# top_ext_pipe2
right_ext_pipe = gmsh.model.geo.addLine(7, 8)  #
bottom_ext_pipe = gmsh.model.geo.addLine(8, 5)  #

#Internal pipe
# top_int_pipe
right_int_pipe = gmsh.model.geo.addLine(11, 12)  #
bottom_int_pipe = gmsh.model.geo.addLine(12, 9)  #
left_int_pipe = gmsh.model.geo.addLine(9, 10)  #

# Add curves

ground_lines = [left_ground,left_ext_pipe,top_ext_pipe1,top_int_pipe,top_ext_pipe2,top_ground,right_ground,bottom_ground]
ground_loop = gmsh.model.geo.addCurveLoop(ground_lines, 1)

ext_pipe_lines = [left_ext_pipe,top_ext_pipe1,top_int_pipe,top_ext_pipe2,right_ext_pipe,bottom_ext_pipe]
ext_pipe_loop = gmsh.model.geo.addCurveLoop(ext_pipe_lines, 2)

int_pipe_lines = [top_int_pipe,right_int_pipe,bottom_int_pipe,left_int_pipe]
int_pipe_loop = gmsh.model.geo.addCurveLoop(int_pipe_lines, 3)

# Add surfaces

ground_surface = gmsh.model.geo.addPlaneSurface([ground_loop,ext_pipe_loop])
ext_pipe_surface = gmsh.model.geo.addPlaneSurface([ext_pipe_loop,int_pipe_loop])
int_pipe_surface = gmsh.model.geo.addPlaneSurface([int_pipe_loop])

# Physical surfaces
gmsh.model.addPhysicalGroup(2, [ground_surface], 1)
gmsh.model.setPhysicalName(2, 1, "ground")
#
gmsh.model.addPhysicalGroup(2, [ext_pipe_surface], 2)
gmsh.model.setPhysicalName(2, 2, "fluid")
#
gmsh.model.addPhysicalGroup(2, [int_pipe_surface], 3)
gmsh.model.setPhysicalName(2, 3, "int_pipe")

# A physical lineS
#
#
gmsh.model.addPhysicalGroup(1, [right_ground,bottom_ground], 101)
gmsh.model.setPhysicalName(1, 101, "ground_HT_boundary")
#
gmsh.model.addPhysicalGroup(1, [top_ext_pipe1], 102)
gmsh.model.setPhysicalName(1, 102, "inlet_HT_boundary")
#
gmsh.model.addPhysicalGroup(1, [top_ext_pipe2], 103)
gmsh.model.setPhysicalName(1, 102, "outlet_HT_boundary")
#
#gmsh.model.addPhysicalGroup(1, [top_int_pipe,top_ground], 103)
#gmsh.model.setPhysicalName(1, 103, "top_boundary")

# Synchronize and generate the mesh
gmsh.model.geo.synchronize()

# Define mesh size (optional)
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc_coarser)  # Coarser mesh
gmsh.model.mesh.setSize([(0, 10),(0, 11),(0, 2),(0, 7)], lc_less_finer)
gmsh.model.mesh.setSize([(0, 9),(0, 12),(0, 5),(0, 8)], lc_finer)
gmsh.model.mesh.generate(2)

# Save the mesh to a file
gmsh.write("Coaxial_Borehole_Heat_Exchanger.msh")

# Finalize the gmsh API
gmsh.finalize()