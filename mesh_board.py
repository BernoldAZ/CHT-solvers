import gmsh

# Initialize GMSH
gmsh.initialize()

# Add a new model
gmsh.model.add("Convection_Cooling_Circuit_Boards")

# Define parameters
board_length = 130  # mm
board_thickness = 2  # mm
ic_size = 20  # mm (length and width of IC)
ic_thickness = 2  # mm
ic_spacing = 10  # mm (vertical spacing between ICs)
num_ics = 4  # Number of ICs on the board
air_space_thickness = 10 #10

lc_finer = 0.3;  # Finer mesh around certain regions
lc_coarser = 1.0;   # Coarser mesh elsewhere

# Create the points
board_LD = gmsh.model.geo.addPoint(0, 0, 0)
board_RD_air_LD = gmsh.model.geo.addPoint(board_thickness, 0, 0)
board_RU_air_LU = gmsh.model.geo.addPoint(board_thickness, board_length, 0)
board_LU = gmsh.model.geo.addPoint(0, board_length, 0)

air_RD = gmsh.model.geo.addPoint(board_thickness + air_space_thickness, 0, 0)
air_RU = gmsh.model.geo.addPoint(board_thickness + air_space_thickness, board_length, 0)

ic1_LD = gmsh.model.geo.addPoint(board_thickness, 10, 0)
ic1_RD = gmsh.model.geo.addPoint(board_thickness + ic_thickness, 10, 0)
ic1_LU = gmsh.model.geo.addPoint(board_thickness + ic_thickness, 30, 0)
ic1_LU = gmsh.model.geo.addPoint(board_thickness, 30, 0)

ic2_LD = gmsh.model.geo.addPoint(board_thickness, 40, 0)
ic2_RD = gmsh.model.geo.addPoint(board_thickness + ic_thickness, 40, 0)
ic2_LU = gmsh.model.geo.addPoint(board_thickness + ic_thickness, 60, 0)
ic2_LU = gmsh.model.geo.addPoint(board_thickness, 60, 0)

ic3_LD = gmsh.model.geo.addPoint(board_thickness, 70, 0)
ic3_RD = gmsh.model.geo.addPoint(board_thickness + ic_thickness, 70, 0)
ic3_LU = gmsh.model.geo.addPoint(board_thickness + ic_thickness, 90, 0)
ic3_LU = gmsh.model.geo.addPoint(board_thickness, 90, 0)

ic4_LD = gmsh.model.geo.addPoint(board_thickness,100, 0,lc_finer)
ic4_RD = gmsh.model.geo.addPoint(board_thickness + ic_thickness, 100, 0)
ic4_LU = gmsh.model.geo.addPoint(board_thickness + ic_thickness, 120, 0)
ic4_LU = gmsh.model.geo.addPoint(board_thickness, 120, 0)


# Create the lines

# board
left_board = gmsh.model.geo.addLine(1, 4)  #
up_board = gmsh.model.geo.addLine(4, 3)  #
#right part
right1 = gmsh.model.geo.addLine(3, 22)  #
right2 = gmsh.model.geo.addLine(22, 19)  #
right3 = gmsh.model.geo.addLine(19, 18)  #
right4 = gmsh.model.geo.addLine(18, 15)  #
right5 = gmsh.model.geo.addLine(15, 14)  #
right6 = gmsh.model.geo.addLine(14, 11)  #
right7 = gmsh.model.geo.addLine(11, 10)  #
right8 = gmsh.model.geo.addLine(10, 7)  #
right9 = gmsh.model.geo.addLine(7, 2)  #

down_board = gmsh.model.geo.addLine(2, 1)  #

# air
up_air = gmsh.model.geo.addLine(6, 3)
# right1 to right9
down_air = gmsh.model.geo.addLine(2, 5)
right_air = gmsh.model.geo.addLine(5, 6)

# ICs - from down to up
# 1
ic1_down = gmsh.model.geo.addLine(7, 8)  #
ic1_right = gmsh.model.geo.addLine(8, 9)  #
ic1_up = gmsh.model.geo.addLine(9, 10)  #
# right8

# 2
ic2_down = gmsh.model.geo.addLine(11, 12)  #
ic2_right = gmsh.model.geo.addLine(12, 13)  #
ic2_up = gmsh.model.geo.addLine(13, 14)  #
# right6

# 3
ic3_down = gmsh.model.geo.addLine(15, 16)  #
ic3_right = gmsh.model.geo.addLine(16, 17)  #
ic3_up = gmsh.model.geo.addLine(17, 18)  #
# right4

# 4
ic4_down = gmsh.model.geo.addLine(19, 20)  #
ic4_right = gmsh.model.geo.addLine(20, 21)  #
ic4_up = gmsh.model.geo.addLine(21, 22)  #
# right2

# Create loops and surfaces

air_lines = [up_air,right1,right2,right3,right4,right5,right6,right7,right8,right9,down_air,right_air]
air_loop = gmsh.model.geo.addCurveLoop(air_lines, 1)

board_lines = [left_board,up_board,right1,right2,right3,right4,right5,right6,right7,right8,right9,down_board]
board_loop = gmsh.model.geo.addCurveLoop(board_lines, 2)

ic1_lines = [ic1_down,ic1_right,ic1_up,right8]
ic1_loop = gmsh.model.geo.addCurveLoop(ic1_lines, 3)

ic2_lines = [ic2_down,ic2_right,ic2_up,right6]
ic2_loop = gmsh.model.geo.addCurveLoop(ic2_lines, 4)

ic3_lines = [ic3_down,ic3_right,ic3_up,right4]
ic3_loop = gmsh.model.geo.addCurveLoop(ic3_lines, 5)

ic4_lines = [ic4_down,ic4_right,ic4_up,right2]
ic4_loop = gmsh.model.geo.addCurveLoop(ic4_lines, 6)

board_surface = gmsh.model.geo.addPlaneSurface([board_loop])
air_surface = gmsh.model.geo.addPlaneSurface([air_loop,ic1_loop,ic2_loop,ic3_loop,ic4_loop])
ic1_surface = gmsh.model.geo.addPlaneSurface([ic1_loop])
ic2_surface = gmsh.model.geo.addPlaneSurface([ic2_loop])
ic3_surface = gmsh.model.geo.addPlaneSurface([ic3_loop])
ic4_surface = gmsh.model.geo.addPlaneSurface([ic4_loop])

# Physical surfaces
gmsh.model.addPhysicalGroup(2, [board_surface], 1)
gmsh.model.setPhysicalName(2, 1, "board")
#
gmsh.model.addPhysicalGroup(2, [air_surface], 2)
gmsh.model.setPhysicalName(2, 2, "air")
#
gmsh.model.addPhysicalGroup(2, [ic1_surface,ic2_surface,ic3_surface,ic4_surface], 3)
gmsh.model.setPhysicalName(2, 3, "ic")

# A physical lineS
#
#
gmsh.model.addPhysicalGroup(1, [down_air], 201)
gmsh.model.setPhysicalName(1, 201, "air_inlet")

# Synchronize and generate the mesh
gmsh.model.geo.synchronize()

# Define mesh size (optional)
gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc_coarser)  # Coarser mesh
gmsh.model.mesh.setSize([(0, 2),(0, 7),(0, 8),(0, 9),(0, 10),(0, 11),(0, 12),(0, 13),(0, 14),(0, 15),(0, 16),(0, 17),(0, 18),(0, 19),(0, 20),(0, 21),(0, 22),(0, 3)], lc_finer)
gmsh.model.mesh.generate(2)

# Save the mesh to a file
gmsh.write("Convection_Cooling_Circuit_Boards.msh")

# Finalize the gmsh API
gmsh.finalize()