import pygmsh

geom = pygmsh.opencascade.Geometry(
  characteristic_length_min=0.02,
  characteristic_length_max=0.02,
  )

point1 = geom.add_point([-0.04,-0.04,-2.0])
point2 = geom.add_point([0.04,-0.04,-2.0])
point3 = geom.add_point([-0.04,0.04,-2.0])
point4 = geom.add_point([0.04,0.04,-2.0])

line1 = geom.add_line(point1,point2)
line2 = geom.add_line(point2,point4)
line3 = geom.add_line(point4,point3)
line4 = geom.add_line(point3,point1)

lines = [];
lines.append(line1); lines.append(line2); lines.append(line3); lines.append(line4); 
geom.set_transfinite_lines(lines,5)		# (n-1) segments per line
line_loop1 = geom.add_line_loop(lines)

rect = geom.add_plane_surface(line_loop1,holes=None)
geom.set_transfinite_surface(rect) 
geom.add_raw_code('Mesh.RecombineAll = 1;')

disk = geom.add_disk([0.0, 0.0, -2.0], 0.07)
rect2 = geom.add_plane_surface(line_loop1,holes=None)
diskwithhole = geom.boolean_difference([disk],[rect2])
union = geom.boolean_union([diskwithhole,rect])


geom.extrude(union, translation_axis=[0.0, 0.0, 2.0], num_layers=100, recombine=True,)  #10 layers in z direction
#geom.extrude(union, translation_axis=[0.0, 0.0, 0.1], num_layers=5, recombine=True,)  #10 layers in z direction


mesh = pygmsh.generate_mesh(geom)

import meshio
meshio.write("cylinder_msh2.vtk", mesh)
meshio.write("cylinder_msh2.msh",mesh,"gmsh2-ascii")

