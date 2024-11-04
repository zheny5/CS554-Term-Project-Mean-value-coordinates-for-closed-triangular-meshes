For the Color interpolation, 
Press 4, show the control mesh，then press 5, show the model after color interpolation.
If we want to take a look on the control mesh and its color,
press 6, then press 5.

For the deformation,
First we need to load the poly2. Just change the path of this_file2 to dog_control_mesh_deformed
Uncomment the line 96, which is poly->InterpolatePos(*poly2, *poly, meanValueCoordinates);
Other procedure is the same like the color interpolation.

To convert the byte .ply file to ASCII-code .ply file, a Python script is provided.

We can also change the file to bunny, bunny_control mesh(generated using convex hull) or cube, cube_control mesh(provided online)