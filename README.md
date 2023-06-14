Dr. Denis Tsygankov (2023)

Integrative Systems Biology Lab

Wallace H. Coulter Department of Biomedical Engineering, Georgia Institute of Technology and Emory University SOM

If you use any part of the scripts in this package, please cite:

Siarhei Hladyshau, Jorik P Stoop, Kosei Kamada, Shuyi Nie, Denis V Tsygankov

Spatiotemporal coordination of Rac1 and Cdc42 at the whole cell level during cell ruffling

doi: https://doi.org/10.1101/2023.03.31.535147

## Coupled_RDE-Morphodynamics_Model

crop_frame.m - crop cell mask from a large simulation domain

extapolate_C.m - extrapolate values of concentration (A) along the preimeter of the cell mask (Im)

generate_regular_polyhedron.m - generate a binary mask representing a regular polyhedron

geometry_factor.m - copmute probabilities for a geometry factor

laplacian_DT.m - compute discrete laplacian using 5-point stencil with no-flux boundary conditions

outline_8p.m - compute outline of a binary mask as an 8-connected object

U_matrix.m - supplemental finction for computing laplacian

return_frame_to_canvas.m - transform the cropped simulation domain to the large simulation domain

model1: four-component model of GTPase and inhibitor

model2: four-component model of GTPase and inhibitor with regulation of actin factor by GTPase activation rate

model3: six-component model of Rac1 activity induced by Cdc42 with regulation of actin factor by Rac1 activation rate

model4: six-component model of bidirectionally coupled Rac1 and Cdc42 with regulation of actin factor by Rac1 activation rate

model5: eight-component model of Rac1 and Cdc42 activity induced by the upstream effector with regulation of actin factor by Rac1 activation rate

model6: eight-component model of Rac1 and Cdc42 activity (with feedback from Cdc42 to Rac1) induced by the upstream effector with regulation of actin factor by Rac1 activation rate

run_2D_dynamic_cell_model1.m - main script for running simulation of model1

protrude_extrapolate_model1.m - function for computing protrusion probabilities in model1

retract_reduce_model1.m - function for computing retraction probabilities in model1

run_2D_dynamic_cell_model2.m - main script for running simulation of model2

protrude_extrapolate_diff_model2.m - function for computing protrusion probabilities in model2

retract_reduce_diff_model2.m - function for computing retraction probabilities in model2

run_2D_dynamic_cell_model3.m - main script for running simulation of model3

protrude_extrapolate_diff_model3.m - function for computing protrusion probabilities in model3

retract_reduce_diff_model3.m - function for computing retraction probabilities in model3

run_2D_dynamic_cell_model4.m - main script for running simulation of model4

protrude_extrapolate_diff_model4.m - function for computing protrusion probabilities in model4

retract_reduce_diff_model4.m - function for computing retraction probabilities in model4

run_2D_dynamic_cell_model5.m - main script for running simulation of model5

protrude_extrapolate_diff_model5.m - function for computing protrusion probabilities in model5

retract_reduce_diff_model5.m - function for computing retraction probabilities in model5

run_2D_dynamic_cell_model6.m - main script for running simulation of model6

protrude_extrapolate_diff_model6.m - function for computing protrusion probabilities in model6

retract_reduce_diff_model6.m - function for computing retraction probabilities in model6
