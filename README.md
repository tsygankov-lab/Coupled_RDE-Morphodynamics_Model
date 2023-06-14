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

inspect_activator_time_derivative.m - analyze the values of kinetic rates from simulation data

inspect_data.m - analyze simulation data

laplacian_DT.m - compute discrete laplacian using 5-point stencil with no-flux boundary conditions

make_movies_all.m - generate movie from a timeseries images data

outline_8p.m - compute outline of a binary mask as an 8-connected object

plot_colored_all_cropped_outline.m - plot the results of simulation

protrude_extrapolate.m - perform protrusion step in the simulation and extrapolate the values of the model components in the protruded pixels

protrude_extrapolate_diff.m - perform protrusion step in the simulation based on the differentiator assumption (regulation by the rate of concentration change) and extrapolate the values of the model components in the protruded pixels

retract_reduce.m - perform retraction step in the simulation and reduce the values of the model components in the retracted pixels (distribute the concentration along the whole simulation domain)

retract_reduce_diff.m - perform retraction step in the simulation based on the differentiator assumption (regulation by the rate of concentration change) and reduce the values of the model components in the retracted pixels (distribute the concentration along the whole simulation domain)

return_frame_to_canvas.m - transform the cropped simulation domain to the large simulation domain

run_2D_dynamic_cell_GTPase_inhibitor.m - main script for simulation of cell membrane ruffling, four-component model of GTPase and inhibitor

run_2D_dynamic_cell_GTPase_inhibitor_differentiator.m - main script for simulation of the four-component model of GTPase and inhibitor with regulation of actin factor by GTPase activation rate

run_2D_dynamic_cell_Cdc42_act_Rac1_differentiator.m - main script for simulation of the six-component model of Rac1 activity induced by Cdc42 with regulation of actin factor by Rac1 activation rate

run_2D_dynamic_cell_Cdc42_coupled_Rac1_differentiator.m - main script for simulation of the six-component model of bidirectionally coupled Rac1 and Cdc42 with regulation of actin factor by Rac1 activation rate

run_2D_dynamic_cell_Cdc42_Rac1_up_reg_differentiator.m - main script for simulation of the eight-component model of Rac1 and Cdc42 activity induced by the upstream effector with regulation of actin factor by Rac1 activation rate

run_2D_dynamic_cell_Cdc42_Rac1_up_reg_feedback_differentiator.m - main script for simulation of the eight-component model of Rac1 and Cdc42 activity (with feedback from Cdc42 to Rac1) induced by the upstream effector with regulation of actin factor by Rac1 activation rate

U_matrix.m - supplemental finction for computing laplacian

