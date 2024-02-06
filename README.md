# STCutFEMOverlapMesh
Implementations presented in the two papers "Space-Time CutFEM on Overlapping Meshes I &amp; II".

The code is first structured into either one of the two papers, then into either using dG(0) or dG(1) in time. After this the structure is the same for all four cases. 

1) To run simulations for mesh and solution plots: run "main.m" for space-time visuals and "main_PSACT" for animations of spatial results. The domain, its discretization, and some other problem parameters are defined in a "parameters" m-file which should be uncommented in the "main" m-file. Other files for defining problem data are "mu_func.m", "f_func.m", and "u_func.m".

2) To perform error convergence studies: run "ErrConComp_k.m" for k-convergence and "ErrConComp_h.m" for h-convergence. In these files, general convergence study parameters are defined such as interval on the x-axis. Both files call "main_ECC.m" to solve the problem used in the study. This problem is defined just as in point 1) above.
