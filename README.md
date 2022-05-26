# Sums-of-Heterogeneous-Quadratics

------------------------------------------------------------------------------------------

This repository contains all source code and files necessary to reproduce the experiments of 

``A Semidefinite Relaxation for Sums of Heterogeneous Quadratics on the Stiefel Manifold.”

------------------------------------------------------------------------------------------

Organization:

The following MATLAB files produce the experiments in the tables and figures of the main paper and supplement:

shq_sdp_hppca_extended_experiments.m : Generate and save data for Tables 1 and 2.

shq_sdp_randomExperiments_master.m : Generate and save data to test SDP tightness in more problem instances, including the experiment in Table 3.

read_experiment_results.m: Output tables from shq_sdp_hppca_extended_experiments.m (Tables 1 and 2) and shq_sdp_randomExperiments_master.m (Table 3)



shq_orderedASD_master.m : Generate experiments from Figure 1.

shq_hppca_master.m : Generate experiments from Figure 2.

feasibilityCheck_timing.m : Generate experiments from Figure 3. 

------------------------------------------------------------------------------------------

Note that some of the experiments take many hours to compute. Hence, all data generated from our experiments are included in files with the extension .mat, with the source code provided to plot their results.