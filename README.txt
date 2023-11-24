%% readme file for scripts in "Analyse_SUREvers2_database" and "Regressions"

% readme file for Analyse_SUREvers2_database.
All the codes use directly both kinematics (Reverse and Normal) within a loop
workflow for the codes is: 
A1_compute_r.m
A2_union_table_r.m
B1_interpolate_PF_coseismicslip.m
B2_interpolate_R15_coseismicslip.m
C_compute_s.m
P1_prepare_shapefile_distanze_r.m,
P2_prepare_shapefile_distanze_s.m

Some notes are given below.

1) Unzip file SURE-main.zip

2) Copy and paste the folder SURE-main (containing the folder SURE2.0_ruptures and the file SURE2.0_Slip_Obs_matlab.xlsx) into the folder "Analyse_SUREvers2_database"
 Note that point 2) is necessary as the unzip creates extra sub-folder SURE-main.

3)  In B1 we assign Throw = 0 m to the tips. Tips are the two most distant tips define the rupture length (lines 84-89).
We here use scatteredInterpolant to perform interpolation on a 2-D or 3-D data set of scattered data. 
In the future one could test least cost path function for defining PF trace and interpolates along it.

% readme file for Regression.
All the codes use directly both kinematics (Reverse and Normal) within a loop
workflow for the codes is: 
A1_Script_LOGISTICmultislicedimensions_20231019.m
B1_script_fitlme_residual_20231006.m
B1_interpolate_PF_coseismicslip.m
B2_interpolate_R15_coseismicslip.m
C_ratio_DR2lenght_vs_PFlenght.m
D_Montecarloprobability_at_the_site_20231023.m
E_build_logistic_using_sitedimension.m




