Contents
Folders “Figures” and “Table_AIC_BIC_DF_LL_Deviance” contain supplementary material described in the manuscript Visini et al., 2024 EQSpectra.

Folders “Analyse_SUREvers2_database”, “Regressions” and “Model” contains matlab scripts.

License
This material is distributed under the Creative Common License Attribution-NonCommercial- ShareAlike 3.0 Unported (CC BY-NC-SA 3.0). You can share it with others as long as you provide proper credit, but you cannot change it in any way or use it commercially.

Disclaimer
The material shared here is distributed in the hope that it will be useful, but without any warranty: without even the implied warranty of merchantability or fitness for a particular purpose. While every precaution has been taken in the preparation of this material, in no event shall the authors be liable to any party for direct, indirect, special, incidental, or consequential damages, including lost profits, arising out of the use of information contained in this document or from the use of programs and source code that may accompany it, even if the authors have been advised of the possibility of such damage. The authors have no obligations to provide maintenance, support, updates, enhancements, or modifications.

This guide illustrates steps to analyze the SURE2.0 database, perform regressions and calculate PFDH curves.

1) Unzip file SURE-main.zip
2) Copy and paste the folder SURE-main (containing the folder SURE2.0_ruptures and the file SURE2.0_Slip_Obs_matlab.xlsx) into the folder "Analyse_SUREvers2_database"
Note that point 2) is necessary as the unzip creates extra sub-folder SURE-main.

The folder “Analyse_SUREvers2_database” contains conde to anaòyse the database and produce tables.

1)	A1_compute_r.m compute r- distances for normal and reverse from a list of events: list_Normal.txt and list_Revers.txt
2)	A2_union_table_r.m merge the tables calculated by A1_compute_r.m
3)	B1_interpolate_PF_coseismicslip.m and B2_interpolate_R15_coseismicslip.m calculate throw values on the PF and on the R1.5 from values in the database
4)	C_compute_s.m calculate the s-distances
5)	P1_prepare_shapefile_distanze_r.m and P2_prepare_shapefile_distanze_s.m create shapefiles of r and s distances (segments) to be visualized on a GIS

The folder “Regressions” contains codes to perform regressions, calculate median and standard deviation of the expected throw, calculate F ratio and perform Monte Carlo simulations.

1)	A1_Script_LOGISTICmultislicedimensions_20231019 calculates the parametres of the logistic regressions.
2)	B1_script_fitlme_residual_20231006.m computes regression to estimate median throw and standard deviation of DR
3)	C_ratio_DR2lenght_vs_PFlenght_20240729.m calculates the ratio between DR and PF lengths


The folder “Model” contains codes to perform a forward modelling.

1)	A_Montecarloprobability_at_the_site_20231023.m calculates the P of DR rank 2 occurrence at site using a montecarlo approach
2)	B_build_logistic_Comb_A.m, B_build_logistic_Comb_B.m and B_build_logistic_Comb_C.m calculate probailities of DR occurrence for a specific combination
3)	C_Calc_ThrowPFmean.m calculate ThrowPFmean
4)	script_combA_pfdhcurves.m script_combB_pfdhcurves.m and script_combC_pfdhcurves.m calculate curves of probability of exceedance for a specific combination

