% generate test data for copulastat.py

clear;
clc;

% remove the old copulastat_test.mat
delete('copulastat_test.mat')

gauss_ktau_rho_0_3 = copulastat('Gaussian', 0.3, 'type', 'kendall');
gauss_srho_rho_0_3 = copulastat('Gaussian', 0.3, 'type', 'spearman');
gauss_ktau_rho_0_7 = copulastat('Gaussian', 0.7, 'type', 'kendall');
gauss_srho_rho_0_7 = copulastat('Gaussian', 0.7, 'type', 'spearman');
gauss_ktau_rho_1_0 = copulastat('Gaussian', 1.0, 'type', 'kendall');
gauss_srho_rho_1_0 = copulastat('Gaussian', 1.0, 'type', 'spearman');

clayton_ktau_alpha_0_3 = copulastat('clayton', 0.3, 'type', 'kendall');
clayton_srho_alpha_0_3 = copulastat('clayton', 0.3, 'type', 'spearman');
clayton_ktau_alpha_0_7 = copulastat('clayton', 0.7, 'type', 'kendall');
clayton_srho_alpha_0_7 = copulastat('clayton', 0.7, 'type', 'spearman');
clayton_ktau_alpha_1_0 = copulastat('clayton', 1.0, 'type', 'kendall');
clayton_srho_alpha_1_0 = copulastat('clayton', 1.0, 'type', 'spearman');

gumbel_ktau_alpha_1_0 = copulastat('gumbel', 1.0, 'type', 'kendall');
gumbel_srho_alpha_1_0 = copulastat('gumbel', 1.0, 'type', 'spearman');
gumbel_ktau_alpha_3_0 = copulastat('gumbel', 3.0, 'type', 'kendall');
gumbel_srho_alpha_3_0 = copulastat('gumbel', 3.0, 'type', 'spearman');

frank_ktau_alpha_0_3 = copulastat('frank', 0.3, 'type', 'kendall');
frank_srho_alpha_0_3 = copulastat('frank', 0.3, 'type', 'spearman');
frank_ktau_alpha_0_7 = copulastat('frank', 0.7, 'type', 'kendall');
frank_srho_alpha_0_7 = copulastat('frank', 0.7, 'type', 'spearman');
frank_ktau_alpha_1_0 = copulastat('frank', 1.0, 'type', 'kendall');
frank_srho_alpha_1_0 = copulastat('frank', 1.0, 'type', 'spearman');

save('copulastat_test.mat');