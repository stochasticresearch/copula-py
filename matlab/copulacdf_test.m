% Matlab test script which generates copula samples similar to copulacdf.py
% for comparison purposes

% remove the old copulacdf_test.mat
delete('copulacdf_test.mat')

% data which will define where we want to know the value of the Copula
u = linspace(0+eps,1-eps,10);
d = 2;
[U1,U2] = meshgrid(u,u);

% Generate samples of Gaussian copula
rho = 0.8;
Rho = [1 rho; rho 1];
gaussian_copula_cdf = copulacdf('gaussian',[U1(:) U2(:)], Rho);

% Generate samples of T copula

% Generate samples of the Clayton copula
alpha = 0.3;
clayton_copula_cdf = copulacdf('clayton',[U1(:) U2(:)], alpha);

% Generate samples of the Frank copula
frank_copula_cdf = copulacdf('frank',[U1(:) U2(:)], alpha);

% Generate samples of the Gumbel Copula
alpha = 1.5;
gumbel_copula_cdf = copulacdf('gumbel',[U1(:) U2(:)], alpha);

% save them all for testing against python generated data
save('copulacdf_test.mat', ...
        'gaussian_copula_cdf', ...
        'clayton_copula_cdf', ...
        'frank_copula_cdf', ...
        'gumbel_copula_cdf')
