% Matlab test script which generates copula samples similar to copulacdf.py
% for comparison purposes

% Generate samples of gaussian copula
u = linspace(0+eps,1-eps,10);
d = 2;
[U1,U2] = meshgrid(u,u);
rho = 0.8;
Rho = [1 rho; rho 1];
gaussian_copula_cdf = copulacdf('Gaussian',[U1(:) U2(:)], Rho);

% Generate samples of other coupla's here

% save them all for testing against python generated data
save('copula_cdf_test.mat', 'gaussian_copula_cdf')
