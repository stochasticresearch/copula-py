% Matlab test script which generates copula samples similar to copulapdf.py
% for comparison purposes

% remove the old copulapdf_test.mat
delete('copulapdf_test.mat')

% data which will define where we want to know the value of the Copula
u = linspace(0.1,0.9,10);
d = 2;
[U1,U2] = meshgrid(u,u);

% Generate samples of Gaussian copula
rho = 0.8;
Rho = [1 rho; rho 1];
gaussian_copula_pdf = copulapdf('gaussian',[U1(:) U2(:)], Rho);

% Generate samples of T copula
nu = 2;
t_copula_pdf = copulapdf('t',[U1(:) U2(:)], Rho, nu);

u = linspace(0+eps,1-eps,10);
d = 2;
[U1,U2] = meshgrid(u,u);

% Generate samples of the Clayton copula
alpha = 0.3;
clayton_copula_pdf = copulapdf('clayton',[U1(:) U2(:)], alpha);

% Generate samples of the Frank copula
frank_copula_pdf = copulapdf('frank',[U1(:) U2(:)], alpha);

% Generate samples of the Gumbel Copula
alpha = 1.5;
gumbel_copula_pdf = copulapdf('gumbel',[U1(:) U2(:)], alpha);

% save them all for testing against python generated data
save('copulapdf_test.mat', ...
        'gaussian_copula_pdf', ...
        't_copula_pdf', ...
        'clayton_copula_pdf', ...
        'frank_copula_pdf', ...
        'gumbel_copula_pdf')
