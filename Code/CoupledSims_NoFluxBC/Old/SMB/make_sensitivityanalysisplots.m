
% Set actin filament parameters
% 
Da      = 0.1;                  % diffusion coefficient for actin
m0      = 2.0;                  % competition for actin monomers
kon     = 0.0001;
koff    = 0.11;
kfb     = 0.1;

% Set polarity protein parameters
% 
N       = 200;                  % total number of molecules in the cell (conserved)
ron     = 0.001;                % units: 1/min
rfb     = 1.0;                  % units: 1/min
roff    = 0.9;                  % units: 1/min
D       = 0.012;

% sensitivity plot for alpha, beta
% 5 runs with different initial conditions
% 
alpha = [0.01,0.05,0.1,0.2,0.5];
beta = [0.01,0.05,0.1,0.5,1.0];
ccc = zeros(5,5);
ccc(1,1:5) = [0,0,0,0,0];


figure(10);
imagesc(alpha,beta,ccc./5)
set(gca,'YDir','normal');
set(gcf,'color','w');
ylabel('\beta','fontsize',20);
xlabel('\alpha','fontsize',20);
set(gca,'FontSize',20);
xlabel('\alpha','fontsize',20);
ylabel('\beta','fontsize',20);
colorbar
colormap(magma)
colormap(flipud(magma))