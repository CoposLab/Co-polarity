Da      = 0.5;                  % diffusion coefficient for actin
m0      = 2.0;                  % competition for actin monomers
K       = 1.0;  

Tend    = 25.0;
dt      = 0.1; %0.05
Nt      = Tend/dt; % number of time steps
Na      = 101; % number of space steps
dxa     = 5.0/((Na-1)/2);
L       = 10.0;
Xa      = 0:dxa:L;
Xb      = 0:dxa:L;

Da = 0.5; % diffusion coefficient
pa = dt*Da/(dxa^2);

Ka = 1; % branched coeff
Kb = 1; % bundled coeff
alpha = 2;

F = @(U,V) -m0*U.*V;

a1 = 0.1 + 0.9.*rand(Na,1);
b1 = 0.1 + 0.9.*rand(Na,1);

vid = 0;
vidObj = VideoWriter('PDESolver1','MPEG-4');

for i=1:Nt

    % Laplacian difference operator with no flux boundary conditions
    % Crank-Nicolson operators
    II1 = speye(Na,Na);
    Lapdiff1 = spdiags(ones(Na,1)*[0 1 -2 1 0], [-Na+1, -1, 0 , 1, Na-1], Na, Na);
    Lapdiff1(1,1) = -2; Lapdiff1(1,2) = 1; Lapdiff1(1,Na) = 1;
    Lapdiff1(Na,1) = 1; Lapdiff1(Na,Na-1) = 1; Lapdiff1(Na,Na) = -2;
    Hm1 = II1+(pa/2)*Lapdiff1;
    Hs1 = II1-(pa/2)*Lapdiff1;


    diffRHSa1 = Hm1*a1;
    diffRHSb1 = Hm1*b1;


    rxna1 = dt*(F(a1,b1) + K*(a1 - a1.*a1));
    rxnb1 = dt*(F(b1,a1) + K*(b1 - b1.*b1));


    a1 = Hs1\(diffRHSa1+rxna1);
    b1 = Hs1\(diffRHSb1+rxnb1);

    figure(1)
    plot(Xa,a1,'-o','markerfacecolor',[159 219 229]/255,'linewidth',3); hold on;
    plot(Xa,b1,'-ok','markerfacecolor','k','linewidth',3);
    set(gca,'fontname','times','fontsize',20); box on;
    set(gcf,'color','w');
    legend('branched','bundled'); 
    ylim([-1 1])
    hold off

    if vid==1
        currframe = getframe(gcf);
        writeVideo(vidObj,currframe);
    end

end
