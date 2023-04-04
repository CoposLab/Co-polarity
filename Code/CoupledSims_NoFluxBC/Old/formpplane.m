function [] = formpplane(mf,a12hat,a21)

% solution space
M   = 20;
D   = 0.1;  
a21 = 2.0;
af  = linspace(0,1,M);
bf  = linspace(0,1,M);
ar  = linspace(0,1,M);
br  = linspace(0,1,M);

switch ode
    case '1'
    %% no diffusion / no competition
    afprime = af - af.^2 - mf*af.*bf;
    arprime = ar - ar.^2 - a21*ar.*br;
    bfprime = bf - bf.^2 - mf*af.*bf;
    brprime = br - br.^2 - a21*ar.*br;

    ff = @(t,YF) [YF(1)-YF(1)^2-mf*YF(1)*YF(2);YF(2)-YF(2)^2-mf*YF(1)*YF(2)];
    fr = @(t,YR) [YR(1)-YR(1)^2-a21*YR(1)*YR(2);YR(2)-YR(2)^2-a21*YR(1)*YR(2)];

    case '2'
    %% with diffusion / no competition
    afprime = af - af.^2 - mf*af.*bf + d*(ar-af);
    arprime = ar - ar.^2 - a21*ar.*br + d*(br-bf);
    bfprime = bf - bf.^2 - mf*af.*bf + d*(af-ar);
    brprime = br - br.^2 - a21*ar.*br + d*(bf-bf);

    ff = @(t,YF,YR) [YF(1)-YF(1)^2-mf*YF(1)*YF(2)+D*(YR(1)-YF(1));YF(2)-YF(2)^2-mf*YF(1)*YF(2)+D*(YR(2)-YF(2))];
    fr = @(t,YF,YR) [YR(1)-YR(1)^2-a21*YR(1)*YR(2)+D*(YF(1)-YR(1));YR(2)-YR(2)^2-a21*YR(1)*YR(2)+D*(YR(2)-YF(2))];
end

figure(10);
subplot(1,2,1);

%% plot phase portrait
[xf,yf] = meshgrid(af,bf);
uf = zeros(size(xf));
vf = zeros(size(xf));
t = 0;
for i = 1:numel(xf)
    Yprime = ff(t,[xf(i); yf(i)]);
    uf(i) = Yprime(1);
    vf(i) = Yprime(2);
end
quiver(xf,yf,uf,vf,'r'); figure(gcf)
xlabel('front branched actin, a_f')
ylabel('front contractile actin, b_f')
axis tight equal; hold on;
xlim([-0.1 1.1]); ylim([-0.1 1.1]);

switch ode
    case '1'
    %% no diffusion / no competition
    % plot nullclines 
    plot(zeros(length(af),1),af,'-k','linewidth',2); % af nullcline
    plot(af,(1-af)/mf,'-k','linewidth',2); % af nullcline
    plot(bf,zeros(length(bf),1),'--b','linewidth',2); % bf nullcline
    %plot((1-bf)/mf,bf,'--b','linewidth',2); % af nullcline
    plot(af,1-af*mf,'--b','linewidth',2);

    % steady states
    scatter(0,0,100,'ok');
    scatter(0,1,100,'ok');
    scatter(1,0,100,'ok');
    scatter(1/(mf+1),1/(mf+1),100,'ok');
end

end

