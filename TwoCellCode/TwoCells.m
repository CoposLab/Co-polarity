% Simulate competition for resources for actin filaments
% and fully coupled to stochastic biochemistry for polarity proteins
% on a circular cell membrane (periodic BC, 1.5D)
%
% Mechanics -- two actin networks: branched (a) and contractile (b)
% Polarity proteins : Rac (X) and Rho (Y)
%
% Last updated: 1/31/2020
% Calina Copos

clear;
close all;
clc;

vid = 0;
% vidObj1 = VideoWriter('lineplot1.mp4','MPEG-4');
% vidObjCol1 = VideoWriter('colorplot1.mp4','MPEG-4');
% vidObjRR1 = VideoWriter('colorplotrr1.mp4','MPEG-4');
% vidObj2 = VideoWriter('lineplot2.mp4','MPEG-4');
% vidObjCol2 = VideoWriter('colorplot2.mp4','MPEG-4');
% vidObjRR2 = VideoWriter('colorplotrr2.mp4','MPEG-4');

counter_ppp = 1;
ppp = 1;

while (ppp<=1)
    counter_ppp = counter_ppp+1;

    clearvars -except counter_ppp vid vidObj1 ppp vidObjCol1 vidObjRR1 vidObj2 vidObjCol2 vidObjRR2
    % Set actin filament parameters
    %
    Da      = 0.5;                  % diffusion coefficient for actin
    m0      = 2.0;                  % competition for actin monomers

    % Set polarity protein parameters
    %
    N       = 200;                  % total number of molecules in the cell (conserved)
    ron     = 0.001;                % spontaneous association
    rfb     = 1.0;                  % enhanced association
    roff    = 0.9;                  % disaassociation
    D       = Da;                   % diffusion coefficient for membrane-bound particles

    % Set feedback (or coupling) strength
    %
    alpha = 2;
    beta = 2;

    % Set discretization
    %
    L      = 10.0;                  % cell length
    dt     = 0.01;                  % temporal discretization
    Na     = 101;                   % number of space steps
    dxa    = 5.0/((Na-1)/2);        % spatial discretization
    Xa     = 0:dxa:L;
    pa     = dt*Da/(dxa^2);
    Tend   = 5.0;                  % total simulation time
    Nt     = Tend/dt;
    dx     = sqrt(2*D*dt);
    tplot  = 100;

    posx1 = zeros(N,Nt);              % array of positions of X(t) cell 1
    posy1 = zeros(N,Nt);              % array of positions of Y(t) cell 1

    posx2 = zeros(N,Nt);              % array of positions of X(t) cell 2
    posy2 = zeros(N,Nt);              % array of positions of Y(t) cell 2

    % Actin reaction term
    %
    F = @(U,V) -U.*U - m0*U.*V;

    % Set initial conditions for actin distribution
    %
    ictype = 2;
    %    1 = step in the middle
    %    2 = random
    %    3 = sigmoidal curve
    %    4 = odd condition #1: branched peak in the middle, contractile peaks at
    %    the front and rear cell
    %    5 = odd condition #2: both peaks in the middle (w/ noise)

    a1       = zeros(N,1);
    anew1    = zeros(N,1);
    b1       = zeros(N,1);
    bnew1    = zeros(N,1);

    a2       = zeros(N,1);
    anew2    = zeros(N,1);
    b2       = zeros(N,1);
    bnew2    = zeros(N,1);

    % (1) pulse in middle
    if (ictype==1)
        a1   = ones(N,1);
        anew1= ones(N,1);

        b1   = 0.5*ones(N,1);
        bnew1= 0.5*ones(N,1);

        a2   = ones(N,1);
        anew2= ones(N,1);

        b2   = 0.5*ones(N,1);
        bnew2= 0.5*ones(N,1);

        % (2) random
    elseif (ictype==2)
        a1 = 0.1 + 0.9.*rand(length(Xa),1);
        b1 = 0.1 + 0.9.*rand(length(Xa),1);

        a2 = 0.1 + 0.9.*rand(length(Xa),1);
        b2 = 0.1 + 0.9.*rand(length(Xa),1);

        % (3) arctangent
    elseif (ictype==3)
        steepness = 20;
        a1 = (tanh(steepness*(X1-0.375)) - tanh(steepness*(X1-1.125)) + 0.2)/2.2;
        b1 = (2 - tanh(steepness*(X1-0.375)) + tanh(steepness*(X1-1.125)) + 0.2)/2.2;

        a2 = (tanh(steepness*(X2-0.375)) - tanh(steepness*(X2-1.125)) + 0.2)/2.2;
        b2 = (2 - tanh(steepness*(X2-0.375)) + tanh(steepness*(X2-1.125)) + 0.2)/2.2;

        %a = (tanh(steepness*(X-0.5)) - tanh(steepness*(X-1.5)) + 0.2)/2.2;
        %b = (2 - tanh(steepness*(X-0.5)) + tanh(steepness*(X-1.5)) +0.2)/2.2;
    elseif (ictype==4)
        % (4) odd condition #1 (multiple peaks)
        steepness = 10;
        a1 = (1-cos(3*Xa*pi/5))/2; a1=a1';
        b1 = (tanh(steepness*(Xa-7.5))+1)/2 + (1-tanh(steepness*(Xa-2.5)))/2; b1=b1';

        a2 = (1-cos(3*Xa*pi/5))/2; a2=a2';
        b2 = (tanh(steepness*(Xa-7.5))+1)/2 + (1-tanh(steepness*(Xa-2.5)))/2; b2=b2';

    elseif (ictype==5)
        % (4) odd condition #2
        steepness = 10;
        b1 = (1-cos(Xa*pi/5))/2; b1=b1';
        a1 = (tanh(steepness*(Xa-7.5))+1)/2 + (1-tanh(steepness*(Xa-2.5)))/2; a1=a1';

        b2 = (1-cos(Xa*pi/5))/2; b2=b2';
        a2 = (tanh(steepness*(Xa-7.5))+1)/2 + (1-tanh(steepness*(Xa-2.5)))/2; a2=a2';

    elseif (ictype==6)
        % (5) odd condition #3
        mu = 1.8; sigma = 0.1;
        a1 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); a1=a1';
        a2 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); a2=a2';
        mu = 1.9; sigma = 0.1;
        b1 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); b1=b1';
        b2 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); b2=b2';
    end

    % Laplacian difference operator with no flux boundary conditions
    % Crank-Nicolson operators
    II1 = speye(Na,Na);
    Lapdiff1 = spdiags(ones(Na,1)*[0 1 -2 1 0], [-Na+1, -1, 0 , 1, Na-1], Na, Na);
    Lapdiff1(1,1) = -2; Lapdiff1(1,2) = 1; Lapdiff1(1,Na) = 1;
    Lapdiff1(Na,1) = 1; Lapdiff1(Na,Na-1) = 1; Lapdiff1(Na,Na) = -2;
    Hm1 = II1+(pa/2)*Lapdiff1;
    Hs1 = II1-(pa/2)*Lapdiff1;

    II2 = speye(Na,Na);
    Lapdiff2 = spdiags(ones(Na,1)*[0 1 -2 1 0], [-Na+1, -1, 0 , 1, Na-1], Na, Na);
    Lapdiff2(1,1) = -2; Lapdiff2(1,2) = 1; Lapdiff2(1,Na) = 1;
    Lapdiff2(Na,1) = 1; Lapdiff2(Na,Na-1) = 1; Lapdiff2(Na,Na) = -2;
    Hm2 = II2+(pa/2)*Lapdiff2;
    Hs2 = II2-(pa/2)*Lapdiff2;

    % Setup polarity concentrations
    %
    MAX_OUTPUT_LENGTH = 10000;
    nx1   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
    ny1   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
    Tx1   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for X(t)
    Ty1   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for Y(t)
    X1    = zeros(MAX_OUTPUT_LENGTH,1);       % number of X(t) molecules on the membrane
    Y1    = zeros(MAX_OUTPUT_LENGTH,1);       % number of Y(t) molecules on the membrane
    NNx1  = zeros(Nt,1);
    NNy1  = zeros(Nt,1);

    nx2   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
    ny2   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
    Tx2   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for X(t)
    Ty2   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for Y(t)
    X2    = zeros(MAX_OUTPUT_LENGTH,1);       % number of X(t) molecules on the membrane
    Y2    = zeros(MAX_OUTPUT_LENGTH,1);       % number of Y(t) molecules on the membrane
    NNx2  = zeros(Nt,1);
    NNy2  = zeros(Nt,1);

    % Set initial conditions for polarity molecules distribution
    %
    rxn_count_x1       = 1;
    rxn_count_y1       = 1;

    rxn_count_x2       = 1;
    rxn_count_y2       = 1;

    X1(1)              = 0.1*N;                 % # of particles on membrane
    Y1(1)              = 0.1*N;                 % # of particles on membrane
    NNx1(1)            = X1(1);
    NNy1(1)            = Y1(1);
    Tx1(1)             = 0.0;
    Ty1(1)             = 0.0;
    nx1(1:X1(1),1)      = 1;                     % activate mem-bound particles
    ny1(1:X1(1),1)      = 1;
    r1 = randperm(ceil(L/(0.0102)),X1(1)+Y1(1))*0.0102;
    posx1(1:X1(1),1)=r1(1:X1(1));
    posy1(1:Y1(1),1)=r1(X1(1)+1:end);

    X2(1)              = 0.1*N;                 % # of particles on membrane
    Y2(1)              = 0.1*N;                 % # of particles on membrane
    NNx2(1)            = X2(1);
    NNy2(1)            = Y2(1);
    Tx2(1)             = 0.0;
    Ty2(1)             = 0.0;
    nx2(1:X2(1),1)      = 1;                     % activate mem-bound particles
    ny2(1:X2(1),1)      = 1;
    r2 = randperm(ceil(L/(0.0102)),X2(1)+Y2(1))*0.0102;
    posx2(1:X2(1),1)=r2(1:X2(1));
    posy2(1:Y2(1),1)=r2(X2(1)+1:end);

    % Sample concentration at actin filament spatial scale
    %
    [s1,xC1,yC1] = resamplePolarityMolecules(posx1(1:NNx1(1),1),posy1(1:NNy1(1),1),NNx1(1),NNy1(1),L,Na);
    [s2,xC2,yC2] = resamplePolarityMolecules(posx2(1:NNx2(1),1),posy2(1:NNy2(1),1),NNx2(1),NNy2(1),L,Na);

    aic1 = a1;
    bic1 = b1;

    aic2 = a2;
    bic2 = b2;

    % Setup convergence checks for actin quantities
    %
    conv1 = zeros(Nt,2);
    conv2 = zeros(Nt,2);
    convsup = zeros(Nt,2);

    % Set movie making
    %
    if vid==1
        vidObj1.FrameRate = 5;
        vidObj1.Quality = 75;
        open(vidObj1);

        vidObjCol1.FrameRate = 5;
        vidObjCol1.Quality = 75;
        open(vidObjCol1);

        vidObjRR1.FrameRate = 5;
        vidObjRR1.Quality = 75;
        open(vidObjRR1);

        vidObj2.FrameRate = 5;
        vidObj2.Quality = 75;
        open(vidObj2);

        vidObjCol2.FrameRate = 5;
        vidObjCol2.Quality = 75;
        open(vidObjCol2);

        vidObjRR2.FrameRate = 5;
        vidObjRR2.Quality = 75;
        open(vidObjRR2);
    end

    % Plot the initial condition
    figure(ppp);
    subplot(1,2,1); %Cell 1
    plot(Xa,a1,'-o','markerfacecolor',[159 219 229]/255,'linewidth',3); hold on;
    plot(Xa,b1,'-ok','markerfacecolor','k','linewidth',3);
    plot(s1,xC1,'-.','color',[0 0.45 0.75],'linewidth',1);
    plot(s1,yC1,'-.k','linewidth',1);
    % xlim([0 10]); ylim([0 2]);
    %title('Time = 0');
    set(gca,'fontname','times','fontsize',20); box on;
    lgd = legend('Branched network','Bundled network','Rac','Rho','Location','northeast');
    lgd.NumColumns = 2;
    set(gcf,'color','w');
    title('Cell 1')
    hold off;
    %keyboard
    pause(1.0);

    if vid==1
        currFrame = getframe(gcf);
        writeVideo(vidObj1,currFrame);
    end

    subplot(1,2,2); %Cell 2
    plot(Xa,a2,'-o','markerfacecolor',[159 219 229]/255,'linewidth',3); hold on;
    plot(Xa,b2,'-ok','markerfacecolor','k','linewidth',3);
    plot(s2,xC2,'-.','color',[0 0.45 0.75],'linewidth',1);
    plot(s2,yC2,'-.k','linewidth',1);
    % xlim([0 10]); ylim([0 2]);
    %title('Time = 0');
    set(gca,'fontname','times','fontsize',20); box on;
    lgd = legend('Branched network','Bundled network','Rac','Rho','Location','northeast');
    lgd.NumColumns = 2;
    set(gcf,'color','w');
    title('Cell 2')
    hold off;
    %keyboard
    pause(1.0);

    if vid==1
        currFrame = getframe(gcf);
        writeVideo(vidObj2,currFrame);
    end

    %Plot on circle
    %Define colors
    colorLength = 50;
    red = [1,0,0];
    blue = [0,0.5,1];
    myColors = [linspace(red(1),blue(1),colorLength)',linspace(red(2),blue(2),colorLength)',linspace(red(3),blue(3),colorLength)'];

    [th,rad] = meshgrid((0:3.6:360)*pi/180,0.96:0.01:1);
    [Xcol,Ycol] = pol2cart(th,rad);
    ZBranch1 = [a1 a1 a1 a1 a1]';
    ZBund1 = [b1 b1 b1 b1 b1]';

    ZBranch2 = [a2 a2 a2 a2 a2]';
    ZBund2 = [b2 b2 b2 b2 b2]';

    % figure(11); %Branched
    % % subplot(1,2,1);
    % contourf(Xcol,Ycol,Z1,100,'LineStyle','none')
    % axis square
    % col=colorbar;
    % colormap(myColors);
    % clim([0,max(max(a),max(b))])
    % title('Branched Actin')
    % ylabel(col,'Concentration')
    % set(gca,'XTick',[], 'YTick', [])
    % 
    % figure(12); %Bundled
    % % subplot(1,2,2);
    % contourf(Xcol,Ycol,Z2,100,'LineStyle','none')
    % axis square
    % col=colorbar;
    % colormap(flip(myColors));
    % clim([0,max(max(a),max(b))])
    % title('Bundled Actin')
    % ylabel(col,'Concentration')
    % set(gca,'XTick',[], 'YTick', [])

    %Plot both circles together
    figure(13);
    colFrame = subplot(1,2,1);
    contourf(Xcol,Ycol,ZBranch1-ZBund1,100,'LineStyle','none')
    axis square
    col=colorbar;
    colormap(myColors);
    clim([-max(max(a1),max(b1)),max(max(a1),max(b1))])
    title('Cell 1 Combined: Blue=Branched, Red=Bundled')
    ylabel(col,'Concentration Branched - Concentration Bundled')
    set(gca,'XTick',[], 'YTick', [])

    [Ma,Ia] = max(a1);
    [Mb,Ib] = max(b1);

    hold on;
    plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
    plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
    plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
    hold off;

    if vid==1
        currFrame = getframe(colFrame);
        writeVideo(vidObjCol1,currFrame);
    end

    colFrame = subplot(1,2,2);
    contourf(Xcol,Ycol,ZBranch2-ZBund2,100,'LineStyle','none')
    axis square
    col=colorbar;
    colormap(myColors);
    clim([-max(max(a2),max(b2)),max(max(a2),max(b2))])
    title('Cell 2 Combined: Blue=Branched, Red=Bundled')
    ylabel(col,'Concentration Branched - Concentration Bundled')
    set(gca,'XTick',[], 'YTick', [])

    [Ma,Ia] = max(a2);
    [Mb,Ib] = max(b2);

    hold on;
    plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
    plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
    plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
    hold off;

    if vid==1
        currFrame = getframe(colFrame);
        writeVideo(vidObjCol2,currFrame);
    end

    % Plot Rac and Rho
    % figure(14);
    % colRRFrame = subplot(1,2,1);
    % ZRac1 = [xC1 xC1 xC1 xC1 xC1]';
    % ZRho1 = [yC1 yC1 yC1 yC1 yC1]';
    % contourf(Xcol,Ycol,ZRac1-ZRho1,100,'LineStyle','none');
    % axis square
    % col=colorbar;
    % colormap(myColors);
    % clim([-max(max(xC1),max(yC1)),max(max(xC1),max(yC1))])
    % title('Cell 1 Combined: Blue=Rac, Red=Rho')
    % ylabel(col,'Concentration Rac - Concentration Rho')
    % set(gca,'XTick',[], 'YTick', [])
    % 
    % if vid==1
    %     currFrame = getframe(colRRFrame);
    %     writeVideo(vidObjRR1,currFrame);
    % end
    % 
    % colRRFrame = subplot(1,2,2);
    % ZRac2 = [xC2 xC2 xC2 xC2 xC2]';
    % ZRho2 = [yC2 yC2 yC2 yC2 yC2]';
    % contourf(Xcol,Ycol,ZRac2-ZRho2,100,'LineStyle','none');
    % axis square
    % col=colorbar;
    % colormap(myColors);
    % clim([-max(max(xC2),max(yC2)),max(max(xC2),max(yC2))])
    % title('Cell 2 Combined: Blue=Rac, Red=Rho')
    % ylabel(col,'Concentration Rac - Concentration Rho')
    % set(gca,'XTick',[], 'YTick', [])

    if vid==1
        currFrame = getframe(colRRFrame);
        writeVideo(vidObjRR2,currFrame);
    end

    %% Run simulation
    %
    tic
    quit_cond = 0;
    cond = 0;
    for t=1:(Nt-1)

        %% Run biochemisty
        [Konx,Kony,Kfbx,Kfby,Koffx,Koffy] = spatialrates(ron,rfb,roff,a1,b1,s1,beta,cond); % set rates

        %Cell 1
        if((t-1)*dt<Tx1(rxn_count_x1))
            NNx1(t+1) = X1(rxn_count_x1-1);
        else
            nnx = X1(rxn_count_x1);
            taux = zeros(nnx,1);
            dn = zeros(nnx,1);
            r1 = rand(nnx,1);

            if(nnx==0)
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nnx          % all agents
                konx = interp1(s1,Konx,posx1(j,t));
                koffx = interp1(s1,Koffx,posx1(j,t));
                kfbx = interp1(s1,Kfbx,posx1(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffx + (konx+kfbx*nnx/N)*(N/nnx-1);
                taux(j) = -log(r1(j))/a0;
                rr = rand(1,1);
                dn(j) = (rr<((konx+kfbx*nnx/N)*(N/nnx-1)/a0))*1.0 + (rr>=((konx+kfbx*nnx/N)*(N/nnx-1)/a0))*(-1.0);
            end

            [mintaux1,minidx1] = min(taux(1:j));       % find first chemical rxn
            Tx1(rxn_count_x1+1) = Tx1(rxn_count_x1) + mintaux1;
            X1(rxn_count_x1+1) = nnx + dn(minidx1);
            rxn_count_x1 = rxn_count_x1 + 1;
            NNx1(t+1) = X1(rxn_count_x1-1);
        end

        %Cell 2
        if((t-1)*dt<Tx2(rxn_count_x2))
            NNx2(t+1) = X2(rxn_count_x2-1);
        else
            nnx = X2(rxn_count_x2);
            taux = zeros(nnx,1);
            dn = zeros(nnx,1);
            r2 = rand(nnx,1);

            if(nnx==0)
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nnx          % all agents
                konx = interp1(s2,Konx,posx2(j,t));
                koffx = interp1(s2,Koffx,posx2(j,t));
                kfbx = interp1(s2,Kfbx,posx2(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffx + (konx+kfbx*nnx/N)*(N/nnx-1);
                taux(j) = -log(r2(j))/a0;
                rr = rand(1,1);
                dn(j) = (rr<((konx+kfbx*nnx/N)*(N/nnx-1)/a0))*1.0 + (rr>=((konx+kfbx*nnx/N)*(N/nnx-1)/a0))*(-1.0);
            end

            [mintaux2,minidx2] = min(taux(1:j));       % find first chemical rxn
            Tx2(rxn_count_x2+1) = Tx2(rxn_count_x2) + mintaux2;
            X2(rxn_count_x2+1) = nnx + dn(minidx2);
            rxn_count_x2 = rxn_count_x2 + 1;
            NNx2(t+1) = X2(rxn_count_x2-1);
        end

        %Cell 1
        if((t-1)*dt<Ty1(rxn_count_y1))
            NNy1(t+1) = Y1(rxn_count_y1-1);
        else
            nny = Y1(rxn_count_y1);
            tauy = zeros(nny,1);
            dn = zeros(nny,1);
            r1 = rand(nny,1);

            if(nny==0)
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nny          % all agents
                kony = interp1(s1,Kony,posy1(j,t));
                koffy = interp1(s1,Koffy,posy1(j,t));
                kfby = interp1(s1,Kfby,posy1(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffy + (kony+kfby*nny/N)*(N/nny-1);
                tauy(j) = -log(r1(j))/a0;
                rr = rand(1,1);
                dn(j) = (rr<((kony+kfby*nny/N)*(N/nny-1)/a0))*1.0 + (rr>=((kony+kfby*nny/N)*(N/nny-1)/a0))*(-1.0);
            end

            [mintauy1,minidy1] = min(tauy(1:j));       % find first chemical rxn
            Ty1(rxn_count_y1+1) = Ty1(rxn_count_y1) + mintauy1;
            Y1(rxn_count_y1+1) = nny + dn(minidy1);
            rxn_count_y1 = rxn_count_y1 + 1;
            NNy1(t+1) = Y1(rxn_count_y1-1);
        end

        %Cell 2
        if((t-1)*dt<Ty2(rxn_count_y2))
            NNy2(t+1) = Y2(rxn_count_y2-1);
        else
            nny = Y2(rxn_count_y2);
            tauy = zeros(nny,1);
            dn = zeros(nny,1);
            r2 = rand(nny,1);

            if(nny==0)
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nny          % all agents
                kony = interp1(s2,Kony,posy2(j,t));
                koffy = interp1(s2,Koffy,posy2(j,t));
                kfby = interp1(s2,Kfby,posy2(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffy + (kony+kfby*nny/N)*(N/nny-1);
                tauy(j) = -log(r2(j))/a0;
                rr = rand(1,1);
                dn(j) = (rr<((kony+kfby*nny/N)*(N/nny-1)/a0))*1.0 + (rr>=((kony+kfby*nny/N)*(N/nny-1)/a0))*(-1.0);
            end

            [mintauy2,minidy2] = min(tauy(1:j));       % find first chemical rxn
            Ty2(rxn_count_y2+1) = Ty2(rxn_count_y2) + mintauy2;
            Y2(rxn_count_y2+1) = nny + dn(minidy2);
            rxn_count_y2 = rxn_count_y2 + 1;
            NNy2(t+1) = Y2(rxn_count_y2-1);
        end

        if (quit_cond==1)
            sprintf("It's happening at kk = %d, ppp = %d\n",kk,ppp)
            ppp = counter_ppp-1;
            break
        end

        %% Run diffusion of membrane-bound polarity proteins
        p  = 0.5;                  % probability of hoping left or right

        % Fetch the number of particles at this time
        K1_1 = NNx1(t+1);
        K2_1 = NNy1(t+1);

        K1_2 = NNx2(t+1);
        K2_2 = NNy2(t+1);

        % Between reactions, perform Brownian motion with periodic BC
        r1 = rand(K1_1,1);    % coin flip
        nx1(1:K1_1,t+1) = 1;
        posx1(1:K1_1,t+1) = posx1(1:K1_1,t) + dx*((r1<p)*1.0 + (r1>(1-p))*(-1.0));

        r2 = rand(K1_2,1);    % coin flip
        nx2(1:K1_2,t+1) = 1;
        posx2(1:K1_2,t+1) = posx2(1:K1_2,t) + dx*((r2<p)*1.0 + (r2>(1-p))*(-1.0));

        r1 = rand(K2_1,1);    % coin flip
        ny1(1:K2_1,t+1) = 1;
        posy1(1:K2_1,t+1) = posy1(1:K2_1,t) + dx*((r1<p)*1.0 + (r1>(1-p))*(-1.0));

        r2 = rand(K2_2,1);    % coin flip
        ny2(1:K2_2,t+1) = 1;
        posy2(1:K2_2,t+1) = posy2(1:K2_2,t) + dx*((r2<p)*1.0 + (r2>(1-p))*(-1.0));

        % Check for collision(s) and resolve any collisions
        % Resolution strategy: No one advances
        %
        firstcoll = sum(ismembertol(posx1(1:K1_1,t+1),posy1(1:K2_1,t+1),0.005,'DataScale',1));
        if firstcoll~=0
            % Get indices of collisions
            aa = ismembertol(posx1(1:K1_1,t+1),posy1(1:K2_1,t+1),0.005,'DataScale',1);
            list_idx = find(aa~=0);
            bb = ismembertol(posy1(1:K2_1,t+1),posx1(1:K1_1,t+1),0.005,'DataScale',1);
            list_idy = find(bb~=0);

            posx1(list_idx,t+1) = posx1(list_idx,t);
            posy1(list_idy,t+1) = posy1(list_idy,t);
        end

        firstcoll = sum(ismembertol(posx2(1:K1_2,t+1),posy2(1:K2_2,t+1),0.005,'DataScale',1));
        if firstcoll~=0
            % Get indices of collisions
            aa = ismembertol(posx2(1:K1_2,t+1),posy2(1:K2_2,t+1),0.005,'DataScale',1);
            list_idx = find(aa~=0);
            bb = ismembertol(posy2(1:K2_2,t+1),posx2(1:K1_2,t+1),0.005,'DataScale',1);
            list_idy = find(bb~=0);

            posx2(list_idx,t+1) = posx2(list_idx,t);
            posy2(list_idy,t+1) = posy2(list_idy,t);
        end

        % Enforce periodic boundary conditions
        posx1(1:K1_1,t+1) = posx1(1:K1_1,t+1) + (-L).*(posx1(1:K1_1,t+1)>L) + (L).*(posx1(1:K1_1,t+1)<0.0);
        posy1(1:K2_1,t+1) = posy1(1:K2_1,t+1) + (-L).*(posy1(1:K2_1,t+1)>L) + (L).*(posy1(1:K2_1,t+1)<0.0);

        posx2(1:K1_2,t+1) = posx2(1:K1_2,t+1) + (-L).*(posx2(1:K1_2,t+1)>L) + (L).*(posx2(1:K1_2,t+1)<0.0);
        posy1(1:K2_2,t+1) = posy2(1:K2_2,t+1) + (-L).*(posy2(1:K2_2,t+1)>L) + (L).*(posy2(1:K2_2,t+1)<0.0);

        % Enforce no-flux boundary conditions
        %posx(1:K1,t+1) = posx(1:K1,t+1) + (posx(1:K1,t)-posx(1:K1,t+1)).*(posx(1:K1,t+1)>L) + (posx(1:K1,t)-posx(1:K1,t+1)).*(posx(1:K1,t+1)<0.0);
        %posy(1:K2,t+1) = posy(1:K2,t+1) + (posy(1:K2,t)-posy(1:K2,t+1)).*(posy(1:K2,t+1)>L) + (posy(1:K2,t)-posy(1:K2,t+1)).*(posy(1:K2,t+1)<0.0);

        %% Determine if a biochemical rxn has occured - update positions

        % Find spontaneous association location cell 1
        ss1 = sort(posx1(1:K1_1,t));
        [ijk1] = find(ss1==posx1(minidx1,t),1);
        prevind1 = (ijk1-1)*(ijk1>1) + (K1_1)*(ijk1==1);
        nextind1 = (ijk1+1)*(ijk1<K1_1) + 1*(ijk1==K1_1);
        x2 = posx1(minidx1,t)+(ss1(prevind1)-posx1(minidx1,t))/2;
        x1 = posx1(minidx1,t)+(ss1(nextind1)-posx1(minidx1,t))/2;
        locx = (x2-x1).*rand(1,1) + x1; % random location halfway between the closest left/right particles
        ss1 = sort(posy1(1:K2_1,t));
        [ijk1] = find(ss1==posy1(minidy1,t),1);
        prevind1 = (ijk1-1)*(ijk1>1) + (K2_1)*(ijk1==1);
        nextind1 = (ijk1+1)*(ijk1<K2_1) + 1*(ijk1==K2_1);
        y2 = posy1(minidy1,t)+(ss1(prevind1)-posy1(minidy1,t))/2;
        y1 = posy1(minidy1,t)+(ss1(nextind1)-posy1(minidy1,t))/2;
        locy = (y2-y1).*rand(1,1) + y1; % random location halfway between the closest left/right particles

        ponx1 = ron/(ron+rfb*(N-K1_1));
        pony1 = ron/(ron+rfb*(N-K2_1));

        % Find spontaneous association location cell 2
        ss2 = sort(posx2(1:K1_2,t));
        [ijk2] = find(ss2==posx2(minidx2,t),1);
        prevind2 = (ijk2-1)*(ijk2>1) + (K1_2)*(ijk2==1);
        nextind2 = (ijk2+1)*(ijk2<K1_2) + 1*(ijk2==K1_2);
        x2 = posx2(minidx2,t)+(ss2(prevind2)-posx2(minidx2,t))/2;
        x1 = posx2(minidx2,t)+(ss2(nextind2)-posx2(minidx2,t))/2;
        locx = (x2-x1).*rand(1,1) + x1; % random location halfway between the closest left/right particles
        ss2 = sort(posy2(1:K2_2,t));
        [ijk2] = find(ss2==posy2(minidy2,t),1);
        prevind2 = (ijk2-1)*(ijk2>1) + (K2_2)*(ijk2==1);
        nextind2 = (ijk2+1)*(ijk2<K2_2) + 1*(ijk2==K2_2);
        y2 = posy2(minidy2,t)+(ss2(prevind2)-posy2(minidy2,t))/2;
        y1 = posy2(minidy2,t)+(ss2(nextind2)-posy2(minidy2,t))/2;
        locy = (y2-y1).*rand(1,1) + y1; % random location halfway between the closest left/right particles

        ponx2 = ron/(ron+rfb*(N-K1_2));
        pony2 = ron/(ron+rfb*(N-K2_2));

        %Cell 1
        if(NNx1(t+1) < NNx1(t))                % diassociation event (particle off)
            oldcol = posx1(minidx1,1:end);
            othercols = posx1([1:minidx1-1,minidx1+1:K1_1],1:end);
            otherothercols = posx1(K1_1+1:end,1:end);
            newpos = [othercols;oldcol;otherothercols];
            posx1 = newpos;
            nx1(K1_1,t+1) = 0;
        elseif(NNx1(t+1) > NNx1(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posx1(K1_1,t+1) = posx1(K1_1,t)+(rr<ponx1)*locx;              % on event
            posx1(K1_1,t+1) = posx1(K1_1,t)+(rr>=ponx1)*posx1(minidx1,t);   % recruitment event
            nx1(K1_1+1,t+1) = 1;
        end

        %Cell 2
        if(NNx2(t+1) < NNx2(t))                % diassociation event (particle off)
            oldcol = posx2(minidx1,1:end);
            othercols = posx2([1:minidx1-1,minidx1+1:K1_2],1:end);
            otherothercols = posx2(K1_2+1:end,1:end);
            newpos = [othercols;oldcol;otherothercols];
            posx2 = newpos;
            nx2(K1_2,t+1) = 0;
        elseif(NNx2(t+1) > NNx2(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posx2(K1_2,t+1) = posx2(K1_2,t)+(rr<ponx2)*locx;              % on event
            posx2(K1_2,t+1) = posx2(K1_2,t)+(rr>=ponx2)*posx2(minidx1,t);   % recruitment event
            nx2(K1_2+1,t+1) = 1;
        end

        %Cell 1
        if (NNy1(t+1) < NNy1(t))                % diassociation event (particle off)
            oldcol = posy1(minidy1,1:end);
            othercols = posy1([1:minidy1-1,minidy1+1:K2_1],1:end);
            otherothercols = posy1(K2_1+1:end,1:end);
            newpos = [othercols;oldcol;otherothercols];
            posy1 = newpos;
            ny1(K2_1,t+1) = 0;
        elseif(NNy1(t+1) > NNy1(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posy1(K2_1,t+1) = posy1(K2_1,t)+(rr<pony1)*locy;               % on event
            posy1(K2_1,t+1) = posy1(K2_1,t)+(rr>=pony1)*posy1(minidy1,t);    % recruitment event
            ny1(K2_1,t+1) = 1;
        end

        %Cell 2
        if (NNy2(t+1) < NNy2(t))                % diassociation event (particle off)
            oldcol = posy2(minidy1,1:end);
            othercols = posy2([1:minidy1-1,minidy1+1:K2_2],1:end);
            otherothercols = posy2(K2_2+1:end,1:end);
            newpos = [othercols;oldcol;otherothercols];
            posy2 = newpos;
            ny2(K2_2,t+1) = 0;
        elseif(NNy2(t+1) > NNy2(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posy2(K2_2,t+1) = posy2(K2_2,t)+(rr<pony2)*locy;               % on event
            posy2(K2_2,t+1) = posy2(K2_2,t)+(rr>=pony2)*posy2(minidy1,t);    % recruitment event
            ny2(K2_2,t+1) = 1;
        end

        [s1,xC1,yC1] = resamplePolarityMolecules(posx1(1:K1_1,t+1),posy1(1:K2_1,t+1),K1_1,K2_1,L,Na);

        [s2,xC2,yC2] = resamplePolarityMolecules(posx2(1:K1_2,t+1),posy2(1:K2_2,t+1),K1_2,K2_2,L,Na);

        %% Update actin filaments
        diffRHSa1 = Hm1*a1;
        diffRHSb1 = Hm1*b1;

        diffRHSa2 = Hm1*a2;
        diffRHSb2 = Hm1*b2;

        rxna1 = dt*( F(a1,b1) + a1.*(1+alpha*xC1));
        rxnb1 = dt*( F(b1,a1) + b1.*(1+alpha*yC1));

        rxna2 = dt*( F(a2,b2) + a2.*(1+alpha*xC2));
        rxnb2 = dt*( F(b2,a2) + b2.*(1+alpha*yC2));

        a1 = Hs1\(diffRHSa1+rxna1);
        b1 = Hs1\(diffRHSb1+rxnb1);

        a2 = Hs1\(diffRHSa2+rxna2);
        b2 = Hs1\(diffRHSb2+rxnb2);

        %% Plot the solution(s)
        % if mod(t,tplot) == 0
        if t==(Nt-1)
            figure(ppp);
            subplot(1,2,1); %Cell 1
            plot(Xa,a1,'-o','markerfacecolor',[159 219 229]/255,'linewidth',3); hold on;
            plot(Xa,b1,'-ok','markerfacecolor','k','linewidth',3);
            plot(s1,xC1,'-.','color',[0 0.45 0.75],'linewidth',1);
            plot(s1,yC1,'-.k','linewidth',1);
            % xlim([0 10]); ylim([0 2]);
            %title('Time = 0');
            set(gca,'fontname','times','fontsize',20); box on;
            lgd = legend('Branched network','Bundled network','Rac','Rho','Location','northeast');
            lgd.NumColumns = 2;
            set(gcf,'color','w');
            title('Cell 1')
            hold off;
            %keyboard
            pause(1.0);

            if vid==1
                currFrame = getframe(gcf);
                writeVideo(vidObj1,currFrame);
            end

            subplot(1,2,2); %Cell 2
            plot(Xa,a2,'-o','markerfacecolor',[159 219 229]/255,'linewidth',3); hold on;
            plot(Xa,b2,'-ok','markerfacecolor','k','linewidth',3);
            plot(s2,xC2,'-.','color',[0 0.45 0.75],'linewidth',1);
            plot(s2,yC2,'-.k','linewidth',1);
            % xlim([0 10]); ylim([0 2]);
            %title('Time = 0');
            set(gca,'fontname','times','fontsize',20); box on;
            lgd = legend('Branched network','Bundled network','Rac','Rho','Location','northeast');
            lgd.NumColumns = 2;
            set(gcf,'color','w');
            title('Cell 2')
            hold off;
            %keyboard
            pause(1.0);

            if vid==1
                currFrame = getframe(gcf);
                writeVideo(vidObj2,currFrame);
            end

            %Plot on circle
            %Define colors
            colorLength = 50;
            red = [1,0,0];
            blue = [0,0.5,1];
            myColors = [linspace(red(1),blue(1),colorLength)',linspace(red(2),blue(2),colorLength)',linspace(red(3),blue(3),colorLength)'];

            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.96:0.01:1);
            [Xcol,Ycol] = pol2cart(th,rad);
            ZBranch1 = [a1 a1 a1 a1 a1]';
            ZBund1 = [b1 b1 b1 b1 b1]';

            ZBranch2 = [a2 a2 a2 a2 a2]';
            ZBund2 = [b2 b2 b2 b2 b2]';

            % figure(11); %Branched
            % % subplot(1,2,1);
            % contourf(Xcol,Ycol,Z1,100,'LineStyle','none')
            % axis square
            % col=colorbar;
            % colormap(myColors);
            % clim([0,max(max(a),max(b))])
            % title('Branched Actin')
            % ylabel(col,'Concentration')
            % set(gca,'XTick',[], 'YTick', [])
            %
            % figure(12); %Bundled
            % % subplot(1,2,2);
            % contourf(Xcol,Ycol,Z2,100,'LineStyle','none')
            % axis square
            % col=colorbar;
            % colormap(flip(myColors));
            % clim([0,max(max(a),max(b))])
            % title('Bundled Actin')
            % ylabel(col,'Concentration')
            % set(gca,'XTick',[], 'YTick', [])

            %Plot both circles together
            figure(13);
            colFrame = subplot(1,2,1);
            contourf(Xcol,Ycol,ZBranch1-ZBund1,100,'LineStyle','none')
            axis square
            col=colorbar;
            colormap(myColors);
            clim([-max(max(a1),max(b1)),max(max(a1),max(b1))])
            title('Cell 1 Combined: Blue=Branched, Red=Bundled')
            ylabel(col,'Concentration Branched - Concentration Bundled')
            set(gca,'XTick',[], 'YTick', [])

            [Ma,Ia] = max(a1);
            [Mb,Ib] = max(b1);

            hold on;
            plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
            plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
            plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
            hold off;

            if vid==1
                currFrame = getframe(colFrame);
                writeVideo(vidObjCol1,currFrame);
            end

            colFrame = subplot(1,2,2);
            contourf(Xcol,Ycol,ZBranch2-ZBund2,100,'LineStyle','none')
            axis square
            col=colorbar;
            colormap(myColors);
            clim([-max(max(a2),max(b2)),max(max(a2),max(b2))])
            title('Cell 2 Combined: Blue=Branched, Red=Bundled')
            ylabel(col,'Concentration Branched - Concentration Bundled')
            set(gca,'XTick',[], 'YTick', [])

            [Ma,Ia] = max(a2);
            [Mb,Ib] = max(b2);

            hold on;
            plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
            plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
            plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
            hold off;

            if vid==1
                currFrame = getframe(colFrame);
                writeVideo(vidObjCol2,currFrame);
            end

            % Plot Rac and Rho
            % figure(14);
            % colRRFrame = subplot(1,2,1);
            % ZRac1 = [xC1 xC1 xC1 xC1 xC1]';
            % ZRho1 = [yC1 yC1 yC1 yC1 yC1]';
            % contourf(Xcol,Ycol,ZRac1-ZRho1,100,'LineStyle','none');
            % axis square
            % col=colorbar;
            % colormap(myColors);
            % clim([-max(max(xC1),max(yC1)),max(max(xC1),max(yC1))])
            % title('Cell 1 Combined: Blue=Rac, Red=Rho')
            % ylabel(col,'Concentration Rac - Concentration Rho')
            % set(gca,'XTick',[], 'YTick', [])
            % 
            % if vid==1
            %     currFrame = getframe(colRRFrame);
            %     writeVideo(vidObjRR1,currFrame);
            % end
            % 
            % colRRFrame = subplot(1,2,2);
            % ZRac2 = [xC2 xC2 xC2 xC2 xC2]';
            % ZRho2 = [yC2 yC2 yC2 yC2 yC2]';
            % contourf(Xcol,Ycol,ZRac2-ZRho2,100,'LineStyle','none');
            % axis square
            % col=colorbar;
            % colormap(myColors);
            % clim([-max(max(xC2),max(yC2)),max(max(xC2),max(yC2))])
            % title('Cell 2 Combined: Blue=Rac, Red=Rho')
            % ylabel(col,'Concentration Rac - Concentration Rho')
            % set(gca,'XTick',[], 'YTick', [])

            if vid==1
                currFrame = getframe(colRRFrame);
                writeVideo(vidObjRR2,currFrame);
            end
        end
    end

    % measure of polarized state (1 if polarized and 0 otherwise)
    %st = 1*( (abs(a(1)-b(end))<1e-3 || abs(a(end)-b(1))<1e-3 ) && abs(a(1)-a(end))>1e-3 && abs(b(1)-b(end))>1e-3 );

    if vid==1
        close(vidObj1);
        close(vidObjCol1);
        close(vidObjRR1);

        close(vidObj2);
        close(vidObjCol2);
        close(vidObjRR2);
    end
    sprintf('Simulation %d done',ppp)
    toc
    if(quit_cond==0)
        ppp = ppp + 1;
    end
end

%% Plot all particle trajectories
% ccx = [0 0 255]/256.*ones(Nt,1);     % blue
% ccy = [255 219 88]/256.*ones(Nt,1);  % mustard yellow
% time = linspace(0,Tend,Nt);
% figure(10);
% subplot(1,2,1);
% for j=1:2:max(max([NNx1,NNy1]))
%     hold on;
%     %plot(linspace(0,Tend,Nt),pos(j,:));
%     scatter(linspace(0,Tend,Nt),posx1(j,:),1,ccx);
%     scatter(linspace(0,Tend,Nt),posy1(j,:),1,ccy);
%     box on;
%     set(gca,'Color','k','fontsize',20,'fontname','times');
%     pbaspect([3 1 1]);
%     set(gcf,'color','w');
%     title('Cell 1')
%     %xlabel('Time');
%     %ylabel('Location on cell membrane');
%     %yticks([0 0.5 1 1.5 2]);
%     %yticklabels({'0','0.25','0.5','0.75','1'});
% end
% 
% subplot(1,2,2);
% for j=1:2:max(max([NNx2,NNy2]))
%     hold on;
%     %plot(linspace(0,Tend,Nt),pos(j,:));
%     scatter(linspace(0,Tend,Nt),posx2(j,:),1,ccx);
%     scatter(linspace(0,Tend,Nt),posy2(j,:),1,ccy);
%     box on;
%     set(gca,'Color','k','fontsize',20,'fontname','times');
%     pbaspect([3 1 1]);
%     set(gcf,'color','w');
%     title('Cell 2')
%     %xlabel('Time');
%     %ylabel('Location on cell membrane');
%     %yticks([0 0.5 1 1.5 2]);
%     %yticklabels({'0','0.25','0.5','0.75','1'});
% end