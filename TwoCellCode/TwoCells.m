% Simulate competition for resources for actin filaments
% and fully coupled to stochastic biochemistry for polarity proteins
% on a circular cell membrane (periodic BC, 1.5D)
%
% Mechanics -- two actin networks: branched (a) and contractile (b)
% Polarity proteins : Rac (X) and Rho (Y)
%
% Last updated: 1/31/2020
% Calina Copos
addpath('/freeze_colors')
addpath('../SingleCellCode_Published')

clear;
close all;
clc;

savefigs=0;
if savefigs==1
    setnum='20';
    % filenameC1=strcat('savedgraphs/doubleRhoOnCell1_',setnum);
    % filenameC2=strcat('savedgraphs/doubleRhoOnCell2_',setnum);
    filenameCells=strcat('savedgraphs3/100timesRacOnCells_',setnum);
    filenameScatter=strcat('savedgraphs3/100timesRacOnScatter_',setnum);
end

vid = 0;
vidObj1 = VideoWriter('videos1/uncoupledScatter_1.mp4','MPEG-4');
vidObjCol1 = VideoWriter('videos1/uncoupledColor_1.mp4','MPEG-4');
% vidObjRR1 = VideoWriter('colorplotrr1.mp4','MPEG-4');
% vidObj2 = VideoWriter('lineplot2.mp4','MPEG-4');
% vidObjCol2 = VideoWriter('colorplot2.mp4','MPEG-4');
% vidObjRR2 = VideoWriter('colorplotrr2.mp4','MPEG-4');

counter_ppp = 1;
ppp = 1;

while (ppp<=1)
    counter_ppp = counter_ppp+1;

    clearvars -except counter_ppp vid vidObj1 ppp vidObjCol1 vidObjRR1 vidObj2 vidObjCol2 vidObjRR2 savefigs filenameC1 filenameC2 filenameScatter filenameCells
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
    bper   = 0.25;                  % percent overlap
    dt     = 0.01;                  % temporal discretization
    Na     = 101;                   % number of space steps
    dxa    = 5.0/((Na-1)/2);        % spatial discretization
    Xa     = 0:dxa:L;
    Xb     = 0:dxa:L;
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

        a2 = 0.1 + 0.9.*rand(length(Xb),1);
        b2 = 0.1 + 0.9.*rand(length(Xb),1);

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

        a2 = (1-cos(3*Xb*pi/5))/2; a2=a2';
        b2 = (tanh(steepness*(Xb-7.5))+1)/2 + (1-tanh(steepness*(Xb-2.5)))/2; b2=b2';

    elseif (ictype==5)
        % (4) odd condition #2
        steepness = 10;
        b1 = (1-cos(Xa*pi/5))/2; b1=b1';
        a1 = (tanh(steepness*(Xa-7.5))+1)/2 + (1-tanh(steepness*(Xa-2.5)))/2; a1=a1';

        b2 = (1-cos(Xb*pi/5))/2; b2=b2';
        a2 = (tanh(steepness*(Xb-7.5))+1)/2 + (1-tanh(steepness*(Xb-2.5)))/2; a2=a2';

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

    % Boundary between cells
    blen = floor(bper * L); % length of overlap between cell membranes
    boundC1 = (floor(3*Na/4 - ((Na-1)*bper)/2)):((floor(3*Na/4 - ((Na-1)*bper)/2))+(Na-1)*bper); %boundary region in cell 1
    boundC2 = (floor(Na/4 - ((Na-1)*bper)/2)):((floor(Na/4 - ((Na-1)*bper)/2))+(Na-1)*bper); %boundary region in cell 2
    

    % Setup convergence checks for actin quantities
    %
    conv1_1 = zeros(Nt,2);
    conv2_1 = zeros(Nt,2);
    convsup_1 = zeros(Nt,2);

    conv1_2 = zeros(Nt,2);
    conv2_2 = zeros(Nt,2);
    convsup_2 = zeros(Nt,2);

    % Set movie making
    %
    if vid==1
        vidObj1.FrameRate = 5;
        vidObj1.Quality = 75;
        open(vidObj1);

        vidObjCol1.FrameRate = 5;
        vidObjCol1.Quality = 75;
        open(vidObjCol1);

        % vidObjRR1.FrameRate = 5;
        % vidObjRR1.Quality = 75;
        % open(vidObjRR1);
        % 
        % vidObj2.FrameRate = 5;
        % vidObj2.Quality = 75;
        % open(vidObj2);
        % 
        % vidObjCol2.FrameRate = 5;
        % vidObjCol2.Quality = 75;
        % open(vidObjCol2);
        % 
        % vidObjRR2.FrameRate = 5;
        % vidObjRR2.Quality = 75;
        % open(vidObjRR2);
    end

    % Plot the initial condition
    scatterFrame = figure(ppp);
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

    % if vid==1
    %     currFrame = getframe(gcf);
    %     writeVideo(vidObj1,currFrame);
    % end

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
        currFrame = getframe(scatterFrame);
        writeVideo(vidObj1,currFrame);
    end

    %Plot on circle
    %Define colors
    colorLength = 50;
    white = [1,1,1];
    red = [1,0,0];
    blue = [0,0.5,1];
    whitered = [linspace(white(1),red(1),colorLength)',linspace(white(2),red(2),colorLength)',linspace(white(3),red(3),colorLength)'];
    whiteblue = [linspace(white(1),blue(1),colorLength)',linspace(white(2),blue(2),colorLength)',linspace(white(3),blue(3),colorLength)'];
    myColors = [linspace(red(1),blue(1),colorLength)',linspace(red(2),blue(2),colorLength)',linspace(red(3),blue(3),colorLength)'];
    redblue = abs(whiteblue+whitered)./2;
    redwhiteblue = [flip(whitered); whiteblue];

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
    % colormap(whitered);
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
    % colormap(whiteblue);
    % clim([0,max(max(a),max(b))])
    % title('Bundled Actin')
    % ylabel(col,'Concentration')
    % set(gca,'XTick',[], 'YTick', [])

    %Plot both circles together
    % figure(13);
    % colFrame = subplot(1,2,1);
    % contourf(Xcol,Ycol,ZBranch1-ZBund1,100,'LineStyle','none')
    % axis square
    % col=colorbar;
    % colormap(redwhiteblue);
    % clim([-max(max(a1),max(b1)),max(max(a1),max(b1))])
    % title('Cell 1 Combined: Blue=Branched, Red=Bundled')
    % ylabel(col,'Concentration Branched - Concentration Bundled')
    % set(gca,'XTick',[], 'YTick', [])
    % 
    % [Ma,Ia] = max(a1);
    % [Mb,Ib] = max(b1);

    % hold on;
    % plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
    % plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
    % plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
    % hold off;

    % if vid==1
    %     currFrame = getframe(colFrame);
    %     writeVideo(vidObjCol1,currFrame);
    % end

    % colFrame = subplot(1,2,2);
    % contourf(Xcol,Ycol,ZBranch2-ZBund2,100,'LineStyle','none')
    % axis square
    % col=colorbar;
    % colormap(redwhiteblue);
    % clim([-max(max(a2),max(b2)),max(max(a2),max(b2))])
    % title('Cell 2 Combined: Blue=Branched, Red=Bundled')
    % ylabel(col,'Concentration Branched - Concentration Bundled')
    % set(gca,'XTick',[], 'YTick', [])
    % 
    % [Ma,Ia] = max(a2);
    % [Mb,Ib] = max(b2);

    % hold on;
    % plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
    % plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
    % plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
    % hold off;

    % if vid==1
    %     currFrame = getframe(colFrame);
    %     writeVideo(vidObjCol2,currFrame);
    % end

    % Concentric circles
    [th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
    [Xsm,Ysm] = pol2cart(th,rad);
    [th,rad] = meshgrid((0:3.6:360)*pi/180,0.86:0.01:0.9);
    [Xmid,Ymid] = pol2cart(th,rad);

    % % Cell 1
    % colFrame=figure(16);
    % subplot(2,1,1)
    % surf(Xcol,Ycol,ZBranch1);
    % view(2)
    % colormap(whiteblue)
    % freezeColors;
    % freezeColors(colorbar);
    % clim([0,max(max(a1),max(b1))])
    % shading interp
    % hold on;
    % surf(Xmid,Ymid,ZBund1);
    % colormap(whitered)
    % freezeColors;
    % freezeColors(jicolorbar);
    % clim([0,max(b1)])
    % shading interp
    % grid off
    % axis square
    % set(gca,'XTick',[], 'YTick', [])
    % title('Cell 1 Combined: Blue=Branched, Red=Bundled')
    % hold off;
    % 
    % [Ma1,Ia1] = max(a1);
    % [Mb1,Ib1] = max(b1);
    % p1=[0 0];
    % p2=[Xcol(3,Ia1) Ycol(3,Ia1)];
    % dp=p2-p1;
    % 
    % hold on;
    % quiver(p1(1),p1(2),dp(1),dp(2),0,'color',[0 0 0])
    % % plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
    % % plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
    % % plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
    % hold off;
    % 
    % % Cell 2
    % % figure(17)
    % subplot(2,1,2)
    % surf(Xcol,Ycol,ZBranch2);
    % view(2)
    % colormap(whiteblue)
    % freezeColors;
    % freezeColors(colorbar);
    % clim([0,max(max(a2),max(b2))])
    % shading interp
    % hold on;
    % surf(Xmid,Ymid,ZBund2);
    % colormap(whitered)
    % freezeColors;
    % freezeColors(jicolorbar);
    % clim([0,max(b2)])
    % shading interp
    % grid off
    % axis square
    % set(gca,'XTick',[], 'YTick', [])
    % title('Cell 2 Combined: Blue=Branched, Red=Bundled')
    % hold off;
    % 
    % [Ma2,Ia2] = max(a2);
    % [Mb2,Ib2] = max(b2);
    % p1=[0 0];
    % p2=[Xcol(3,Ia2) Ycol(3,Ia2)];
    % dp=p2-p1;
    % 
    % hold on;
    % quiver(p1(1),p1(2),dp(1),dp(2),0,'color',[0 0 0])
    % % plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
    % % plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
    % % plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
    % hold off;
    % 
    % if vid==1
    %     currframe = getframe(colFrame);
    %     writeVideo(vidObjCol1,currFrame);
    % end

    % Plot Rac and Rho
    % figure(14);
    % colRRFrame = subplot(1,2,1);
    % ZRac1 = [xC1 xC1 xC1 xC1 xC1]';
    % ZRho1 = [yC1 yC1 yC1 yC1 yC1]';
    % contourf(Xcol,Ycol,ZRac1-ZRho1,100,'LineStyle','none');
    % axis square
    % col=colorbar;
    % colormap(redwhiteblue);
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
    % colormap(redwhiteblue);
    % clim([-max(max(xC2),max(yC2)),max(max(xC2),max(yC2))])
    % title('Cell 2 Combined: Blue=Rac, Red=Rho')
    % ylabel(col,'Concentration Rac - Concentration Rho')
    % set(gca,'XTick',[], 'YTick', [])

    % if vid==1
    %     currFrame = getframe(colRRFrame);
    %     writeVideo(vidObjRR2,currFrame);
    % end

    %% Run simulation
    %
    tic
    quit_cond = 0;
    cond = 0;
    for t=1:(Nt-1)

        %% Run biochemisty
        [Konx1,Kony1,Kfbx1,Kfby1,Koffx1,Koffy1] = spatialrates(ron,rfb,roff,a1,b1,s1,beta,cond); % set rates
        [Konx2,Kony2,Kfbx2,Kfby2,Koffx2,Koffy2] = spatialrates(ron,rfb,roff,a2,b2,s2,beta,cond);

        % Set konx and kony
        Konx1(boundC1)=Konx1(boundC1)*100;
        Konx2(boundC2)=Konx2(boundC2)*100;

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
                konx1 = interp1(s1,Konx1,posx1(j,t));
                koffx1 = interp1(s1,Koffx1,posx1(j,t));
                kfbx1 = interp1(s1,Kfbx1,posx1(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffx1 + (konx1+kfbx1*nnx/N)*(N/nnx-1);
                taux(j) = -log(r1(j))/a0;
                rr = rand(1,1);
                dn(j) = (rr<((konx1+kfbx1*nnx/N)*(N/nnx-1)/a0))*1.0 + (rr>=((konx1+kfbx1*nnx/N)*(N/nnx-1)/a0))*(-1.0);
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
                konx2 = interp1(s2,Konx2,posx2(j,t));
                koffx2 = interp1(s2,Koffx2,posx2(j,t));
                kfbx2 = interp1(s2,Kfbx2,posx2(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffx2 + (konx2+kfbx2*nnx/N)*(N/nnx-1);
                taux(j) = -log(r2(j))/a0;
                rr = rand(1,1);
                dn(j) = (rr<((konx2+kfbx1*nnx/N)*(N/nnx-1)/a0))*1.0 + (rr>=((konx2+kfbx2*nnx/N)*(N/nnx-1)/a0))*(-1.0);
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
                kony1 = interp1(s1,Kony1,posy1(j,t));
                koffy1 = interp1(s1,Koffy1,posy1(j,t));
                kfby1 = interp1(s1,Kfby1,posy1(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffy1 + (kony1+kfby1*nny/N)*(N/nny-1);
                tauy(j) = -log(r1(j))/a0;
                rr = rand(1,1);
                dn(j) = (rr<((kony1+kfby1*nny/N)*(N/nny-1)/a0))*1.0 + (rr>=((kony1+kfby1*nny/N)*(N/nny-1)/a0))*(-1.0);
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
                kony2 = interp1(s2,Kony2,posy2(j,t));
                koffy2 = interp1(s2,Koffy2,posy2(j,t));
                kfby2 = interp1(s2,Kfby2,posy2(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffy2 + (kony2+kfby2*nny/N)*(N/nny-1);
                tauy(j) = -log(r2(j))/a0;
                rr = rand(1,1);
                dn(j) = (rr<((kony2+kfby2*nny/N)*(N/nny-1)/a0))*1.0 + (rr>=((kony2+kfby2*nny/N)*(N/nny-1)/a0))*(-1.0);
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
        % Cell 1
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

        % Cell 2
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
        posy2(1:K2_2,t+1) = posy2(1:K2_2,t+1) + (-L).*(posy2(1:K2_2,t+1)>L) + (L).*(posy2(1:K2_2,t+1)<0.0);

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
        locx1 = (x2-x1).*rand(1,1) + x1; % random location halfway between the closest left/right particles
        ss1 = sort(posy1(1:K2_1,t));
        [ijk1] = find(ss1==posy1(minidy1,t),1);
        prevind1 = (ijk1-1)*(ijk1>1) + (K2_1)*(ijk1==1);
        nextind1 = (ijk1+1)*(ijk1<K2_1) + 1*(ijk1==K2_1);
        y2 = posy1(minidy1,t)+(ss1(prevind1)-posy1(minidy1,t))/2;
        y1 = posy1(minidy1,t)+(ss1(nextind1)-posy1(minidy1,t))/2;
        locy1 = (y2-y1).*rand(1,1) + y1; % random location halfway between the closest left/right particles

        ponx1 = ron/(ron+rfb*(N-K1_1));
        pony1 = ron/(ron+rfb*(N-K2_1));

        % Find spontaneous association location cell 2
        ss2 = sort(posx2(1:K1_2,t));
        [ijk2] = find(ss2==posx2(minidx2,t),1);
        prevind2 = (ijk2-1)*(ijk2>1) + (K1_2)*(ijk2==1);
        nextind2 = (ijk2+1)*(ijk2<K1_2) + 1*(ijk2==K1_2);
        x2 = posx2(minidx2,t)+(ss2(prevind2)-posx2(minidx2,t))/2;
        x1 = posx2(minidx2,t)+(ss2(nextind2)-posx2(minidx2,t))/2;
        locx2 = (x2-x1).*rand(1,1) + x1; % random location halfway between the closest left/right particles
        ss2 = sort(posy2(1:K2_2,t));
        [ijk2] = find(ss2==posy2(minidy2,t),1);
        prevind2 = (ijk2-1)*(ijk2>1) + (K2_2)*(ijk2==1);
        nextind2 = (ijk2+1)*(ijk2<K2_2) + 1*(ijk2==K2_2);
        y2 = posy2(minidy2,t)+(ss2(prevind2)-posy2(minidy2,t))/2;
        y1 = posy2(minidy2,t)+(ss2(nextind2)-posy2(minidy2,t))/2;
        locy2 = (y2-y1).*rand(1,1) + y1; % random location halfway between the closest left/right particles

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
            posx1(K1_1,t+1) = posx1(K1_1,t)+(rr<ponx1)*locx1;              % on event
            posx1(K1_1,t+1) = posx1(K1_1,t)+(rr>=ponx1)*posx1(minidx1,t);   % recruitment event
            nx1(K1_1+1,t+1) = 1;
        end

        %Cell 2
        if(NNx2(t+1) < NNx2(t))                % diassociation event (particle off)
            oldcol = posx2(minidx2,1:end);
            othercols = posx2([1:minidx2-1,minidx2+1:K1_2],1:end);
            otherothercols = posx2(K1_2+1:end,1:end);
            newpos = [othercols;oldcol;otherothercols];
            posx2 = newpos;
            nx2(K1_2,t+1) = 0;
        elseif(NNx2(t+1) > NNx2(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posx2(K1_2,t+1) = posx2(K1_2,t)+(rr<ponx2)*locx2;              % on event
            posx2(K1_2,t+1) = posx2(K1_2,t)+(rr>=ponx2)*posx2(minidx2,t);   % recruitment event
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
            posy1(K2_1,t+1) = posy1(K2_1,t)+(rr<pony1)*locy1;               % on event
            posy1(K2_1,t+1) = posy1(K2_1,t)+(rr>=pony1)*posy1(minidy1,t);    % recruitment event
            ny1(K2_1,t+1) = 1;
        end

        %Cell 2
        if (NNy2(t+1) < NNy2(t))                % diassociation event (particle off)
            oldcol = posy2(minidy2,1:end);
            othercols = posy2([1:minidy2-1,minidy2+1:K2_2],1:end);
            otherothercols = posy2(K2_2+1:end,1:end);
            newpos = [othercols;oldcol;otherothercols];
            posy2 = newpos;
            ny2(K2_2,t+1) = 0;
        elseif(NNy2(t+1) > NNy2(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posy2(K2_2,t+1) = posy2(K2_2,t)+(rr<pony2)*locy2;               % on event
            posy2(K2_2,t+1) = posy2(K2_2,t)+(rr>=pony2)*posy2(minidy2,t);    % recruitment event
            ny2(K2_2,t+1) = 1;
        end

        [s1,xC1,yC1] = resamplePolarityMolecules(posx1(1:K1_1,t+1),posy1(1:K2_1,t+1),K1_1,K2_1,L,Na);

        [s2,xC2,yC2] = resamplePolarityMolecules(posx2(1:K1_2,t+1),posy2(1:K2_2,t+1),K1_2,K2_2,L,Na);

        %% Update actin filaments
        diffRHSa1 = Hm1*a1;
        diffRHSb1 = Hm1*b1;

        diffRHSa2 = Hm2*a2;
        diffRHSb2 = Hm2*b2;

        rxna1 = dt*( F(a1,b1) + a1.*(1+alpha*xC1));
        rxnb1 = dt*( F(b1,a1) + b1.*(1+alpha*yC1));

        rxna2 = dt*( F(a2,b2) + a2.*(1+alpha*xC2));
        rxnb2 = dt*( F(b2,a2) + b2.*(1+alpha*yC2));

        a1 = Hs1\(diffRHSa1+rxna1);
        b1 = Hs1\(diffRHSb1+rxnb1);

        a2 = Hs2\(diffRHSa2+rxna2);
        b2 = Hs2\(diffRHSb2+rxnb2);

        %% Plot the solution(s)
        % if mod(t,tplot) == 0
        if t==(Nt-1)
            scatplot=figure(ppp);
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
            % pause(1.0);

            % if vid==1
            %     currFrame = getframe(gcf);
            %     writeVideo(vidObj1,currFrame);
            % end

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
            % pause(1.0);

            % if vid==1
            %     currFrame = getframe(gcf);
            %     writeVideo(vidObj2,currFrame);
            % end

            if vid==1
                currframe = getframe(scatplot);
                writeVideo(vidObj1,currFrame);
            end

            %Plot on circle
            %Define colors
            colorLength = 50;
            white = [1,1,1];
            red = [1,0,0];
            blue = [0,0.5,1];
            maroon = [0.4,0,0];
            navy = [0,0,0.3];
            whitered = [linspace(white(1),red(1),colorLength)',linspace(white(2),red(2),colorLength)',linspace(white(3),red(3),colorLength)'];
            redmaroon = [linspace(red(1),maroon(1),colorLength)',linspace(red(2),maroon(2),colorLength)',linspace(red(3),maroon(3),colorLength)'];
            whiteredmaroon = [whitered;redmaroon];
            whiteblue = [linspace(white(1),blue(1),colorLength)',linspace(white(2),blue(2),colorLength)',linspace(white(3),blue(3),colorLength)'];
            bluenavy = [linspace(blue(1),navy(1),colorLength)',linspace(blue(2),navy(2),colorLength)',linspace(blue(3),navy(3),colorLength)'];
            whitebluenavy = [whiteblue; bluenavy];
            myColors = [linspace(red(1),blue(1),colorLength)',linspace(red(2),blue(2),colorLength)',linspace(red(3),blue(3),colorLength)'];
            redblue = abs(whiteblue+whitered)./2;
            redwhiteblue = [flip(whitered); whiteblue];

            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.93:0.01:1);
            [Xcol,Ycol] = pol2cart(th,rad);
            ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1]';
            ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1]';

            ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2]';
            ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2]';

            % figure(11); %Branched
            % % subplot(1,2,1);
            % contourf(Xcol,Ycol,Z1,100,'LineStyle','none')
            % axis square
            % col=colorbar;
            % colormap(whitered);
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
            % colormap(whiteblue);
            % clim([0,max(max(a),max(b))])
            % title('Bundled Actin')
            % ylabel(col,'Concentration')
            % set(gca,'XTick',[], 'YTick', [])

            %Plot both circles together
            % figure(13);
            % colFrame = subplot(1,2,1);
            % contourf(Xcol,Ycol,ZBranch1-ZBund1,100,'LineStyle','none')
            % axis square
            % col=colorbar;
            % colormap(redwhiteblue);
            % clim([-max(max(a1),max(b1)),max(max(a1),max(b1))])
            % title('Cell 1 Combined: Blue=Branched, Red=Bundled')
            % ylabel(col,'Concentration Branched - Concentration Bundled')
            % set(gca,'XTick',[], 'YTick', [])
            % 
            % [Ma,Ia] = max(a1);
            % [Mb,Ib] = max(b1);

            % hold on;
            % plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
            % plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
            % plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
            % hold off;

            % if vid==1
            %     currFrame = getframe(colFrame);
            %     writeVideo(vidObjCol1,currFrame);
            % end

            % colFrame = subplot(1,2,2);
            % contourf(Xcol,Ycol,ZBranch2-ZBund2,100,'LineStyle','none')
            % axis square
            % col=colorbar;
            % colormap(redwhiteblue);
            % clim([-max(max(a2),max(b2)),max(max(a2),max(b2))])
            % title('Cell 2 Combined: Blue=Branched, Red=Bundled')
            % ylabel(col,'Concentration Branched - Concentration Bundled')
            % set(gca,'XTick',[], 'YTick', [])
            % 
            % [Ma,Ia] = max(a2);
            % [Mb,Ib] = max(b2);

            % hold on;
            % plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
            % plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
            % plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
            % hold off;

            % if vid==1
            %     currFrame = getframe(colFrame);
            %     writeVideo(vidObjCol2,currFrame);
            % end

            % Concentric circles
            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
            [Xsm,Ysm] = pol2cart(th,rad);
            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.83:0.01:0.9);
            [Xmid,Ymid] = pol2cart(th,rad);

            % Cell 1
            % figc1=figure(16);
            figcells=figure(16);
            clf
            subplot(2,1,1)
            surf(Xcol,Ycol,ZBranch1);
            view(2)
            colormap(whitebluenavy)
            freezeColors;
            freezeColors(colorbar);
            clim([0,max(b1)])
            shading interp
            hold on;
            surf(Xmid,Ymid,ZBund1);
            colormap(whiteredmaroon)
            freezeColors;
            freezeColors(jicolorbar);
            clim([0,max(a1)])
            shading interp
            grid off
            axis square
            set(gca,'XTick',[], 'YTick', [])
            title('Cell 1 Combined: Blue=Branched, Red=Bundled')
            scatter(Xsm(boundC1),Ysm(boundC1),'black');
            hold off;
            box on;
            box off;
            set(gca,'XColor','w')
            set(gca,'YColor','w')
            set(gcf,'color','w');

            [Ma1,Ia1] = max(a1);
            [Mb1,Ib1] = max(b1);

            hold on;
            quiver(0,0,Xsm(1,Ia1),Ysm(1,Ia1),0,'color',[0 0 0])
            % plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
            % plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
            % plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
            hold off;
            
            % sinsum1=0;
            % cossum1=0;
            % for i=1:length(a1)
            %     sinsum1=sinsum1 + a1(i)*sin(th(i));
            %     cossum1=cossum1 + a1(i)*cos(th(i));
            % end
            % sinsum2=0;
            % cossum2=0;
            % for i=1:length(a2)
            %     sinsum2=sinsum2 + a2(i)*sin(th(i));
            %     cossum2=cossum2 + a2(i)*cos(th(i));
            % end
            % avg1=atan2(sinsum1,cossum1); %atan(sinsum1/cossum1);
            % avg2=atan2(sinsum2,cossum2); %atan(sinsum2/cossum2);

            % hold on;
            % quiver(0,0,0.8*cos(avg1),0.8*sin(avg1),0,'color',[1 0 0])
            % hold off;


            a1New = a1;
            a1New(a1New<1)=0;
            th2=th(1,:);
            if (a1New(1)~=0 && a1New(length(a1New))~=0)
                th2(ceil(length(th2)/2):length(th2))=-flip(th2(1:ceil(length(th2)/2)));
                % zeroIndex=find(a1New==0,1,'first');
                % a1New(zeroIndex:length(a1New))=-a1New(zeroIndex:length(a1New));
            end
            % avg1 = sum(a1New.*th(1,:)')/length(a1New);
            avg1=sum(a1New.*th2')/length(a1New);
            % hold on;
            % quiver(0,0,0.8*cos(avg1),0.8*sin(avg1),0,'color',[1 0 0])
            % hold off;

            if (a1New(1)~=0 && a1New(length(a1New))~=0)
                zeroInd1=find(a1New==0,1,'first');
                zeroInd2=find(a1New==0,1,'last');
                dirIndex=ceil((zeroInd1+zeroInd2)/2) - 50;
            else
                ind1=find(a1New~=0,1,'first');
                ind2=find(a1New~=0,1,'last');
                dirIndex=ceil((ind1+ind2)/2);
            end
            if dirIndex<1
                dirIndex=dirIndex+101;
            end
            hold on;
            quiver(0,0,Xsm(dirIndex),Ysm(dirIndex),0,'color',[1 0 0])
            hold off;

            % Cell 2
            % figc2=figure(17);
            % clf
            subplot(2,1,2)
            surf(Xcol,Ycol,ZBranch2);
            view(2)
            colormap(whitebluenavy)
            freezeColors;
            freezeColors(colorbar);
            clim([0,max(b2)])
            shading interp
            hold on;
            surf(Xmid,Ymid,ZBund2);
            colormap(whiteredmaroon)
            freezeColors;
            freezeColors(jicolorbar);
            clim([0,max(a2)])
            shading interp
            grid off
            axis square
            set(gca,'XTick',[], 'YTick', [])
            title('Cell 2 Combined: Blue=Branched, Red=Bundled')
            scatter(Xsm(boundC2),Ysm(boundC2),'black');
            hold off;
            box off;
            set(gca,'XColor','w')
            set(gca,'YColor','w')
            set(gcf,'color','w');

            [Ma2,Ia2] = max(a2);
            [Mb2,Ib2] = max(b2);
            p1=[0 0];
            p2=[Xsm(1,Ia2) Ysm(1,Ia2)];
            dp=p2-p1;

            hold on;
            quiver(p1(1),p1(2),dp(1),dp(2),0,'color',[0 0 0])
            % plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
            % plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
            % plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
            hold off;

            a2New = a2;
            a2New(a2New<1)=0;
            th2=th(1,:);
            if (a2New(1)~=0 && a2New(length(a1New))~=0)
                th2(ceil(length(th2)/2):length(th2))=-flip(th2(1:ceil(length(th2)/2)));
                % zeroIndex=find(a1New==0,1,'first');
                % a1New(zeroIndex:length(a1New))=-a1New(zeroIndex:length(a1New));
            end
            % avg2 = sum(a2New.*th(1,:)')/length(a2New);
            avg2=sum(a2New.*th2')/length(a2New);
            % hold on;
            % quiver(0,0,0.8*cos(avg2),0.8*sin(avg2),0,'color',[1 0 0])
            % hold off;

            % hold on;
            % quiver(0,0,0.8*cos(avg2),0.8*sin(avg2),0,'color',[1 0 0])
            % hold off;

            if (a2New(1)~=0 && a2New(length(a2New))~=0)
                zeroInd1=find(a2New==0,1,'first');
                zeroInd2=find(a2New==0,1,'last');
                dirIndex=ceil((zeroInd1+zeroInd2)/2) - 50;
            else
                ind1=find(a2New~=0,1,'first');
                ind2=find(a2New~=0,1,'last');
                dirIndex=ceil((ind1+ind2)/2);
            end
            if dirIndex<1
                dirIndex=dirIndex+101;
            end
            hold on;
            quiver(0,0,Xsm(dirIndex),Ysm(dirIndex),0,'color',[1 0 0])
            hold off;

            if vid==1
                currframe = getframe(figcells);
                writeVideo(vidObjCol1,currFrame);
            end

            % Calculate difference in direction angles
            ang1 = th(1,Ia1);
            ang2 = th(1,Ia2);
            angdiff = abs(ang1-ang2);
            if angdiff < pi/4
                samedirection='yes';
            elseif (abs(ang1-3*pi/2)<pi/4 && abs(ang2-pi/2)<pi/4)
                samedirection='strong no; collision';
            else
                samedirection='no';
            end
            sprintf('Angle difference: %d\nSame direction? %s',angdiff,samedirection)

            figure(16)
            hold on;
            
            % Plot Rac and Rho
            % figure(14);
            % colRRFrame = subplot(1,2,1);
            % ZRac1 = [xC1 xC1 xC1 xC1 xC1]';
            % ZRho1 = [yC1 yC1 yC1 yC1 yC1]';
            % contourf(Xcol,Ycol,ZRac1-ZRho1,100,'LineStyle','none');
            % axis square
            % col=colorbar;
            % colormap(redwhiteblue);
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
            % colormap(redwhiteblue);
            % clim([-max(max(xC2),max(yC2)),max(max(xC2),max(yC2))])
            % title('Cell 2 Combined: Blue=Rac, Red=Rho')
            % ylabel(col,'Concentration Rac - Concentration Rho')
            % set(gca,'XTick',[], 'YTick', [])

            % if vid==1
            %     currFrame = getframe(colRRFrame);
            %     writeVideo(vidObjRR2,currFrame);
            % end
        end
    end

    % measure of polarized state (1 if polarized and 0 otherwise)
    %st = 1*( (abs(a(1)-b(end))<1e-3 || abs(a(end)-b(1))<1e-3 ) && abs(a(1)-a(end))>1e-3 && abs(b(1)-b(end))>1e-3 );

    if vid==1
        close(vidObj1);
        close(vidObjCol1);
        % close(vidObjRR1);
        % 
        % close(vidObj2);
        % close(vidObjCol2);
        % close(vidObjRR2);
    end
    sprintf('Simulation %d done',ppp)
    toc
    if(quit_cond==0)
        ppp = ppp + 1;
    end
end

if savefigs==1
    % savefig(figc1,filenameC1);
    % savefig(figc2,filenameC2);
    savefig(figcells,filenameCells);
    savefig(scatplot,filenameScatter);
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