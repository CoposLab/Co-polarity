% Simulate competition for resources for actin filaments
% and fully coupled to stochastic biochemistry for polarity proteins
% on a circular cell membrane (periodic BC, 1.5D) with two cells
%
% Mechanics -- two actin networks: branched (a) and contractile (b)
% Polarity proteins : Rac (X) and Rho (Y)
%
% Last updated: 7/10/2023
% Katie Levandosky
% Calina Copos
addpath('./freeze_colors')
addpath('../SingleCellCode_Published')

clear;
close all;
clc;

savefigs=1;
setnum='6';
savelocation='./results2/movies1/firstCadherinImpl';
if savefigs==1
    filenameCells=strcat(savelocation,'Cells_',setnum);
    filenameScatter=strcat(savelocation,'Scatter_',setnum);
    filenameCads=strcat(savelocation,'Cadherins_',setnum);
end

vid = 1;
vidObj1 = VideoWriter(strcat(savelocation,'ScatterVid_',setnum,'.mp4'),'MPEG-4');
vidObjCol1 = VideoWriter(strcat(savelocation,'ColorVid_',setnum,'.mp4'),'MPEG-4');
vidObjCads = VideoWriter(strcat(savelocation, 'CadherinsVid_',setnum,'.mp4'),'MPEG-4');

counter_ppp = 1;
ppp = 1;

while (ppp<=1)
    counter_ppp = counter_ppp+1;

    clearvars -except counter_ppp vid vidObj1 ppp vidObjCol1 savefigs filenameScatter filenameCells vidObjCads filenameCads
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
    alpha = [2,2];
    beta = [2,2]; %first argument is away from overlap, second is on overlap

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
    Tend   = 25.0;                  % total simulation time
    Nt     = Tend/dt;
    dx     = sqrt(2*D*dt);
    tplot  = 100;

    posx1 = zeros(N,Nt);              % array of positions of X(t) cell 1 (rac)
    posy1 = zeros(N,Nt);              % array of positions of Y(t) cell 1 (rho)
    posz1 = zeros(N,Nt);              % array of positions of Z(t) cell 1 (cadherins)

    posx2 = zeros(N,Nt);              % array of positions of X(t) cell 2 (rac)
    posy2 = zeros(N,Nt);              % array of positions of Y(t) cell 2 (rho)
    posz2 = zeros(N,Nt);              % array of positions of Z(t) cell 2 (cadherins)

    epsilon=0.5; % distance to detect other molecules (finding nearby rac/rho to remove)
    numToRemove=0;
    counter1=0;
    counter2=0;

    % Boundary between cells
    blen = floor(bper * L); % length of overlap between cell membranes

    boundC1 = (floor((Na-1)*3/4 - floor((Na-1)*bper/2)))+1:(floor((Na-1)*3/4 + floor((Na-1)*bper/2)))+1;
    boundC2 = (floor((Na-1)*1/4 - floor((Na-1)*bper/2)))+1:(floor((Na-1)*1/4 + floor((Na-1)*bper/2)))+1;

    boundzC1 = 51:101;
    boundzC2 = 1:51;


    % Competition for limited resource (actin monomers) term
    %
    %F = @(U,V) -U.*U - m0*U.*V;
    F = @(U,V) -m0*U.*V;

    branchedConst1 = 1.0;
    bundledConst1 = 1.0;
    branchedConst2 = 1.0;
    bundledConst2 = 1.0;

    Ka1=ones(Na,1);
    Kb1=ones(Na,1);
    Ka2=ones(Na,1);
    Kb2=ones(Na,1);

    Ka1(boundC1) = branchedConst1*Ka1(boundC1);
    Kb1(boundC1) = bundledConst1*Kb1(boundC1);
    Ka2(boundC2) = branchedConst2*Ka2(boundC2);
    Kb2(boundC2) = bundledConst2*Kb2(boundC2);

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
    nx1   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active) (rac)
    ny1   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active) (rho)
    nz1   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active) (cadherins)
    Tx1   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for X(t) (rac)
    Ty1   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for Y(t) (rho)
    Tz1   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for Z(t) (cadherins)
    X1    = zeros(MAX_OUTPUT_LENGTH,1);       % number of X(t) molecules on the membrane (rac)
    Y1    = zeros(MAX_OUTPUT_LENGTH,1);       % number of Y(t) molecules on the membrane (rho)
    Z1    = zeros(MAX_OUTPUT_LENGTH,1);       % number of Z(t) molecules on the membrane (cadherins)
    NNx1  = zeros(Nt,1);
    NNy1  = zeros(Nt,1);
    NNz1  = zeros(Nt,1);

    nx2   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active) (rac)
    ny2   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active) (rho)
    nz2   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active) (cadherins)
    Tx2   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for X(t) (rac)
    Ty2   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for Y(t) (rho)
    Tz2   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for Z(t) (cadherins)
    X2    = zeros(MAX_OUTPUT_LENGTH,1);       % number of X(t) molecules on the membrane (rac)
    Y2    = zeros(MAX_OUTPUT_LENGTH,1);       % number of Y(t) molecules on the membrane (rho)
    Z2    = zeros(MAX_OUTPUT_LENGTH,1);       % number of Z(t) molecules on the membrane (cadherins)
    NNx2  = zeros(Nt,1);
    NNy2  = zeros(Nt,1);
    NNz2  = zeros(Nt,1);

    % Set initial conditions for polarity molecules distribution
    %
    rxn_count_x1       = 1;
    rxn_count_y1       = 1;
    rxn_count_z1       = 1;

    rxn_count_x2       = 1;
    rxn_count_y2       = 1;
    rxn_count_z2       = 1;

    X1(1)              = 0.1*N;                 % # of particles on membrane (rac)
    Y1(1)              = 0.1*N;                 % # of particles on membrane (rho)
    Z1(1)              = 0.1*N;                 % # of particles on membrane (cadherins)
    NNx1(1)            = X1(1);
    NNy1(1)            = Y1(1);
    NNz1(1)            = Y1(1);
    Tx1(1)             = 0.0;
    Ty1(1)             = 0.0;
    Tz1(1)             = 0.0;
    nx1(1:X1(1),1)      = 1;                     % activate mem-bound particles
    ny1(1:Y1(1),1)      = 1;
    nz1(1:Z1(1),1)      = 1;
    r1 = randperm(ceil(L/(0.0102)),X1(1)+Y1(1)+Z1(1))*0.0102;
    posx1(1:X1(1),1)=r1(1:X1(1));
    posy1(1:Y1(1),1)=r1(X1(1)+1:X1(1)+Y1(1));
    posz1(1:Z1(1),1)=r1(X1(1)+Y1(1)+1:end);
    % posz1(posz1*Na/L<boundzC1(1) | posz1*Na/L>boundzC1(end)) = 0;
    z_stable_idx1 = zeros(length(posz1),1); % 1 means stable, 0 means unstable; corresponds to indices of posz1

    X2(1)              = 0.1*N;                 % # of particles on membrane (rac)
    Y2(1)              = 0.1*N;                 % # of particles on membrane (rho)
    Z2(1)              = 0.1*N;                 % # of particles on membrane (cadherins)
    NNx2(1)            = X2(1);
    NNy2(1)            = Y2(1);
    NNz2(1)            = Z2(1);
    Tx2(1)             = 0.0;
    Ty2(1)             = 0.0;
    Tz2(1)             = 0.0;
    nx2(1:X2(1),1)      = 1;                     % activate mem-bound particles
    ny2(1:Y2(1),1)      = 1;
    nz2(1:Z2(1),1)      = 1;
    r2 = randperm(ceil(L/(0.0102)),X2(1)+Y2(1)+Z2(1))*0.0102;
    posx2(1:X2(1),1)=r2(1:X2(1));
    posy2(1:Y2(1),1)=r2(X2(1)+1:X2(1)+Y2(1));
    posz2(1:Z2(1),1)=r2(X2(1)+Y2(1)+1:end);
    % posz2(posz2*Na/L<boundzC2(1) | posz2*Na/L>boundzC2(end)) = 0;
    z_stable_idx2 = zeros(length(posz2),1); % 1 means stable, 0 means unstable; corresponds to indices of posz1


    % Sample concentration at actin filament spatial scale
    %
    % [s1,xC1,yC1] = resamplePolarityMolecules(posx1(1:NNx1(1),1),posy1(1:NNy1(1),1),NNx1(1),NNy1(1),L,Na);
    % [s2,xC2,yC2] = resamplePolarityMolecules(posx2(1:NNx2(1),1),posy2(1:NNy2(1),1),NNx2(1),NNy2(1),L,Na);

    [s1,xC1,yC1,zC1] = resamplePolarityMoleculesCad(posx1(1:NNx1(1),1),posy1(1:NNy1(1),1),posz1(1:NNz1(1),1),NNx1(1),NNy1(1),NNz1(1),L,Na);
    [s2,xC2,yC2,zC2] = resamplePolarityMoleculesCad(posx2(1:NNx2(1),1),posy2(1:NNy2(1),1),posz2(1:NNz2(1),1),NNx2(1),NNy2(1),NNz2(1),L,Na);

    % Check for stable cadherins
    epsilon2 = 0.1;
    for i=1:NNz1
        find(abs(posz1(i,1)-posz2(1:NNz2,1))<=epsilon2);
    end


    aic1 = a1;
    bic1 = b1;

    aic2 = a2;
    bic2 = b2;


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

        vidObjCads.FrameRate = 5;
        vidObjCads.Quality = 75;
        open(vidObjCads);
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
    % pause(1.0);


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

    if vid==1
        currframe = getframe(scatterFrame);
        writeVideo(vidObj1,currframe);
    end

    %Define colors
    colorLength = 50;
    white = [1,1,1];
    red = [1,0,0];
    blue = [143/256,177/256,221/256];
    maroon = [0.4,0,0];
    navy = [33/256,81/256,127/256];
    yellow = [1,0.9,0];
    darkyellow = [227/256,180/256,76/256];
    whitered = [linspace(white(1),red(1),colorLength)',linspace(white(2),red(2),colorLength)',linspace(white(3),red(3),colorLength)'];
    redmaroon = [linspace(red(1),maroon(1),colorLength)',linspace(red(2),maroon(2),colorLength)',linspace(red(3),maroon(3),colorLength)'];
    whiteredmaroon = [whitered;redmaroon];
    whiteblue = [linspace(white(1),blue(1),colorLength)',linspace(white(2),blue(2),colorLength)',linspace(white(3),blue(3),colorLength)'];
    bluenavy = [linspace(blue(1),navy(1),colorLength)',linspace(blue(2),navy(2),colorLength)',linspace(blue(3),navy(3),colorLength)'];
    whitebluenavy = [whiteblue; bluenavy];
    myColors = [linspace(red(1),blue(1),colorLength)',linspace(red(2),blue(2),colorLength)',linspace(red(3),blue(3),colorLength)'];
    redblue = abs(whiteblue+whitered)./2;
    redwhiteblue = [flip(whitered); whiteblue];
    whiteyellow = [linspace(white(1),yellow(1),colorLength)',linspace(white(2),yellow(2),colorLength)',linspace(white(3),yellow(3),colorLength)'];
    yellowdarkyellow = [linspace(yellow(1),darkyellow(1),colorLength)',linspace(yellow(2),darkyellow(2),colorLength)',linspace(yellow(3),darkyellow(3),colorLength)'];
    whitedarkyellow = [whiteyellow;yellowdarkyellow];

    % Define circles
    [th,rad] = meshgrid((0:3.6:360)*pi/180,0.93:0.01:1);
    [Xcol,Ycol] = pol2cart(th,rad);
    ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1]';
    ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1]';
    ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2]';
    ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2]';
    [th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
    [Xsm,Ysm] = pol2cart(th,rad);
    [th,rad] = meshgrid((0:3.6:360)*pi/180,0.83:0.01:0.9);
    [Xmid,Ymid] = pol2cart(th,rad);


    if vid == 1
        % Concentric circles
        % Cell 1
        figcells=figure(16);
        surf(Xcol,Ycol,ZBranch1);
        view(2)
        colormap(whitebluenavy)
        freezeColors;
        freezeColors(colorbar('Location','westoutside'));
        clim([0,max(max(b1),max(b2))])
        shading interp
        hold on;
        surf(Xmid,Ymid,ZBund1);
        colormap(whitedarkyellow)
        clim([0,max(max(a1),max(a2))])
        freezeColors;
        freezeColors(jicolorbar);
        shading interp
        grid off
        set(gca,'XTick',[], 'YTick', [])
        % scatter(Xsm(boundC1),Ysm(boundC1),'black');
        set(gca,'XColor','w')
        set(gca,'YColor','w')
        set(gcf,'color','w');



        % Cell 2
        surf(Xcol,Ycol-2,ZBranch2);
        view(2)
        colormap(whitebluenavy)
        freezeColors;
        freezeColors(colorbar('Location','westoutside'));
        clim([0,max(max(b1),max(b2))])
        shading interp
        surf(Xmid,Ymid-2,ZBund2);
        colormap(whitedarkyellow)
        freezeColors;
        freezeColors(jicolorbar);
        clim([0,max(max(a1),max(a2))])
        shading interp
        grid off
        axis equal
        set(gca,'XTick',[], 'YTick', [])
        title('Blue=Branched, Yellow=Bundled')
        % scatter(Xsm(boundC2),Ysm(boundC2)-2,'black');

        allmax = max(max(max(a1),max(a2)),max(max(b1),max(b2)));

        flipc2 = flip(boundC2);
        for i=1:length(boundC1)
            plot3([Xcol(end,boundC1(i)) Xcol(end,flipc2(i))], [Ycol(end,boundC1(i)) Ycol(end,flipc2(i))-2],[allmax+1,allmax+1],'black')
        end

        hold off;
        box off;
        set(gca,'XColor','w')
        set(gca,'YColor','w')
        set(gcf,'color','w')



        % Find median for cell 1
        a1New = a1;
        a1New(a1New<1)=0;
        if (a1New(1)~=0 && a1New(length(a1New))~=0)
            zeroInd1=find(a1New==0,1,'first');
            zeroInd2=find(a1New==0,1,'last');
            dirIndex1=ceil((zeroInd1+zeroInd2)/2) - 50;
        else
            ind1=find(a1New~=0,1,'first');
            ind2=find(a1New~=0,1,'last');
            dirIndex1=ceil((ind1+ind2)/2);
        end
        if dirIndex1<1
            dirIndex1=dirIndex1+101;
        end
        if ~isempty(dirIndex1)
            hold on;
            quiver(0,0,Xsm(dirIndex1),Ysm(dirIndex1),0,'color',[0 0 0])
            hold off;
        end

        % Find median for cell 2
        a2New = a2;
        a2New(a2New<1)=0;
        if (a2New(1)~=0 && a2New(length(a2New))~=0)
            zeroInd1=find(a2New==0,1,'first');
            zeroInd2=find(a2New==0,1,'last');
            dirIndex2=ceil((zeroInd1+zeroInd2)/2) - 50;
        else
            ind1=find(a2New~=0,1,'first');
            ind2=find(a2New~=0,1,'last');
            dirIndex2=ceil((ind1+ind2)/2);
        end
        if dirIndex2<1
            dirIndex2=dirIndex2+101;
        end
        if ~isempty(dirIndex2)
            hold on;
            quiver(0,-2,Xsm(dirIndex2),Ysm(dirIndex2),0,'color',[0 0 0])
            hold off;
        end

        currframe = getframe(figcells);
        writeVideo(vidObjCol1,currframe);

    end

    % Plot cadherins
    zC1s = zeros(length(zC1),1);
    zC2s = zeros(length(zC2),1);
    boundzC2flip = flip(boundzC2);
    for i=1:length(boundzC1)
        if zC1(boundzC1(i))>0.0001 && zC2(boundzC2flip(i))>0.0001
            zC1s(boundzC1(i))=zC1(boundzC1(i));
            zC2s(boundzC2flip(i))=zC2(boundzC2flip(i));
            counter1=counter1+1;
        else
            counter2=counter2+1;
        end
    end

    Zcads1 = [zC1 zC1 zC1 zC1 zC1 zC1 zC1 zC1]';
    Zcads2 = [zC2 zC2 zC2 zC2 zC2 zC2 zC2 zC2]';
    Zcads1s = [zC1s zC1s zC1s zC1s zC1s zC1s zC1s zC1s]';
    Zcads2s = [zC2s zC2s zC2s zC2s zC2s zC2s zC2s zC2s]';
    figcads=figure(2);
    clf
    hold on;
    surf(Xcol,Ycol,Zcads1s)
    surf(Xcol,Ycol-2,Zcads2s)
    colormap(whitedarkyellow)
    freezeColors;
    freezeColors(colorbar('Location','westoutside'));
    clim([0,max(max(max(zC1),max(zC2)),0.1)]);
    shading interp
    surf(Xmid,Ymid,Zcads1)
    surf(Xmid,Ymid-2,Zcads2)
    colormap(whitebluenavy)
    freezeColors;
    freezeColors(jicolorbar);
    clim([0,max(max(max(zC1s),max(zC2s)),0.1)]);
    hold off;
    view(2)
    shading interp
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gcf,'color','w')
    axis equal

    if vid==1
        currframe = getframe(figcads);
        writeVideo(vidObjCads,currframe);
    end

    figure(3)
    subplot(1,2,1)
    plot(s1,zC1)
    hold on
    scatter(posz1(:,1),ones(length(posz1(:,1)),1))
    hold off
    subplot(1,2,2)
    plot(s2,zC2)
    hold on
    scatter(posz2(:,1),ones(length(posz2(:,1)),1))
    hold off

    %% Run simulation
    %
    tic
    quit_cond = 0;
    cond = 0;
    for t=1:(Nt-1)

        %% Run biochemistry
        % [Konx1,Kony1,Kfbx1,Kfby1,Koffx1,Koffy1] = spatialrates(ron,rfb,roff,a1,b1,s1,beta,cond,boundC1); % set rates
        % [Konx2,Kony2,Kfbx2,Kfby2,Koffx2,Koffy2] = spatialrates(ron,rfb,roff,a2,b2,s2,beta,cond,boundC2);
        [Konx1,Kony1,Konz1,Kfbx1,Kfby1,Kfbz1,Koffx1,Koffy1,Koffz1] = spatialratesCad(ron,rfb,roff,a1,b1,s1,beta,cond,boundC1,boundzC1); % set rates
        [Konx2,Kony2,Konz2,Kfbx2,Kfby2,Kfbz2,Koffx2,Koffy2,Koffz2] = spatialratesCad(ron,rfb,roff,a2,b2,s2,beta,cond,boundC2,boundzC2);

        % Set konx and kony
        % Konx1(boundC1)=Konx1(boundC1)*100;
        % Konx2(boundC2)=Konx2(boundC2)*100;

        % Kony1(boundC1)=Kony1(boundC1)*10;
        % Kony2(boundC2)=Kony2(boundC2)*10;

        % Koffx1(boundC1)=Koffx1(boundC1)*10;
        % Koffx2(boundC2)=Koffx2(boundC2)*10;

        % Koffy1(boundC1)=Koffy1(boundC1)/100;
        % Koffy2(boundC2)=Koffy2(boundC2)*100;

        % Kfbx1(boundC1)=Kfbx1(boundC1)*100;
        % Kfbx2(boundC2)=Kfbx2(boundC2)*100;

        % Kfby1(boundC1)=Kfby1(boundC1)*100;
        % Kfby2(boundC2)=Kfby2(boundC2)*100;

        % Konx1(setdiff(1:length(Konx1),boundC1)) = Konx1(setdiff(1:length(Konx1),boundC1))*1000;
        % Kony2(setdiff(1:length(Kony2),boundC2)) = Kony2(setdiff(1:length(Kony2),boundC2))*1000;


        epsilon1 = 0.1;
        flipc2=flip(boundC2);
        scaledC1 = (L*boundC1/Na);
        scaledC2 = L*flipc2/Na;
        % for i=1:length(boundC1)
        %     sumx1 = sum(abs(posx1(:,t)-scaledC1(i))<=epsilon1);
        %     sumx2 = sum(abs(posx2(:,t)-scaledC2(i))<=epsilon1);
        %     sumy1 = sum(abs(posy1(:,t)-scaledC1(i))<=epsilon1);
        %     sumy2 = sum(abs(posy2(:,t)-scaledC2(i))<=epsilon1);
        %     if sumx1>0
        %         % Konx2(flipc2(i)) = Konx2(flipc2(i))*(sumx1*100);
        %         % Koffx2(flipc2(i)) = Koffx2(flipc2(i))*(sumx1*10);
        %         % Konx1(boundC1(i)) = Konx1(boundC1(i))/(sumx1*100);
        %         % Koffx1(boundC1(i)) = Koffx1(boundC1(i))*(sumx1*100);
        %         % Kony2(flipc2(i)) = Kony2(flipc2(i))*(sumx1*100);
        %         Kony1(boundC1(i)) = Kony1(boundC1(i))/(sumx1*100);
        %     end
        %     if sumx2>0
        %         % Konx1(boundC1(i)) = Konx1(boundC1(i))*(sumx2*100);
        %         % Koffx1(boundC1(i)) = Koffx1(boundC1(i))*(sumx2*10);
        %         % Konx2(flipc2(i)) = Konx2(flipc2(i))/(sumx2*100);
        %         % Koffx2(flipc2(i)) = Koffx2(flipc2(i))*(sumx2*100);
        %         % Kony1(boundC1(i)) = Kony1(boundC1(i))*(sumx2*100);
        %         Kony2(flipc2(i)) = Kony2(flipc2(i))/(sumx2*100);
        %     end
        %     if sumy1>0
        %         % Kony2(flipc2(i)) = Kony2(flipc2(i))/(sumy1*100);
        %         % Koffy2(flipc2(i)) = Koffy2(flipc2(i))*(sumy1*10);
        %         % Kony1(boundC1(i)) = Kony1(boundC1(i))/(sumy1*100);
        %         % Koffy1(boundC1(i)) = Koffy1(boundC1(i))*(sumy1*100);
        %     end
        %     if sumy2>0
        %         % Kony1(boundC1(i)) = Kony1(boundC1(i))/(sumy2*100);
        %         % Koffy1(boundC1(i)) = Koffy1(boundC1(i))*(sumy2*10);
        %         % Kony2(flipc2(i)) = Kony2(flipc2(i))/(sumy2*100);
        %         % Koffy2(flipc2(i)) = Koffy2(flipc2(i))*(sumy2*100);
        %     end
        % end

        % if max(a2)>0
        %     Kb1(boundC1) = 2*a2(flipc2)/max(a2)+1; % change bundled coeff in cell 1 proportionally to branched in cell 2
        % end
        % if max(a1)>0
        %     Kb2(flipc2) = 2*a1(boundC1)/max(a1)+1; % change bundled coeff in cell 2 proportionally to branched in cell 1
        % end
        % if max(b2)>0
        %     Ka1(boundC1) = 3*b2(flipc2)/max(b2)+1; % change branched coeff in cell 1 proportionally to bundled in cell 2
        % end
        % if max(b1)>0
        %     Ka2(flipc2) = 3*b1(boundC1)/max(b1)+1; % change branched coeff in cell 2 proportionally to bundled in cell 1
        % end
        %
        % Ka1(Ka1==0)=1;
        % Ka2(Ka2==0)=1;
        % Kb1(Kb1==0)=1;
        % Kb2(Kb2==0)=1;


        %Cell 1 rac
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

        %Cell 2 rac
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

        %Cell 1 rho
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

        %Cell 2 rho
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

        %Cell 1 cadherins
        if((t-1)*dt<Tz1(rxn_count_z1))
            NNz1(t+1) = Z1(rxn_count_z1-1);
        else
            nnz = Z1(rxn_count_z1);
            tauz = zeros(nnz,1);
            dn = zeros(nnz,1);
            r1 = rand(nnz,1);

            if(nnz==0)
                counter_ppp = ppp;
                quit_cond = 1;
                sprintf('here 1')
                break
            end

            for j=1:nnz          % all agents
                konz1 = interp1(s1,Konz1,posz1(j,t));
                koffz1 = interp1(s1,Koffz1,posz1(j,t));
                kfbz1 = interp1(s1,Kfbz1,posz1(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffz1 + (konz1+kfbz1*nnz/N)*(N/nnz-1);
                tauz(j) = -log(r1(j))/a0;
                rr = rand(1,1);
                dn(j) = (rr<((konz1+kfbz1*nnz/N)*(N/nnz-1)/a0))*1.0 + (rr>=((konz1+kfbz1*nnz/N)*(N/nnz-1)/a0))*(-1.0);
            end

            [mintauz1,minidz1] = min(tauz(1:j));       % find first chemical rxn
            Tz1(rxn_count_z1+1) = Tz1(rxn_count_z1) + mintauz1;
            Z1(rxn_count_z1+1) = nnz + dn(minidz1);
            rxn_count_z1 = rxn_count_z1 + 1;
            NNz1(t+1) = Z1(rxn_count_z1-1);
        end

        %Cell 2 cadherins
        if((t-1)*dt<Tz2(rxn_count_z2))
            NNz2(t+1) = Z2(rxn_count_z2-1);
        else
            nnz = Z2(rxn_count_z2);
            tauz = zeros(nnz,1);
            dn = zeros(nnz,1);
            r2 = rand(nnz,1);

            if(nnz==0)
                counter_ppp = ppp;
                quit_cond = 1;
                sprintf('here 2')
                break
            end

            for j=1:nnz          % all agents
                konz2 = interp1(s2,Konz2,posz2(j,t));
                koffz2 = interp1(s2,Koffz2,posz2(j,t));
                kfbz2 = interp1(s2,Kfbz2,posz2(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffz2 + (konz2+kfbz2*nnz/N)*(N/nnz-1);
                tauz(j) = -log(r2(j))/a0;
                rr = rand(1,1);
                dn(j) = (rr<((konz2+kfbz2*nnz/N)*(N/nnz-1)/a0))*1.0 + (rr>=((konz2+kfbz2*nnz/N)*(N/nnz-1)/a0))*(-1.0);
            end

            [mintauz2,minidz2] = min(tauz(1:j));       % find first chemical rxn
            Tz2(rxn_count_z2+1) = Tz2(rxn_count_z2) + mintauz2;
            Z2(rxn_count_z2+1) = nnz + dn(minidz2);
            rxn_count_z2 = rxn_count_z2 + 1;
            NNz2(t+1) = Z2(rxn_count_z2-1);
        end

        if (quit_cond==1)
            sprintf("It's happening at kk = %d, ppp = %d\n",kk,ppp)
            ppp = counter_ppp-1;
            break
        end

        %% Run diffusion of membrane-bound polarity proteins
        p  = 0.5;                  % probability of hoping left or right

        % Fetch the number of particles at this time
        K1_1 = NNx1(t+1); % rac
        K2_1 = NNy1(t+1); % rho
        K3_1 = NNz1(t+1); % cadherins

        K1_2 = NNx2(t+1); % rac
        K2_2 = NNy2(t+1); % rho
        K3_2 = NNz2(t+1); % cadherins

        % Between reactions, perform Brownian motion with periodic BC
        r1 = rand(K1_1,1);    % coin flip rac c1
        nx1(1:K1_1,t+1) = 1;
        posx1(1:K1_1,t+1) = posx1(1:K1_1,t) + dx*((r1<p)*1.0 + (r1>(1-p))*(-1.0));

        r2 = rand(K1_2,1);    % coin flip rac c2
        nx2(1:K1_2,t+1) = 1;
        posx2(1:K1_2,t+1) = posx2(1:K1_2,t) + dx*((r2<p)*1.0 + (r2>(1-p))*(-1.0));

        r1 = rand(K2_1,1);    % coin flip rho c1
        ny1(1:K2_1,t+1) = 1;
        posy1(1:K2_1,t+1) = posy1(1:K2_1,t) + dx*((r1<p)*1.0 + (r1>(1-p))*(-1.0));

        r2 = rand(K2_2,1);    % coin flip rh c2
        ny2(1:K2_2,t+1) = 1;
        posy2(1:K2_2,t+1) = posy2(1:K2_2,t) + dx*((r2<p)*1.0 + (r2>(1-p))*(-1.0));

        r1 = rand(K3_1,1);    % coin flip cadherins c1
        nz1(1:K3_1,t+1) = 1;
        posz1(1:K3_1,t+1) = posz1(1:K3_1,t) + dx*((r1<p)*1.0 + (r1>(1-p))*(-1.0));

        r2 = rand(K3_2,1);    % coin flip cadherins c2
        nz2(1:K3_2,t+1) = 1;
        posz2(1:K3_2,t+1) = posz2(1:K3_2,t) + dx*((r2<p)*1.0 + (r2>(1-p))*(-1.0));

        % Check for collision(s) and resolve any collisions
        % Resolution strategy: No one advances
        %
        % Cell 1
        firstcoll = sum(ismembertol(posx1(1:K1_1,t+1),posy1(1:K2_1,t+1),0.005,'DataScale',1)); % collisions between rac and rho
        if firstcoll~=0
            % Get indices of collisions
            aa = ismembertol(posx1(1:K1_1,t+1),posy1(1:K2_1,t+1),0.005,'DataScale',1);
            list_idx = find(aa~=0);
            bb = ismembertol(posy1(1:K2_1,t+1),posx1(1:K1_1,t+1),0.005,'DataScale',1);
            list_idy = find(bb~=0);

            posx1(list_idx,t+1) = posx1(list_idx,t);
            posy1(list_idy,t+1) = posy1(list_idy,t);
        end
        firstcoll = sum(ismembertol(posx1(1:K1_1,t+1),posz1(1:K3_1,t+1),0.005,'DataScale',1)); % collisions between rac and cadherins
        if firstcoll~=0
            % Get indices of collisions
            aa = ismembertol(posx1(1:K1_1,t+1),posz1(1:K3_1,t+1),0.005,'DataScale',1);
            list_idx = find(aa~=0);
            bb = ismembertol(posz1(1:K3_1,t+1),posx1(1:K1_1,t+1),0.005,'DataScale',1);
            list_idz = find(bb~=0);

            posx1(list_idx,t+1) = posx1(list_idx,t);
            posz1(list_idz,t+1) = posz1(list_idz,t);
        end
        firstcoll = sum(ismembertol(posy1(1:K2_1,t+1),posz1(1:K3_1,t+1),0.005,'DataScale',1)); % collisions between rho and cadherins
        if firstcoll~=0
            % Get indices of collisions
            aa = ismembertol(posy1(1:K2_1,t+1),posz1(1:K3_1,t+1),0.005,'DataScale',1);
            list_idy = find(aa~=0);
            bb = ismembertol(posz1(1:K3_1,t+1),posy1(1:K2_1,t+1),0.005,'DataScale',1);
            list_idz = find(bb~=0);

            posy1(list_idy,t+1) = posy1(list_idy,t);
            posz1(list_idz,t+1) = posz1(list_idz,t);
        end

        % Cell 2
        firstcoll = sum(ismembertol(posx2(1:K1_2,t+1),posy2(1:K2_2,t+1),0.005,'DataScale',1)); % collisions between rac and rho
        if firstcoll~=0
            % Get indices of collisions
            aa = ismembertol(posx2(1:K1_2,t+1),posy2(1:K2_2,t+1),0.005,'DataScale',1);
            list_idx = find(aa~=0);
            bb = ismembertol(posy2(1:K2_2,t+1),posx2(1:K1_2,t+1),0.005,'DataScale',1);
            list_idy = find(bb~=0);

            posx2(list_idx,t+1) = posx2(list_idx,t);
            posy2(list_idy,t+1) = posy2(list_idy,t);
        end
        firstcoll = sum(ismembertol(posx2(1:K1_2,t+1),posz2(1:K3_2,t+1),0.005,'DataScale',1)); % collisions between rac and cadherins
        if firstcoll~=0
            % Get indices of collisions
            aa = ismembertol(posx2(1:K1_2,t+1),posz2(1:K3_2,t+1),0.005,'DataScale',1);
            list_idx = find(aa~=0);
            bb = ismembertol(posz2(1:K3_2,t+1),posx2(1:K1_2,t+1),0.005,'DataScale',1);
            list_idz = find(bb~=0);

            posx2(list_idx,t+1) = posx2(list_idx,t);
            posz2(list_idz,t+1) = posz2(list_idz,t);
        end
        firstcoll = sum(ismembertol(posy2(1:K2_2,t+1),posz2(1:K3_2,t+1),0.005,'DataScale',1)); % collisions between rho and cadherins
        if firstcoll~=0
            % Get indices of collisions
            aa = ismembertol(posy2(1:K2_2,t+1),posz2(1:K3_2,t+1),0.005,'DataScale',1);
            list_idy = find(aa~=0);
            bb = ismembertol(posz2(1:K3_2,t+1),posy2(1:K2_2,t+1),0.005,'DataScale',1);
            list_idz = find(bb~=0);

            posy2(list_idy,t+1) = posy2(list_idy,t);
            posz2(list_idz,t+1) = posz2(list_idz,t);
        end

        % Enforce periodic boundary conditions
        posx1(1:K1_1,t+1) = posx1(1:K1_1,t+1) + (-L).*(posx1(1:K1_1,t+1)>L) + (L).*(posx1(1:K1_1,t+1)<0.0);
        posy1(1:K2_1,t+1) = posy1(1:K2_1,t+1) + (-L).*(posy1(1:K2_1,t+1)>L) + (L).*(posy1(1:K2_1,t+1)<0.0);
        posz1(1:K3_1,t+1) = posz1(1:K3_1,t+1) + (-L).*(posz1(1:K3_1,t+1)>L) + (L).*(posz1(1:K3_1,t+1)<0.0);

        posx2(1:K1_2,t+1) = posx2(1:K1_2,t+1) + (-L).*(posx2(1:K1_2,t+1)>L) + (L).*(posx2(1:K1_2,t+1)<0.0);
        posy2(1:K2_2,t+1) = posy2(1:K2_2,t+1) + (-L).*(posy2(1:K2_2,t+1)>L) + (L).*(posy2(1:K2_2,t+1)<0.0);
        posz2(1:K3_2,t+1) = posz2(1:K3_2,t+1) + (-L).*(posz2(1:K3_2,t+1)>L) + (L).*(posz2(1:K3_2,t+1)<0.0);

        % Enforce no-flux boundary conditions
        %posx(1:K1,t+1) = posx(1:K1,t+1) + (posx(1:K1,t)-posx(1:K1,t+1)).*(posx(1:K1,t+1)>L) + (posx(1:K1,t)-posx(1:K1,t+1)).*(posx(1:K1,t+1)<0.0);
        %posy(1:K2,t+1) = posy(1:K2,t+1) + (posy(1:K2,t)-posy(1:K2,t+1)).*(posy(1:K2,t+1)>L) + (posy(1:K2,t)-posy(1:K2,t+1)).*(posy(1:K2,t+1)<0.0);

        %% Determine if a biochemical rxn has occured - update positions

        % Find spontaneous association location cell 1
        ss1 = sort(posx1(1:K1_1,t)); % rac
        [ijk1] = find(ss1==posx1(minidx1,t),1);
        prevind1 = (ijk1-1)*(ijk1>1) + (K1_1)*(ijk1==1);
        nextind1 = (ijk1+1)*(ijk1<K1_1) + 1*(ijk1==K1_1);
        x2 = posx1(minidx1,t)+(ss1(prevind1)-posx1(minidx1,t))/2;
        x1 = posx1(minidx1,t)+(ss1(nextind1)-posx1(minidx1,t))/2;
        locx1 = (x2-x1).*rand(1,1) + x1; % random location halfway between the closest left/right particles
        
        ss1 = sort(posy1(1:K2_1,t)); % rho
        [ijk1] = find(ss1==posy1(minidy1,t),1);
        prevind1 = (ijk1-1)*(ijk1>1) + (K2_1)*(ijk1==1);
        nextind1 = (ijk1+1)*(ijk1<K2_1) + 1*(ijk1==K2_1);
        y2 = posy1(minidy1,t)+(ss1(prevind1)-posy1(minidy1,t))/2;
        y1 = posy1(minidy1,t)+(ss1(nextind1)-posy1(minidy1,t))/2;
        locy1 = (y2-y1).*rand(1,1) + y1; % random location halfway between the closest left/right particles

        ss1 = sort(posz1(1:K3_1,t)); % cadherins
        [ijk1] = find(ss1==posz1(minidz1,t),1);
        prevind1 = (ijk1-1)*(ijk1>1) + (K3_1)*(ijk1==1);
        nextind1 = (ijk1+1)*(ijk1<K3_1) + 1*(ijk1==K3_1);
        z2 = posz1(minidz1,t)+(ss1(prevind1)-posz1(minidz1,t))/2;
        z1 = posz1(minidz1,t)+(ss1(nextind1)-posz1(minidz1,t))/2;
        locz1 = (z2-z1).*rand(1,1) + z1; % random location halfway between the closest left/right particles

        ponx1 = ron/(ron+rfb*(N-K1_1));
        pony1 = ron/(ron+rfb*(N-K2_1));
        ponz1 = ron/(ron+rfb*(N-K3_1));

        % Find spontaneous association location cell 2
        ss2 = sort(posx2(1:K1_2,t)); % rac
        [ijk2] = find(ss2==posx2(minidx2,t),1);
        prevind2 = (ijk2-1)*(ijk2>1) + (K1_2)*(ijk2==1);
        nextind2 = (ijk2+1)*(ijk2<K1_2) + 1*(ijk2==K1_2);
        x2 = posx2(minidx2,t)+(ss2(prevind2)-posx2(minidx2,t))/2;
        x1 = posx2(minidx2,t)+(ss2(nextind2)-posx2(minidx2,t))/2;
        locx2 = (x2-x1).*rand(1,1) + x1; % random location halfway between the closest left/right particles

        ss2 = sort(posy2(1:K2_2,t)); % rho
        [ijk2] = find(ss2==posy2(minidy2,t),1);
        prevind2 = (ijk2-1)*(ijk2>1) + (K2_2)*(ijk2==1);
        nextind2 = (ijk2+1)*(ijk2<K2_2) + 1*(ijk2==K2_2);
        y2 = posy2(minidy2,t)+(ss2(prevind2)-posy2(minidy2,t))/2;
        y1 = posy2(minidy2,t)+(ss2(nextind2)-posy2(minidy2,t))/2;
        locy2 = (y2-y1).*rand(1,1) + y1; % random location halfway between the closest left/right particles

        ss2 = sort(posz2(1:K3_2,t)); % cadherins
        [ijk2] = find(ss2==posz2(minidz2,t),1);
        prevind2 = (ijk2-1)*(ijk2>1) + (K3_2)*(ijk2==1);
        nextind2 = (ijk2+1)*(ijk2<K3_2) + 1*(ijk2==K3_2);
        z2 = posz2(minidz2,t)+(ss2(prevind2)-posz2(minidz2,t))/2;
        z1 = posz2(minidz2,t)+(ss2(nextind2)-posz2(minidz2,t))/2;
        locz2 = (z2-z1).*rand(1,1) + z1; % random location halfway between the closest left/right particles

        ponx2 = ron/(ron+rfb*(N-K1_2));
        pony2 = ron/(ron+rfb*(N-K2_2));
        ponz2 = ron/(ron+rfb*(N-K3_2));

        %Cell 1 rac
        if(NNx1(t+1) < NNx1(t))                % diassociation event (particle off)
            oldcol = posx1(minidx1,1:end); % Find the particle to be removed
            othercols = posx1([1:minidx1-1,minidx1+1:K1_1],1:end); % Gather other "on" particles
            otherothercols = posx1(K1_1+1:end,1:end); % Gather "off" particles
            newpos = [othercols;oldcol;otherothercols]; % Put removed particle at the end of "on" particles
            posx1 = newpos;
            nx1(K1_1,t+1) = 0; % Set the removed particle to inactive
        elseif(NNx1(t+1) > NNx1(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posx1(K1_1,t+1) = posx1(K1_1,t)+(rr<ponx1)*locx1; % on event
            posx1(K1_1,t+1) = posx1(K1_1,t)+(rr>=ponx1)*posx1(minidx1,t);   % recruitment event
            nx1(K1_1,t+1) = 1;
            % Look for nearby rho (posy1), take them off
            % posx1(K1_1,t+1)=location of rac binding
            if numToRemove>0
                boundC1Scaled=(L*boundC1/Na);
                locRemovey1 = find(abs(posy1(:,t+1)-posx1(K1_1,t+1))<epsilon,numToRemove);
                numFound = length(locRemovey1);
                if ~isempty(locRemovey1) && boundC1Scaled(1)<=posx1(K1_1,t+1) && boundC1Scaled(end)>=posx1(K1_1,t+1)
                    % posy1(locRemovey1,t+1)=0;
                    oldcol = posy1(locRemovey1,1:end); % Find the particle(s) to be removed
                    othercols = posy1(setdiff(1:K2_1,locRemovey1),1:end); % Gather other "on" particles
                    otherothercols = posy1(K2_1+1:end,1:end); % Gather "off" particles
                    newpos = [othercols;oldcol;otherothercols]; % Put removed particle at the end of "on" particles
                    posy1 = newpos;
                    ny1(K2_1-numFound+1:K2_1,t+1) = 0;
                    counter1=counter1+numFound;
                end
            end
        end

        %Cell 2 rac
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
            nx2(K1_2,t+1) = 1;
            % Look for nearby rho (posy2), take them off
            % locx2=location of rac binding
            if numToRemove>0
                boundC2Scaled=(L*boundC2/Na);
                locRemovey2 = find(abs(posy2(:,t+1)-posx2(K1_2,t+1))<epsilon,numToRemove);
                numFound = length(locRemovey2);
                if ~isempty(locRemovey2) && boundC2Scaled(1)<=posx2(K1_2,t+1) && boundC2Scaled(end)>=posx2(K1_2,t+1)
                    % posy2(locRemovey2,t+1)=0;
                    oldcol = posy2(locRemovey2,1:end); % Find the particle to be removed
                    othercols = posy2(setdiff(1:K2_2,locRemovey2),1:end); % Gather other "on" particles
                    otherothercols = posy2(K2_2+1:end,1:end); % Gather "off" particles
                    newpos = [othercols;oldcol;otherothercols]; % Put removed particle at the end of "on" particles
                    posy2 = newpos;
                    ny2(K2_2-numFound+1:K2_2,t+1) = 0;
                    counter2=counter2+numFound;
                end
            end
        end

        %Cell 1 rho
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

        %Cell 2 rho
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

        %Cell 1 cadherins
        if (NNz1(t+1) < NNz1(t))                % diassociation event (particle off)
            oldcol = posz1(minidz1,1:end);
            othercols = posz1([1:minidz1-1,minidz1+1:K3_1],1:end);
            otherothercols = posz1(K3_1+1:end,1:end);
            newpos = [othercols;oldcol;otherothercols];
            posz1 = newpos;
            nz1(K3_1,t+1) = 0;
        elseif(NNz1(t+1) > NNz1(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posz1(K3_1,t+1) = posz1(K3_1,t)+(rr<ponz1)*locz1;               % on event
            posz1(K3_1,t+1) = posz1(K3_1,t)+(rr>=ponz1)*posz1(minidz1,t);    % recruitment event
            nz1(K3_1,t+1) = 1;
        end

        %Cell 2 cadherins
        if (NNz2(t+1) < NNz2(t))                % diassociation event (particle off)
            oldcol = posz2(minidz2,1:end);
            othercols = posz2([1:minidz2-1,minidz2+1:K3_2],1:end);
            otherothercols = posz2(K3_2+1:end,1:end);
            newpos = [othercols;oldcol;otherothercols];
            posz2 = newpos;
            nz2(K3_2,t+1) = 0;
        elseif(NNz2(t+1) > NNz2(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posz2(K3_2,t+1) = posz2(K3_2,t)+(rr<ponz2)*locz2;               % on event
            posz2(K3_2,t+1) = posz2(K3_2,t)+(rr>=ponz2)*posz2(minidz2,t);    % recruitment event
            nz2(K3_2,t+1) = 1;
        end

        % [s1,xC1,yC1] = resamplePolarityMolecules(posx1(1:K1_1,t+1),posy1(1:K2_1,t+1),K1_1,K2_1,L,Na);
        % [s2,xC2,yC2] = resamplePolarityMolecules(posx2(1:K1_2,t+1),posy2(1:K2_2,t+1),K1_2,K2_2,L,Na);
        [s1,xC1,yC1,zC1] = resamplePolarityMoleculesCad(posx1(1:K1_1,t+1),posy1(1:K2_1,t+1),posz1(1:K3_1,t+1),K1_1,K2_1,K3_1,L,Na);
        [s2,xC2,yC2,zC2] = resamplePolarityMoleculesCad(posx2(1:K1_2,t+1),posy2(1:K2_2,t+1),posz2(1:K3_2,t+1),K1_2,K2_2,K3_2,L,Na);


        %% Update actin filaments
        diffRHSa1 = Hm1*a1;
        diffRHSb1 = Hm1*b1;

        diffRHSa2 = Hm2*a2;
        diffRHSb2 = Hm2*b2;

        rxna1 = dt*( F(a1,b1) + Ka1.*(a1.*(1+alpha(1)*xC1) - a1.*a1)); %Cell 1 branched
        % rxna1(boundC1) = dt*( F(a1(boundC1),b1(boundC1)) + branchedConst1*(a1(boundC1).*(1+alpha(2)*xC1(boundC1)) - a1(boundC1).*a1(boundC1))); %Cell 1 branched on overlap

        rxnb1 = dt*( F(b1,a1) + Kb1.*(b1.*(1+alpha(1)*yC1) - b1.*b1)); %Cell 1 bundled
        % rxnb1(boundC1) = dt*( F(b1(boundC1),a1(boundC1)) + bundledConst1*(b1(boundC1).*(1+alpha(2)*yC1(boundC1)) - b1(boundC1).*b1(boundC1))); %Cell 1 bundled on overlap

        rxna2 = dt*( F(a2,b2) + Ka2.*(a2.*(1+alpha(1)*xC2) - a2.*a2)); %Cell 2 branched
        % rxna2(boundC2) = dt*( F(a2(boundC2),b2(boundC2)) + branchedConst2*(a2(boundC2).*(1+alpha(2)*xC2(boundC2)) - a2(boundC2).*a2(boundC2))); %Cell 2 branched on overlap

        rxnb2 = dt*( F(b2,a2) + Kb2.*(b2.*(1+alpha(1)*yC2) - b2.*b2)); %Cell 2 bundled
        % rxnb2(boundC2) = dt*( F(b2(boundC2),a2(boundC2)) + bundledConst2*(b2(boundC2).*(1+alpha(2)*yC2(boundC2)) - b2(boundC2).*b2(boundC2))); %Cell 2 bundled = dt*( F(b2,a2) + bundledConst2*(b2.*(1+alpha*yC2) - b2.*b2)); %Cell 2 bundled on overlap

        a1 = Hs1\(diffRHSa1+rxna1);
        b1 = Hs1\(diffRHSb1+rxnb1);

        a2 = Hs2\(diffRHSa2+rxna2);
        b2 = Hs2\(diffRHSb2+rxnb2);

        %% Plot the solution(s)
        if mod(t,tplot) == 0
        % if t==(Nt-1)
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

            if vid==1
                currframe = getframe(scatplot);
                writeVideo(vidObj1,currframe);
            end


            %Define colors
            colorLength = 50;
            white = [1,1,1];
            red = [1,0,0];
            blue = [143/256,177/256,221/256];
            maroon = [0.4,0,0];
            navy = [33/256,81/256,127/256];
            yellow = [1,0.9,0];
            darkyellow = [227/256,180/256,76/256];
            whitered = [linspace(white(1),red(1),colorLength)',linspace(white(2),red(2),colorLength)',linspace(white(3),red(3),colorLength)'];
            redmaroon = [linspace(red(1),maroon(1),colorLength)',linspace(red(2),maroon(2),colorLength)',linspace(red(3),maroon(3),colorLength)'];
            whiteredmaroon = [whitered;redmaroon];
            whiteblue = [linspace(white(1),blue(1),colorLength)',linspace(white(2),blue(2),colorLength)',linspace(white(3),blue(3),colorLength)'];
            bluenavy = [linspace(blue(1),navy(1),colorLength)',linspace(blue(2),navy(2),colorLength)',linspace(blue(3),navy(3),colorLength)'];
            whitebluenavy = [whiteblue; bluenavy];
            myColors = [linspace(red(1),blue(1),colorLength)',linspace(red(2),blue(2),colorLength)',linspace(red(3),blue(3),colorLength)'];
            redblue = abs(whiteblue+whitered)./2;
            redwhiteblue = [flip(whitered); whiteblue];
            whiteyellow = [linspace(white(1),yellow(1),colorLength)',linspace(white(2),yellow(2),colorLength)',linspace(white(3),yellow(3),colorLength)'];
            yellowdarkyellow = [linspace(yellow(1),darkyellow(1),colorLength)',linspace(yellow(2),darkyellow(2),colorLength)',linspace(yellow(3),darkyellow(3),colorLength)'];
            whitedarkyellow = [whiteyellow;yellowdarkyellow];

            % Define circles
            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.93:0.01:1);
            [Xcol,Ycol] = pol2cart(th,rad);
            ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1]';
            ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1]';
            ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2]';
            ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2]';
            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
            [Xsm,Ysm] = pol2cart(th,rad);
            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.83:0.01:0.9);
            [Xmid,Ymid] = pol2cart(th,rad);
            
            % Concentric circles
            % Cell 1
            figcells=figure(16);
            clf
            surf(Xcol,Ycol,ZBranch1);
            view(2)
            colormap(whitebluenavy)
            freezeColors;
            freezeColors(colorbar('Location','westoutside'));
            clim([0,max(max(b1),max(b2))])
            shading interp
            hold on;
            surf(Xmid,Ymid,ZBund1);
            colormap(whitedarkyellow)
            clim([0,max(max(a1),max(a2))])
            freezeColors;
            freezeColors(jicolorbar);
            shading interp
            grid off
            set(gca,'XTick',[], 'YTick', [])
            % scatter(Xsm(boundC1),Ysm(boundC1),'black');
            set(gca,'XColor','w')
            set(gca,'YColor','w')
            set(gcf,'color','w')


            % Cell 2
            surf(Xcol,Ycol-2,ZBranch2);
            view(2)
            colormap(whitebluenavy)
            freezeColors;
            freezeColors(colorbar('Location','westoutside'));
            clim([0,max(max(b1),max(b2))])
            shading interp
            surf(Xmid,Ymid-2,ZBund2);
            colormap(whitedarkyellow)
            freezeColors;
            freezeColors(jicolorbar);
            clim([0,max(max(a1),max(a2))])
            shading interp
            grid off
            axis equal
            set(gca,'XTick',[], 'YTick', [])
            title('Blue=Branched, Yellow=Bundled')
            % scatter(Xsm(boundC2),Ysm(boundC2)-2,'black');

            allmax = max(max(max(a1),max(a2)),max(max(b1),max(b2)));

            flipc2 = flip(boundC2);
            for i=1:length(boundC1)
                plot3([Xcol(end,boundC1(i)) Xcol(end,flipc2(i))], [Ycol(end,boundC1(i)) Ycol(end,flipc2(i))-2],[allmax+1,allmax+1],'black')
            end

            hold off;
            box off;
            set(gca,'XColor','w')
            set(gca,'YColor','w')
            set(gcf,'color','w')



            % Find median for cell 1
            a1New = a1;
            a1New(a1New<1)=0;
            if (a1New(1)~=0 && a1New(length(a1New))~=0)
                zeroInd1=find(a1New==0,1,'first');
                zeroInd2=find(a1New==0,1,'last');
                dirIndex1=ceil((zeroInd1+zeroInd2)/2) - 50;
            else
                ind1=find(a1New~=0,1,'first');
                ind2=find(a1New~=0,1,'last');
                dirIndex1=ceil((ind1+ind2)/2);
            end
            if dirIndex1<1
                dirIndex1=dirIndex1+101;
            end
            if ~isempty(dirIndex1)
                hold on;
                quiver(0,0,Xsm(dirIndex1),Ysm(dirIndex1),0,'color',[0 0 0])
                hold off;
            end



            % Find median for cell 2
            a2New = a2;
            a2New(a2New<1)=0;
            if (a2New(1)~=0 && a2New(length(a2New))~=0)
                zeroInd1=find(a2New==0,1,'first');
                zeroInd2=find(a2New==0,1,'last');
                dirIndex2=ceil((zeroInd1+zeroInd2)/2) - 50;
            else
                ind1=find(a2New~=0,1,'first');
                ind2=find(a2New~=0,1,'last');
                dirIndex2=ceil((ind1+ind2)/2);
            end
            if dirIndex2<1
                dirIndex2=dirIndex2+101;
            end
            if ~isempty(dirIndex2)
                hold on;
                quiver(0,-2,Xsm(dirIndex2),Ysm(dirIndex2),0,'color',[0 0 0])
                hold off;
            end

            % Add frame to video
            if vid==1
                currframe = getframe(figcells);
                writeVideo(vidObjCol1,currframe);
            end

            % Calculate difference in direction angles
            angTolerance=pi/4;
            strongAngTolerance=pi/5;
            if isempty(dirIndex1) && isempty(dirIndex2)
                samedirection='2NP';
                angdiff=NaN;
            elseif isempty(dirIndex1) || isempty(dirIndex2)
                samedirection='1NP';
                angdiff=NaN;
            else
                medang1 = th(1,dirIndex1);
                medang2 = th(1,dirIndex2);
                angdiff = min(abs(medang1-medang2),abs(2*pi-abs(medang1-medang2)));
                if angdiff < angTolerance
                    samedirection='yes';
                elseif (abs(medang1-3*pi/2)<strongAngTolerance && abs(medang2-pi/2)<strongAngTolerance)
                    samedirection='strong no; collision';
                else
                    samedirection='no';
                end
            end
            sprintf('Median angle difference: %d\nSame direction? %s',angdiff,samedirection)

            % Plot cadherins
            zC1s = zeros(length(zC1),1);
            zC2s = zeros(length(zC2),1);
            boundzC2flip = flip(boundzC2);
            for i=1:length(boundzC1)
                if zC1(boundzC1(i))>0.0001 && zC2(boundzC2flip(i))>0.0001
                    zC1s(boundzC1(i))=zC1(boundzC1(i));
                    zC2s(boundzC2flip(i))=zC2(boundzC2flip(i));
                    counter1=counter1+1;
                else
                    counter2=counter2+1;
                end
            end

            Zcads1 = [zC1 zC1 zC1 zC1 zC1 zC1 zC1 zC1]';
            Zcads2 = [zC2 zC2 zC2 zC2 zC2 zC2 zC2 zC2]';
            Zcads1s = [zC1s zC1s zC1s zC1s zC1s zC1s zC1s zC1s]';
            Zcads2s = [zC2s zC2s zC2s zC2s zC2s zC2s zC2s zC2s]';
            figcads=figure(2);
            clf
            hold on;
            surf(Xcol,Ycol,Zcads1s)
            surf(Xcol,Ycol-2,Zcads2s)
            hold off;
            colormap(whitedarkyellow)
            freezeColors;
            freezeColors(colorbar('Location','westoutside'));
            clim([0,max(max(max(zC1),max(zC2)),1)]);
            shading interp
            hold on;
            surf(Xmid,Ymid,Zcads1)
            surf(Xmid,Ymid-2,Zcads2)
            hold off;
            colormap(whitebluenavy)
            freezeColors;
            freezeColors(jicolorbar);
            clim([0,max(max(max(zC1s),max(zC2s)),1)]);
            view(2)
            shading interp
            set(gca,'XColor','w')
            set(gca,'YColor','w')
            set(gcf,'color','w')
            axis equal

            % flipc2 = flip(boundzC2);
            % for i=1:length(boundzC1)
            %     hold on
            %     plot3([Xcol(end,boundzC1(i)) Xcol(end,flipc2(i))], [Ycol(end,boundzC1(i)) Ycol(end,flipc2(i))-2],[1 1],'black')
            %     hold off
            % end

            if vid==1
                currframe = getframe(figcads);
                writeVideo(vidObjCads,currframe);
            end


            figure(3)
            subplot(1,2,1)
            plot(s1,zC1)
            hold on
            scatter(posz1(:,t),ones(length(posz1(:,t)),1))
            hold off
            subplot(1,2,2)
            plot(s2,zC2)
            hold on
            scatter(posz2(:,t),ones(length(posz2(:,t)),1))
            hold off
        end
    end

    % measure of polarized state (1 if polarized and 0 otherwise)
    %st = 1*( (abs(a(1)-b(end))<1e-3 || abs(a(end)-b(1))<1e-3 ) && abs(a(1)-a(end))>1e-3 && abs(b(1)-b(end))>1e-3 );

    if vid==1
        close(vidObj1);
        close(vidObjCol1);
        close(vidObjCads)
    end
    sprintf('Simulation %d done',ppp)
    sprintf('%d',quit_cond)
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
    savefig(figcads,filenameCads);
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