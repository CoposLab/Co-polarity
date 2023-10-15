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

%c1_vals = [1,10,1000];
%c2_vals = [1,10,1000];
% ka_vals=[-0.9,0,0.9];
% kb_vals=[-0.9,0,0.9];
% kc_vals=[-0.9,0,0.9];
% kd_vals=[-0.9,0,0.9];
%coeff_vals=[1,10,1000];

% all_results_matrix = zeros(length(c1_vals)*length(c2_vals),7);

% for ka_ind=2:2 %-0.9,0,0.9
%     for kb_ind=3:3 %-0.9,0,0.9
%         for kc_ind=1:1 %-0.9,0,0.9
%             for kd_ind=3:3 %-0.9,0,0.9
 %for c1_ind=4:4 %koffx,koffy,konx,kony
    %for c2_ind=1:1 %koffx,koffy,konx,kony
        %for c1coeff_ind=3:3 %1,10,1000
            %for c2coeff_ind=1:1 %1,10,1000
        
% for konx_ind=3:3
%     for koffx_ind=1:3
%         kony_ind=3;

polarize_time=0;
polarize_time_c1=0;
polarize_time_c2=0;
num_polarized=0;
num_pol_c1=0;
num_pol_c2=0;
countpol=0;
writem=0;
res_counters = [0,0,0,0,0,0,0]; %[yes, strong no, 1NP, 2NP, no, LF, dist. effort]

counter_ppp = 1;
ppp = 1;

while (ppp<=1)
    close all;
    savefigs=0;
    setnum=int2str(ppp);
    savelocation='';
    if savefigs==1
        % filenameC1=strcat('savedgraphs/doubleRhoOnCell1_',setnum);
        % filenameC2=strcat('savedgraphs/doubleRhoOnCell2_',setnum);
        filenameCells=strcat(savelocation,'Cells_',setnum);
        filenameScatter=strcat(savelocation,'Scatter_',setnum);
    end
    
    vid = 0;
    % vidObj1 = VideoWriter(strcat(savelocation,'ScatterVid',setnum,'.mp4'),'MPEG-4');
    % vidObjCol1 = VideoWriter(strcat(savelocation,'Colorvid',setnum,'.mp4'),'MPEG-4');
    % vidObjRR1 = VideoWriter('colorplotrr1.mp4','MPEG-4');
    % vidObj2 = VideoWriter('lineplot2.mp4','MPEG-4');
    % vidObjCol2 = VideoWriter('colorplot2.mp4','MPEG-4');
    % vidObjRR2 = VideoWriter('colorplotrr2.mp4','MPEG-4');

    counter_ppp = counter_ppp+1;

    clearvars -except counter_ppp vid vidObj1 ppp vidObjCol1 vidObjRR1 ... 
        vidObj2 vidObjCol2 vidObjRR2 savefigs filenameC1 filenameC2 ... 
        filenameScatter filenameCells res_counters c1_vals c2_vals c1_ind ... 
        c2_ind all_results_matrix polarize_time polarize_time_c1 ... 
        polarize_time_c2 num_polarized num_pol_c1 num_pol_c2 countpol writem ...
        ka_ind kb_ind kc_ind kd_ind ka_vals kb_vals kc_vals kd_vals coeff_vals ...
        konx_ind koffx_ind kony_ind c1coeff_ind c2coeff_ind
    
    rng('shuffle');
    set(0,'DefaultFigureVisible','on')

    polarizedc1=0; %has cell1 polarized yet
    polarizedc2=0;
    samedir_counter=0; %how many timesteps in a row have they polarized in the same direction


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
    tplot  = 50;

    posx1 = zeros(N,Nt);              % array of positions of X(t) cell 1
    posy1 = zeros(N,Nt);              % array of positions of Y(t) cell 1

    posx2 = zeros(N,Nt);              % array of positions of X(t) cell 2
    posy2 = zeros(N,Nt);              % array of positions of Y(t) cell 2

    epsilon=0.5; % distance to detect other molecules (finding nearby rac/rho to remove)
    numRhoToRemove=0;
    numRacToRemove=0;
    counter1=0;
    counter2=0;

    % Boundary between cells
    blen = floor(bper * L); % length of overlap between cell membranes
    % boundC1 = (ceil(3*Na/4 - ((Na-1)*bper)/2)):((ceil(3*Na/4 - ((Na-1)*bper)/2))+(Na-1)*bper); %boundary region in cell 1
    % boundC2 = (floor(Na/4 - ((Na-1)*bper)/2)):((floor(Na/4 - ((Na-1)*bper)/2))+(Na-1)*bper); %boundary region in cell 2

    boundC1 = (floor((Na-1)*3/4 - floor((Na-1)*bper/2)))+1:(floor((Na-1)*3/4 + floor((Na-1)*bper/2)))+1;
    boundC2 = (floor((Na-1)*1/4 - floor((Na-1)*bper/2)))+1:(floor((Na-1)*1/4 + floor((Na-1)*bper/2)))+1;

    % Signal
    signal=0;
    sigper=0.40;
    sigBound1 = (floor((Na-1)*1/8 - floor((Na-1)*sigper/2)))+1:(floor((Na-1)*1/8 + floor((Na-1)*sigper/2)))+1;
    sigBound1(sigBound1<=0)=sigBound1(sigBound1<=0)+Na;
    sigBound1(sigBound1>Na)=sigBound1(sigBound1>Na)-Na;
    sigBound2 = (floor((Na-1)*5/8 - floor((Na-1)*sigper/2)))+1:(floor((Na-1)*5/8 + floor((Na-1)*sigper/2)))+1;
    sigBound2(sigBound2<=0)=sigBound2(sigBound2<=0)+Na;
    sigBound2(sigBound2>Na)=sigBound2(sigBound2>Na)-Na;

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

    % Kb1(setdiff(1:length(Kb1),boundC1)) = 1.5*Kb1(setdiff(1:length(Kb1),boundC1));

    % Set initial conditions for actin distribution
    %
    ictype = 2;
    %    1 = step in the middle
    %    2 = random
    %    3 = sigmoidal curve
    %    4 = odd condition #1: branched peak in the middle, contractile peaks at
    %    the front and rear cell
    %    5 = odd condition #2: both peaks in the middle (w/ noise)

    a1       = zeros(length(Xa),1);
    anew1    = zeros(length(Xa),1);
    b1       = zeros(length(Xa),1);
    bnew1    = zeros(length(Xa),1);

    a2       = zeros(length(Xb),1);
    anew2    = zeros(length(Xb),1);
    b2       = zeros(length(Xb),1);
    bnew2    = zeros(length(Xb),1);


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

    a1all=zeros(length(Xa),Nt);
    a2all=zeros(length(Xb),Nt);
    b1all=zeros(length(Xa),Nt);
    b2all=zeros(length(Xb),Nt);
    a1all(:,1)=a1;
    a2all(:,1)=a2;
    b1all(:,1)=b1;
    b2all(:,1)=b2;

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
    conv1_1 = zeros(Nt,2);
    conv2_1 = zeros(Nt,2);
    convsup_1 = zeros(Nt,2);

    conv1_2 = zeros(Nt,2);
    conv2_2 = zeros(Nt,2);
    convsup_2 = zeros(Nt,2);

    %amount cells have moved
    xshift1=zeros(1,Nt);
    yshift1=zeros(1,Nt);
    xshift2=zeros(1,Nt);
    yshift2=zeros(1,Nt);


    % Set movie making
    %
    if vid==1
        vidObj1.FrameRate = 5;
        vidObj1.Quality = 75;
        open(vidObj1);

        vidObjCol1.FrameRate = 5;
        vidObjCol1.Quality = 75;
        open(vidObjCol1);
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
    yellow2 = [254/256,254/256,98/256];
    pink = [211/256,95/256,183/256];
    darkpink = [141/256,45/256,113/256];
    green = [26,255,26]/256;
    darkgreen = [16,150,16]/256;
    purple = [150,65,240]/256;
    darkpurple = [65,0,136]/256;
    orange = [230,97,0]/256;
    darkorange = [170,27,0]/256;


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
    whiteyellow2 = [linspace(white(1),yellow2(1),colorLength)',linspace(white(2),yellow2(2),colorLength)',linspace(white(3),yellow2(3),colorLength)'];
    yellow2darkyellow = [linspace(yellow2(1),darkyellow(1),colorLength)',linspace(yellow2(2),darkyellow(2),colorLength)',linspace(yellow2(3),darkyellow(3),colorLength)'];
    whitedarkyellow2 = [whiteyellow2;yellow2darkyellow];
    whitepink = [linspace(white(1),pink(1),colorLength)',linspace(white(2),pink(2),colorLength)',linspace(white(3),pink(3),colorLength)'];
    pinkdarkpink = [linspace(pink(1),darkpink(1),colorLength)',linspace(pink(2),darkpink(2),colorLength)',linspace(pink(3),darkpink(3),colorLength)'];
    whitedarkpink = [whitepink;pinkdarkpink];
    whitegreen = [linspace(white(1),green(1),colorLength)',linspace(white(2),green(2),colorLength)',linspace(white(3),green(3),colorLength)'];
    greendarkgreen = [linspace(green(1),darkgreen(1),colorLength)',linspace(green(2),darkgreen(2),colorLength)',linspace(green(3),darkgreen(3),colorLength)'];
    whitedarkgreen = [whitegreen;greendarkgreen];
    whitepurple = [linspace(white(1),purple(1),colorLength)',linspace(white(2),purple(2),colorLength)',linspace(white(3),purple(3),colorLength)'];
    purpledarkpurple = [linspace(purple(1),darkpurple(1),colorLength)',linspace(purple(2),darkpurple(2),colorLength)',linspace(purple(3),darkpurple(3),colorLength)'];
    whitedarkpurple = [whitepurple;purpledarkpurple];
    whiteorange = [linspace(white(1),orange(1),colorLength)',linspace(white(2),orange(2),colorLength)',linspace(white(3),orange(3),colorLength)'];
    orangedarkorange = [linspace(orange(1),darkorange(1),colorLength)',linspace(orange(2),darkorange(2),colorLength)',linspace(orange(3),darkorange(3),colorLength)'];
    whitedarkorange = [whiteorange;orangedarkorange];


    branchedColor = whitedarkpink;
    bundledColor = whitedarkyellow2;
    branchedColName = 'Pink';
    bundledColName = 'Yellow';

    % Define circles
    [th,rad] = meshgrid((0:3.6:360)*pi/180,0.93:0.01:1);
    [Xcol,Ycol] = pol2cart(th,rad);
    ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1]';
    ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1]';
    ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2]';
    ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2]';
    [th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
    [Xsm,Ysm] = pol2cart(th,rad);
    [th,rad] = meshgrid((0:3.6:360)*pi/180,0.86:0.01:0.93);
    [Xmid,Ymid] = pol2cart(th,rad);

    allmax = max(max(max(a1),max(a2)),max(max(b1),max(b2)));

    if vid == 1

        % Concentric circles
        % Cell 1
        figcells=figure(16);
        clf
        surf(Xcol,Ycol,ZBranch1,'AlphaData',ZBranch1+max(0,max(max(ZBranch2))-max(max(ZBranch1))),'FaceAlpha','interp','FaceColor','interp');
        view(2)
        colormap(branchedColor)
        freezeColors;
        freezeColors(colorbar('Location','westoutside'));
        clim([0,allmax])
        shading interp
        hold on;
        surf(Xcol,Ycol,ZBund1,'AlphaData',ZBund1+max(0,max(max(ZBund2))-max(max(ZBund1))),'FaceAlpha','interp','FaceColor','interp');
        colormap(bundledColor)
        clim([0,allmax])
        freezeColors;
        freezeColors(jicolorbar);
        shading interp
        grid off
        set(gca,'XTick',[], 'YTick', [])

        % Cell 2
        surf(Xcol,Ycol-2,ZBranch2,'AlphaData',ZBranch2+max(0,max(max(ZBranch1))-max(max(ZBranch2))),'FaceAlpha','interp','FaceColor','interp');
        view(2)
        colormap(branchedColor)
        freezeColors;
        freezeColors(colorbar('Location','westoutside'));
        clim([0,allmax])
        shading interp
        surf(Xcol,Ycol-2,ZBund2,'AlphaData',ZBund2+max(0,max(max(ZBund1))-max(max(ZBund2))),'FaceAlpha','interp','FaceColor','interp');
        colormap(bundledColor)
        freezeColors;
        freezeColors(jicolorbar);
        clim([0,allmax])
        shading interp
        grid off
        axis equal
        set(gca,'XTick',[], 'YTick', [])
        title(strcat(branchedColName, '=Branched, ', bundledColName, '=Bundled'))

        flipc2 = flip(boundC2);
        for i=1:length(boundC1)
            plot3([Xcol(end,boundC1(i)) Xcol(end,flipc2(i))], [Ycol(end,boundC1(i)) Ycol(end,flipc2(i))-2],[allmax+1,allmax+1],'black')
        end

        hold off;
        box off;
        set(gca,'XColor','w')
        set(gca,'YColor','w')
        set(gcf,'color','w');



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
            quiver(0,0,Xsm(dirIndex1),Ysm(dirIndex1),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5);
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
            quiver(0,-2,Xsm(dirIndex2),Ysm(dirIndex2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5)
            hold off;
        end

        % Plot signal
        if signal==1
            [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
            [Xsig,Ysig] = pol2cart(th,rad);
            if t<=1000
                hold on;
                scatter(Xsig(sigBound2),Ysig(sigBound2)-2,'black','.')
                hold off;
            else
                hold on;
                scatter(Xsig(sigBound1),Ysig(sigBound1),'black','.')
                hold off;
            end
        end

        ohf = findobj(gcf);
        figaxes = findobj(ohf(1), 'Type', 'axes');
        set(figaxes(1),'Fontsize',15)
        set(figaxes(2),'Fontsize',14)
        camroll(90)

        if vid==1
            currframe = getframe(figcells);
            writeVideo(vidObjCol1,currframe);
        end

    end

    % save(strcat('./vid_matfiles/sigswitch_vid_files/vars_t0'),'a1','b1','a2','b2','Xa','s1','s2','xC1','yC1','xC2','yC2','boundC1','boundC2','sigBound1','sigBound2');

    %% Run simulation
    %
    tic
    quit_cond = 0;
    cond = 0;
    for t=1:(Nt-1)

        %% Run biochemistry
        [Konx1,Kony1,Kfbx1,Kfby1,Koffx1,Koffy1] = spatialrates(ron,rfb,roff,a1,b1,s1,beta,cond,boundC1); % set rates
        [Konx2,Kony2,Kfbx2,Kfby2,Koffx2,Koffy2] = spatialrates(ron,rfb,roff,a2,b2,s2,beta,cond,boundC2);


        % Add external signal for cell 1
        % Konx2=ones(Na,1)*ron/8;
        % Kony2=ones(Na,1)*ron;
        % Koffx2=ones(Na,1)*roff;
        % Koffy2=ones(Na,1)*roff/8;
        % Kfbx2=ones(Na,1)*rfb/8;
        % Kfby2=ones(Na,1)*rfb;
        % 
        % Konx2(sigBound) = ron;
        % Kony2(sigBound) = ron/8;
        % Koffx2(sigBound) = roff/8;
        % Koffy2(sigBound) = roff;
        % Kfbx2(sigBound) = rfb;
        % Kfby2(sigBound) = rfb/8;

        % this works
        if signal==1
            steepness = 20;
            Konx2 = (ron*(tanh(steepness*(s2-s2(sigBound2(1)))) - tanh(steepness*(s2-s2(sigBound2(end)))) + 0.2)/2.2)';
            Kony2 = (ron*(2 - tanh(steepness*(s2-s2(sigBound2(1)))) + tanh(steepness*(s2-s2(sigBound2(end)))) + 0.2)/2.2)';
            Kfbx2 = (rfb*(tanh(steepness*(s2-s2(sigBound2(1)))) - tanh(steepness*(s2-s2(sigBound2(end)))) + 0.2)/2.2)';
            Kfby2 = (rfb*(2 - tanh(steepness*(s2-s2(sigBound2(1)))) + tanh(steepness*(s2-s2(sigBound2(end)))) + 0.2)/2.2)';
            Koffx2 = (roff*(2 - tanh(steepness*(s2-s2(sigBound2(1)))) + tanh(steepness*(s2-s2(sigBound2(end)))) + 0.2)/2.2)';
            Koffy2 = (roff*(tanh(steepness*(s2-s2(sigBound2(1)))) - tanh(steepness*(s2-s2(sigBound2(end)))) + 0.2)/2.2)';
        end

        % if signal==1
        %     steepness = 20;
        %     fx1=tanh(steepness*(s2-s2(sigBound(1))));
        %     fx2=tanh(steepness*(s2-s2(sigBound(end))));
        %     Konx2 = ron*ones(length(s2),1); %(ron*fx1+ron)/2 + (-fx2*ron+ron)/2;
        %     Kony2 = ((ron*(-fx1)+ron)/2 + (fx2*ron+ron)/2 + ron)/2;
        %     Kfbx2 = (rfb*fx1+rfb)/2 + (-fx2*rfb+rfb)/2;
        %     Kfby2 = ((rfb*(-fx1)+rfb)/2 + (fx2*rfb+rfb)/2 + rfb)/2;
        %     Koffx2 = ((roff*(-fx1)+roff)/2 + (fx2*roff+roff)/2 + roff)/2;
        %     Koffy2 = (roff*fx1+roff)/2 + (-fx2*roff+roff)/2;
        % end

         % steepness = 20;
         % Konx1 = ron*(tanh(steepness*(s1-1.875)) - tanh(steepness*(s1-5.625)) + 0.2)/2.2;
         % Kony1 = ron*(2 - tanh(steepness*(s1-1.875)) + tanh(steepness*(s1-5.625)) + 0.2)/2.2;
         % Kfbx1 = rfb*(tanh(steepness*(s1-1.875)) - tanh(steepness*(s1-5.625)) + 0.2)/2.2;
         % Kfby1 = rfb*(2 - tanh(steepness*(s1-1.875)) + tanh(steepness*(s1-5.625)) + 0.2)/2.2;
         % Koffx1 = roff*(2 - tanh(steepness*(s1-1.875)) + tanh(steepness*(s1-5.625)) + 0.2)/2.2;
         % Koffy1 = roff*(tanh(steepness*(s1-1.875)) - tanh(steepness*(s1-5.625)) + 0.2)/2.2;

         
         % Konx2=Konx2*10;
         % Kony2=Kony2*10;
         % Koffx2=Koffx2*10;
         % Koffy2=Koffy2*10;
         % Set konx and kony in contact region
         % Konx1(boundC1)=Konx1(boundC1)*1000;
         % Konx2(boundC2)=Konx2(boundC2)*100;
         
         % Kony1(boundC1)=Kony1(boundC1)*10;
         % Kony2(boundC2)=Kony2(boundC2)*1000;
         
         % Koffx1(boundC1)=Koffx1(boundC1)*10;
         % Koffx2(boundC2)=Koffx2(boundC2)*10;

         % Koffy1(boundC1)=Koffy1(boundC1)*1000;
         % Koffy2(boundC2)=Koffy2(boundC2)*1000;

         % Kfbx1(boundC1)=Kfbx1(boundC1)/10;
         % Kfbx2(boundC2)=Kfbx2(boundC2)/10;

         % Kfby1(boundC1)=Kfby1(boundC1)/10;
         % Kfby2(boundC2)=Kfby2(boundC2)/10;

         % Set konx and kony away from contact region
         % Konx1(setdiff(1:length(Konx1),boundC1)) = Konx1(setdiff(1:length(Konx1),boundC1))*1000;
         % Kony2(setdiff(1:length(Kony2),boundC2)) = Kony2(setdiff(1:length(Kony2),boundC2))*1000;
         % Kony1(setdiff(1:length(Kony1),boundC1)) = Kony1(setdiff(1:length(Kony1),boundC1))*100;
         % Konx2(setdiff(1:length(Konx2),boundC2)) = Konx2(setdiff(1:length(Konx2),boundC2))*100;

         % Set konx and kony in entire cell
         % Konx2=Konx2*10;
         % Kony2=Kony2*10;


        % Set konx and kony depending on rac/rho concentrations in contact
        % region
        % epsilon1 = 0.1;
        % flipc2=flip(boundC2);
        % scaledC1 = (L*boundC1/Na);
        % scaledC2 = L*flipc2/Na;
        % for i=1:length(boundC1)
        %     sumx1 = sum(abs(posx1(:,t)-scaledC1(i))<=epsilon1);
        %     sumx2 = sum(abs(posx2(:,t)-scaledC2(i))<=epsilon1);
        %     sumy1 = sum(abs(posy1(:,t)-scaledC1(i))<=epsilon1);
        %     sumy2 = sum(abs(posy2(:,t)-scaledC2(i))<=epsilon1);
        %     if sumx1>0
        %         % Konx2(flipc2(i)) = Konx2(flipc2(i))*(sumx1*1000);
        %         Koffx2(flipc2(i)) = Koffx2(flipc2(i))*(sumx1*1000);
        %         % Konx1(boundC1(i)) = Konx1(boundC1(i))/(sumx1*100);
        %         % Koffx1(boundC1(i)) = Koffx1(boundC1(i))*(sumx1*100);
        %         Kony2(flipc2(i)) = Kony2(flipc2(i))*(sumx1*1000);
        %         % Kony1(boundC1(i)) = Kony1(boundC1(i))/(sumx1*100);
        %     end
        %     if sumx2>0
        %         % Konx1(boundC1(i)) = Konx1(boundC1(i))*(sumx2*1000);
        %         Koffx1(boundC1(i)) = Koffx1(boundC1(i))*(sumx2*1000);
        %         % Konx2(flipc2(i)) = Konx2(flipc2(i))/(sumx2*100);
        %         % Koffx2(flipc2(i)) = Koffx2(flipc2(i))*(sumx2*100);
        %         Kony1(boundC1(i)) = Kony1(boundC1(i))*(sumx2*1000);
        %         % Kony2(flipc2(i)) = Kony2(flipc2(i))/(sumx2*100);
        %     end
        %     if sumy1>0
        %         % Kony2(flipc2(i)) = Kony2(flipc2(i))*(sumy1*1000);
        %         Koffy2(flipc2(i)) = Koffy2(flipc2(i))*(sumy1*1000);
        %         % Kony1(boundC1(i)) = Kony1(boundC1(i))/(sumy1*100);
        %         % Koffy1(boundC1(i)) = Koffy1(boundC1(i))*(sumy1*100);
        %         Konx2(flipc2(i)) = Konx2(flipc2(i))*(sumy1*1000);
        %     end
        %     if sumy2>0
        %         % Kony1(boundC1(i)) = Kony1(boundC1(i))*(sumy2*1000);
        %         Koffy1(boundC1(i)) = Koffy1(boundC1(i))*(sumy2*1000);
        %         % Kony2(flipc2(i)) = Kony2(flipc2(i))/(sumy2*100);
        %         % Koffy2(flipc2(i)) = Koffy2(flipc2(i))*(sumy2*100);
        %         Konx1(flipc2(i)) = Konx1(flipc2(i))*(sumy2*1000);
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

        % Set rac/rho rates depending on branched/bundled concentrations
        Konx1(boundC1) = Konx1(boundC1).*flip(b2(boundC2))*1000;
        % Konx2(boundC2) = Konx2(boundC2).*flip(b1(boundC1))*1000;
        % Kony1(boundC1) = Kony1(boundC1).*flip(a2(boundC2))*1000;
        Kony2(boundC2) = Kony2(boundC2).*flip(a1(boundC1))*1000;


        %Cell 1
        if((t-1)*dt<Tx1(rxn_count_x1))
            NNx1(t+1) = X1(rxn_count_x1-1);
        else
            nnx1 = X1(rxn_count_x1);
            taux1 = zeros(nnx1,1);
            dn1 = zeros(nnx1,1);
            r1 = rand(nnx1,1);

            if(nnx1==0)
                sprintf('here 1rac')
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nnx1          % all agents
                konx1 = interp1(s1,Konx1,posx1(j,t));
                koffx1 = interp1(s1,Koffx1,posx1(j,t));
                kfbx1 = interp1(s1,Kfbx1,posx1(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffx1 + (konx1+kfbx1*nnx1/N)*(N/nnx1-1);
                taux1(j) = -log(r1(j))/a0;
                rr = rand(1,1);
                dn1(j) = (rr<((konx1+kfbx1*nnx1/N)*(N/nnx1-1)/a0))*1.0 + (rr>=((konx1+kfbx1*nnx1/N)*(N/nnx1-1)/a0))*(-1.0);
            end

            [mintaux1,minidx1] = min(taux1(1:j));       % find first chemical rxn
            Tx1(rxn_count_x1+1) = Tx1(rxn_count_x1) + mintaux1;
            X1(rxn_count_x1+1) = nnx1 + dn1(minidx1);
            rxn_count_x1 = rxn_count_x1 + 1;
            NNx1(t+1) = X1(rxn_count_x1-1);
        end

        %Cell 2
        if((t-1)*dt<Tx2(rxn_count_x2))
            NNx2(t+1) = X2(rxn_count_x2-1);
        else
            nnx2 = X2(rxn_count_x2);
            taux2 = zeros(nnx2,1);
            dn2 = zeros(nnx2,1);
            r2 = rand(nnx2,1);

            if(nnx2==0)
                sprintf('here 2rac')
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nnx2          % all agents
                konx2 = interp1(s2,Konx2,posx2(j,t));
                koffx2 = interp1(s2,Koffx2,posx2(j,t));
                kfbx2 = interp1(s2,Kfbx2,posx2(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffx2 + (konx2+kfbx2*nnx2/N)*(N/nnx2-1);
                taux2(j) = -log(r2(j))/a0;
                rr = rand(1,1);
                dn2(j) = (rr<((konx2+kfbx1*nnx2/N)*(N/nnx2-1)/a0))*1.0 + (rr>=((konx2+kfbx2*nnx2/N)*(N/nnx2-1)/a0))*(-1.0);
            end

            [mintaux2,minidx2] = min(taux2(1:j));       % find first chemical rxn
            Tx2(rxn_count_x2+1) = Tx2(rxn_count_x2) + mintaux2;
            X2(rxn_count_x2+1) = nnx2 + dn2(minidx2);
            rxn_count_x2 = rxn_count_x2 + 1;
            NNx2(t+1) = X2(rxn_count_x2-1);
        end

        %Cell 1
        if((t-1)*dt<Ty1(rxn_count_y1))
            NNy1(t+1) = Y1(rxn_count_y1-1);
        else
            nny1 = Y1(rxn_count_y1);
            tauy1 = zeros(nny1,1);
            dn1 = zeros(nny1,1);
            r1 = rand(nny1,1);

            if(nny1==0)
                sprintf('here 1rho')
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nny1          % all agents
                kony1 = interp1(s1,Kony1,posy1(j,t));
                koffy1 = interp1(s1,Koffy1,posy1(j,t));
                kfby1 = interp1(s1,Kfby1,posy1(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffy1 + (kony1+kfby1*nny1/N)*(N/nny1-1);
                tauy1(j) = -log(r1(j))/a0;
                rr = rand(1,1);
                dn1(j) = (rr<((kony1+kfby1*nny1/N)*(N/nny1-1)/a0))*1.0 + (rr>=((kony1+kfby1*nny1/N)*(N/nny1-1)/a0))*(-1.0);
            end

            [mintauy1,minidy1] = min(tauy1(1:j));       % find first chemical rxn
            Ty1(rxn_count_y1+1) = Ty1(rxn_count_y1) + mintauy1;
            Y1(rxn_count_y1+1) = nny1 + dn1(minidy1);
            rxn_count_y1 = rxn_count_y1 + 1;
            NNy1(t+1) = Y1(rxn_count_y1-1);
        end

        %Cell 2
        if((t-1)*dt<Ty2(rxn_count_y2))
            NNy2(t+1) = Y2(rxn_count_y2-1);
        else
            nny2 = Y2(rxn_count_y2);
            tauy2 = zeros(nny2,1);
            dn2 = zeros(nny2,1);
            r2 = rand(nny2,1);

            if(nny2==0)
                sprintf('here 2rho')
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nny2          % all agents
                kony2 = interp1(s2,Kony2,posy2(j,t));
                koffy2 = interp1(s2,Koffy2,posy2(j,t));
                kfby2 = interp1(s2,Kfby2,posy2(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffy2 + (kony2+kfby2*nny2/N)*(N/nny2-1);
                tauy2(j) = -log(r2(j))/a0;
                rr = rand(1,1);
                dn2(j) = (rr<((kony2+kfby2*nny2/N)*(N/nny2-1)/a0))*1.0 + (rr>=((kony2+kfby2*nny2/N)*(N/nny2-1)/a0))*(-1.0);
            end

            [mintauy2,minidy2] = min(tauy2(1:j));       % find first chemical rxn
            Ty2(rxn_count_y2+1) = Ty2(rxn_count_y2) + mintauy2;
            Y2(rxn_count_y2+1) = nny2 + dn2(minidy2);
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
        r1_1 = rand(K1_1,1);    % coin flip
        nx1(1:K1_1,t+1) = 1;
        posx1(1:K1_1,t+1) = posx1(1:K1_1,t) + dx*((r1_1<p)*1.0 + (r1_1>(1-p))*(-1.0));

        r2_1 = rand(K1_2,1);    % coin flip
        nx2(1:K1_2,t+1) = 1;
        posx2(1:K1_2,t+1) = posx2(1:K1_2,t) + dx*((r2_1<p)*1.0 + (r2_1>(1-p))*(-1.0));

        r1_2 = rand(K2_1,1);    % coin flip
        ny1(1:K2_1,t+1) = 1;
        posy1(1:K2_1,t+1) = posy1(1:K2_1,t) + dx*((r1_2<p)*1.0 + (r1_2>(1-p))*(-1.0));

        r2_2 = rand(K2_2,1);    % coin flip
        ny2(1:K2_2,t+1) = 1;
        posy2(1:K2_2,t+1) = posy2(1:K2_2,t) + dx*((r2_2<p)*1.0 + (r2_2>(1-p))*(-1.0));

        % Check for collision(s) and resolve any collisions
        % Resolution strategy: No one advances
        %
        % Cell 1
        firstcoll = sum(ismembertol(posx1(1:K1_1,t+1),posy1(1:K2_1,t+1),0.005,'DataScale',1));
        if firstcoll~=0
            % Get indices of collisions
            aa1 = ismembertol(posx1(1:K1_1,t+1),posy1(1:K2_1,t+1),0.005,'DataScale',1);
            list_idx1 = find(aa1~=0);
            bb1 = ismembertol(posy1(1:K2_1,t+1),posx1(1:K1_1,t+1),0.005,'DataScale',1);
            list_idy1 = find(bb1~=0);

            posx1(list_idx1,t+1) = posx1(list_idx1,t);
            posy1(list_idy1,t+1) = posy1(list_idy1,t);
        end

        % Cell 2
        firstcoll = sum(ismembertol(posx2(1:K1_2,t+1),posy2(1:K2_2,t+1),0.005,'DataScale',1));
        if firstcoll~=0
            % Get indices of collisions
            aa2 = ismembertol(posx2(1:K1_2,t+1),posy2(1:K2_2,t+1),0.005,'DataScale',1);
            list_idx2 = find(aa2~=0);
            bb2 = ismembertol(posy2(1:K2_2,t+1),posx2(1:K1_2,t+1),0.005,'DataScale',1);
            list_idy2 = find(bb2~=0);

            posx2(list_idx2,t+1) = posx2(list_idx2,t);
            posy2(list_idy2,t+1) = posy2(list_idy2,t);
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
        x2_1 = posx1(minidx1,t)+(ss1(prevind1)-posx1(minidx1,t))/2;
        x1_1 = posx1(minidx1,t)+(ss1(nextind1)-posx1(minidx1,t))/2;
        locx1 = (x2_1-x1_1).*rand(1,1) + x1_1; % random location halfway between the closest left/right particles
        ss1 = sort(posy1(1:K2_1,t));
        [ijk1] = find(ss1==posy1(minidy1,t),1);
        prevind1 = (ijk1-1)*(ijk1>1) + (K2_1)*(ijk1==1);
        nextind1 = (ijk1+1)*(ijk1<K2_1) + 1*(ijk1==K2_1);
        y2_1 = posy1(minidy1,t)+(ss1(prevind1)-posy1(minidy1,t))/2;
        y1_1 = posy1(minidy1,t)+(ss1(nextind1)-posy1(minidy1,t))/2;
        locy1 = (y2_1-y1_1).*rand(1,1) + y1_1; % random location halfway between the closest left/right particles

        ponx1 = ron/(ron+rfb*(N-K1_1));
        pony1 = ron/(ron+rfb*(N-K2_1));

        % Find spontaneous association location cell 2
        ss2 = sort(posx2(1:K1_2,t));
        [ijk2] = find(ss2==posx2(minidx2,t),1);
        prevind2 = (ijk2-1)*(ijk2>1) + (K1_2)*(ijk2==1);
        nextind2 = (ijk2+1)*(ijk2<K1_2) + 1*(ijk2==K1_2);
        x2_2 = posx2(minidx2,t)+(ss2(prevind2)-posx2(minidx2,t))/2;
        x1_2 = posx2(minidx2,t)+(ss2(nextind2)-posx2(minidx2,t))/2;
        locx2 = (x2_2-x1_2).*rand(1,1) + x1_2; % random location halfway between the closest left/right particles
        ss2 = sort(posy2(1:K2_2,t));
        [ijk2] = find(ss2==posy2(minidy2,t),1);
        prevind2 = (ijk2-1)*(ijk2>1) + (K2_2)*(ijk2==1);
        nextind2 = (ijk2+1)*(ijk2<K2_2) + 1*(ijk2==K2_2);
        y2_2 = posy2(minidy2,t)+(ss2(prevind2)-posy2(minidy2,t))/2;
        y1_2 = posy2(minidy2,t)+(ss2(nextind2)-posy2(minidy2,t))/2;
        locy2 = (y2_2-y1_2).*rand(1,1) + y1_2; % random location halfway between the closest left/right particles

        ponx2 = ron/(ron+rfb*(N-K1_2));
        pony2 = ron/(ron+rfb*(N-K2_2));

        %Cell 1
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
            if numRhoToRemove>0
                boundC1Scaled=(L*boundC1/Na);
                locRemovey1 = find(abs(posy1(:,t+1)-posx1(K1_1,t+1))<epsilon,numRhoToRemove);
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
            nx2(K1_2,t+1) = 1;
            % Look for nearby rho (posy2), take them off
            % locx2=location of rac binding
            if numRhoToRemove>0
                boundC2Scaled=(L*boundC2/Na);
                locRemovey2 = find(abs(posy2(:,t+1)-posx2(K1_2,t+1))<epsilon,numRhoToRemove);
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

            if numRacToRemove>0
                boundC1Scaled=(L*boundC1/Na);
                locRemovex1 = find(abs(posx1(:,t+1)-posy1(K2_1,t+1))<epsilon,numRacToRemove);
                numFound = length(locRemovex1);
                if ~isempty(locRemovex1) && boundC1Scaled(1)<=posy1(K2_1,t+1) && boundC1Scaled(end)>=posy1(K2_1,t+1)
                    oldcol = posx1(locRemovex1,1:end); % Find the particle(s) to be removed
                    othercols = posx1(setdiff(1:K1_1,locRemovex1),1:end); % Gather other "on" particles
                    otherothercols = posx1(K1_1+1:end,1:end); % Gather "off" particles
                    newpos = [othercols;oldcol;otherothercols]; % Put removed particle at the end of "on" particles
                    posx1 = newpos;
                    nx1(K1_1-numFound+1:K1_1,t+1) = 0;
                    counter1=counter1+numFound;
                end
            end
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

            if numRacToRemove>0
                boundC2Scaled=(L*boundC2/Na);
                locRemovex2 = find(abs(posx2(:,t+1)-posy2(K2_2,t+1))<epsilon,numRacToRemove);
                numFound = length(locRemovex2);
                if ~isempty(locRemovex2) && boundC2Scaled(1)<=posy2(K2_2,t+1) && boundC2Scaled(end)>=posy2(K2_2,t+1)
                    oldcol = posx2(locRemovex2,1:end); % Find the particle(s) to be removed
                    othercols = posx2(setdiff(1:K1_2,locRemovex2),1:end); % Gather other "on" particles
                    otherothercols = posx2(K1_2+1:end,1:end); % Gather "off" particles
                    newpos = [othercols;oldcol;otherothercols]; % Put removed particle at the end of "on" particles
                    posx2 = newpos;
                    nx2(K1_2-numFound+1:K1_2,t+1) = 0;
                    counter2=counter2+numFound;
                end
            end
        end

        [s1,xC1,yC1] = resamplePolarityMolecules(posx1(1:K1_1,t+1),posy1(1:K2_1,t+1),K1_1,K2_1,L,Na);
        [s2,xC2,yC2] = resamplePolarityMolecules(posx2(1:K1_2,t+1),posy2(1:K2_2,t+1),K1_2,K2_2,L,Na);

        %% Update actin filaments
        diffRHSa1 = Hm1*a1;
        diffRHSb1 = Hm1*b1;

        diffRHSa2 = Hm2*a2;
        diffRHSb2 = Hm2*b2;

        kb1=zeros(length(b2),1);
        kb1(boundC1)=0.9*ones(length(boundC1),1);
        kc1=zeros(length(a2),1);
        kc1(boundC1)=0.9*ones(length(boundC1),1);
        kb2=zeros(length(b1),1);
        kb2(boundC2)=0.9*ones(length(boundC2),1);
        kc2=zeros(length(a1),1);
        kc2(boundC2)=0.9*ones(length(boundC2),1);
        abmax=50;

        gamma=1.5;

        rxna1 = dt*( F(a1,b1) + Ka1.*(a1.*(1+alpha(1)*xC1 + 0*kb1)) - a1.*a1); %Cell 1 branched
        rxnb1 = dt*( F(b1,a1) + Kb1.*(b1.*(1+alpha(1)*yC1 + 0*kc1)) - b1.*b1); %Cell 1 bundled
        rxna2 = dt*( F(a2,b2) + Ka2.*(a2.*(1+alpha(1)*xC2 + 0*kb2)) - a2.*a2); %Cell 2 branched
        rxnb2 = dt*( F(b2,a2) + Kb2.*(b2.*(1+alpha(1)*yC2 + 0*kc2)) - b2.*b2); %Cell 2 bundled

        % Growth term maxes out version
        % rxna1 = dt*( F(a1,b1) + Ka1.*(a1.*(1+alpha(1)*xC1 + kb1.* (flip(b2).*(flip(b2)<=abmax) + abmax*(flip(b2)>abmax)) ) - a1.*a1)); %Cell 1 branched
        % rxnb1 = dt*( F(b1,a1) + Kb1.*(b1.*(1+alpha(1)*yC1 + kc1.* (flip(a2).*(flip(a2)<=abmax) + abmax*(flip(a2)>abmax)) ) - b1.*b1)); %Cell 1 bundled
        % rxna2 = dt*( F(a2,b2) + Ka2.*(a2.*(1+alpha(1)*xC2 + kb2.* (flip(b1).*(flip(b1)<=abmax) + abmax*(flip(b1)>abmax)) ) - a2.*a2)); %Cell 2 branched
        % rxnb2 = dt*( F(b2,a2) + Kb2.*(b2.*(1+alpha(1)*yC2 + kc2.* (flip(a1).*(flip(a1)<=abmax) + abmax*(flip(a1)>abmax)) ) - b2.*b2)); %Cell 2 bundled

          % rxna1 = dt*( F(a1,b1) + (a1.*(1+alpha(1)*xC1 ... 
          %    + ka_vals(ka_ind) * kb1.* (flip(a2).*(flip(a2)<=abmax) + abmax*(flip(a2)>abmax)) ...
          %    + kb_vals(kb_ind) * kb1.* (flip(b2).*(flip(b2)<=abmax) + abmax*(flip(b2)>abmax)) ) - a1.*a1)); %Cell 1 branched
          % rxnb1 = dt*( F(b1,a1) + (b1.*(1+alpha(1)*yC1 ...
          %    + kc_vals(kc_ind) * kc1.* (flip(a2).*(flip(a2)<=abmax) + abmax*(flip(a2)>abmax)) ...
          %    + kd_vals(kd_ind) * kc1.* (flip(b2).*(flip(b2)<=abmax) + abmax*(flip(b2)>abmax)) ) - b1.*b1)); %Cell 1 bundled
          % rxna2 = dt*( F(a2,b2) + (a2.*(1+alpha(1)*xC2 ...
          %    + ka_vals(ka_ind) * kb2.* (flip(a1).*(flip(a1)<=abmax) + abmax*(flip(a1)>abmax)) ...
          %    + kb_vals(kb_ind) * kb2.* (flip(b1).*(flip(b1)<=abmax) + abmax*(flip(b1)>abmax)) ) - a2.*a2)); %Cell 2 branched
          % rxnb2 = dt*( F(b2,a2) + (b2.*(1+alpha(1)*yC2 ...
          %    + kc_vals(kc_ind) * kc2.* (flip(a1).*(flip(a1)<=abmax) + abmax*(flip(a1)>abmax)) ...
          %    + kd_vals(kd_ind) * kc2.* (flip(b1).*(flip(b1)<=abmax) + abmax*(flip(b1)>abmax)) ) - b2.*b2)); %Cell 2 bundled

        a1 = Hs1\(diffRHSa1+rxna1);
        b1 = Hs1\(diffRHSb1+rxnb1);

        a2 = Hs2\(diffRHSa2+rxna2);
        b2 = Hs2\(diffRHSb2+rxnb2);

        a1all(:,t)=a1;
        a2all(:,t)=a2;
        b1all(:,t)=b1;
        b2all(:,t)=b2;

        %Calculate direction angles
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
        [th,rad] = meshgrid((0:3.6:360)*pi/180,1);
        
        if ~isempty(dirIndex1)
            xshift1(t+1)=xshift1(t)+cos(th(dirIndex1))*0.0005;
            yshift1(t+1)=yshift1(t)+sin(th(dirIndex1))*0.0005;
        end
        if ~isempty(dirIndex2)
            xshift2(t+1)=xshift2(t)+cos(th(dirIndex2))*0.0005;
            yshift2(t+1)=yshift2(t)+sin(th(dirIndex2))*0.0005;
        end



        %% Plot the solution(s)
         % if mod(t,tplot) == 0
         if t==(Nt-1)

            %Define colors
            colorLength = 50;
            white = [1,1,1];
            red = [1,0,0];
            blue = [143/256,177/256,221/256];
            maroon = [0.4,0,0];
            navy = [33/256,81/256,127/256];
            yellow = [1,0.9,0];
            darkyellow = [227/256,180/256,76/256];
            yellow2 = [254/256,254/256,98/256];
            pink = [211/256,95/256,183/256];
            darkpink = [141/256,45/256,113/256];
            green = [26,255,26]/256;
            darkgreen = [16,150,16]/256;
            purple = [150,65,240]/256;
            darkpurple = [65,0,136]/256;
            orange = [230,97,0]/256;
            darkorange = [170,27,0]/256;


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
            whiteyellow2 = [linspace(white(1),yellow2(1),colorLength)',linspace(white(2),yellow2(2),colorLength)',linspace(white(3),yellow2(3),colorLength)'];
            yellow2darkyellow = [linspace(yellow2(1),darkyellow(1),colorLength)',linspace(yellow2(2),darkyellow(2),colorLength)',linspace(yellow2(3),darkyellow(3),colorLength)'];
            whitedarkyellow2 = [whiteyellow2;yellow2darkyellow];
            whitepink = [linspace(white(1),pink(1),colorLength)',linspace(white(2),pink(2),colorLength)',linspace(white(3),pink(3),colorLength)'];
            pinkdarkpink = [linspace(pink(1),darkpink(1),colorLength)',linspace(pink(2),darkpink(2),colorLength)',linspace(pink(3),darkpink(3),colorLength)'];
            whitedarkpink = [whitepink;pinkdarkpink];
            whitegreen = [linspace(white(1),green(1),colorLength)',linspace(white(2),green(2),colorLength)',linspace(white(3),green(3),colorLength)'];
            greendarkgreen = [linspace(green(1),darkgreen(1),colorLength)',linspace(green(2),darkgreen(2),colorLength)',linspace(green(3),darkgreen(3),colorLength)'];
            whitedarkgreen = [whitegreen;greendarkgreen];
            whitepurple = [linspace(white(1),purple(1),colorLength)',linspace(white(2),purple(2),colorLength)',linspace(white(3),purple(3),colorLength)'];
            purpledarkpurple = [linspace(purple(1),darkpurple(1),colorLength)',linspace(purple(2),darkpurple(2),colorLength)',linspace(purple(3),darkpurple(3),colorLength)'];
            whitedarkpurple = [whitepurple;purpledarkpurple];
            whiteorange = [linspace(white(1),orange(1),colorLength)',linspace(white(2),orange(2),colorLength)',linspace(white(3),orange(3),colorLength)'];
            orangedarkorange = [linspace(orange(1),darkorange(1),colorLength)',linspace(orange(2),darkorange(2),colorLength)',linspace(orange(3),darkorange(3),colorLength)'];
            whitedarkorange = [whiteorange;orangedarkorange];


            branchedColor = whitedarkpink;
            bundledColor = whitedarkyellow2;
            branchedColName = 'Pink';
            bundledColName = 'Yellow';


            % Make scatterplots
            scatplot=figure(ppp);
            subplot(1,2,1); %Cell 1
            plot(Xa,a1,'-o','color',branchedColor(end,:),'linewidth',3); hold on;
            plot(Xa,b1,'-ok','color',bundledColor(end,:),'linewidth',3);
            plot(s1,xC1,'-.','color',branchedColor(end,:),'linewidth',1);
            plot(s1,yC1,'-.k','color',bundledColor(end,:),'linewidth',1);
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
            plot(Xa,a2,'-o','color',branchedColor(end,:),'linewidth',3); hold on;
            plot(Xa,b2,'-ok','color',bundledColor(end,:),'linewidth',3);
            plot(s2,xC2,'-.','color',branchedColor(end,:),'linewidth',1);
            plot(s2,yC2,'-.k','color',bundledColor(end,:),'linewidth',1);
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


            % Define circles
            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.93:0.01:1);
            [Xcol,Ycol] = pol2cart(th,rad);
            ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1]';
            ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1]';
            ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2]';
            ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2]';
            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
            [Xsm,Ysm] = pol2cart(th,rad);
            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.86:0.01:0.93);
            [Xmid,Ymid] = pol2cart(th,rad);

            allmax = max(max(max(a1),max(a2)),max(max(b1),max(b2)));

            % Concentric circles
            % Cell 1
            figcells=figure(ppp+1);
            clf
            surf(Xcol+xshift1(t+1),Ycol+yshift1(t+1),ZBranch1,'AlphaData',ZBranch1+max(0,max(max(ZBranch2))-max(max(ZBranch1))),'FaceAlpha','interp','FaceColor','interp');
            view(2)
            colormap(branchedColor)
            freezeColors;
            freezeColors(colorbar('Location','westoutside'));
            clim([0,allmax])
            shading interp
            hold on;
            surf(Xcol+xshift1(t+1),Ycol+yshift1(t+1),ZBund1,'AlphaData',ZBund1+max(0,max(max(ZBund2))-max(max(ZBund1))),'FaceAlpha','interp','FaceColor','interp');
            colormap(bundledColor)
            clim([0,allmax])
            freezeColors;
            freezeColors(jicolorbar);
            shading interp
            grid off
            set(gca,'XTick',[], 'YTick', [])

            % Cell 2
            surf(Xcol+xshift2(t+1),Ycol-2+yshift2(t+1),ZBranch2,'AlphaData',ZBranch2+max(0,max(max(ZBranch1))-max(max(ZBranch2))),'FaceAlpha','interp','FaceColor','interp');
            view(2)
            colormap(branchedColor)
            freezeColors;
            freezeColors(colorbar('Location','westoutside'));
            clim([0,allmax])
            shading interp
            surf(Xcol+xshift2(t+1),Ycol-2+yshift2(t+1),ZBund2,'AlphaData',ZBund2+max(0,max(max(ZBund1))-max(max(ZBund2))),'FaceAlpha','interp','FaceColor','interp');
            colormap(bundledColor)
            freezeColors;
            freezeColors(jicolorbar);
            clim([0,allmax])
            shading interp
            grid off
            xlim([-3,3])
            ylim([-4,2])
            axis equal
            set(gca,'XTick',[], 'YTick', [])
            title(strcat(branchedColName, '=Branched, ', bundledColName, '=Bundled'))

            % flipc2 = flip(boundC2);
            % for i=1:length(boundC1)
            %     plot3([Xcol(end,boundC1(i)) Xcol(end,flipc2(i))], [Ycol(end,boundC1(i)) Ycol(end,flipc2(i))-2],[allmax+1,allmax+1],'black')
            % end

            hold off;
            box off;
            set(gca,'XColor','w')
            set(gca,'YColor','w')
            set(gcf,'color','w');



            % Plot arrows
            if ~isempty(dirIndex1)
                hold on;
                quiver(0+xshift1(t+1),0+yshift1(t+1),Xsm(dirIndex1),Ysm(dirIndex1),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5);
                hold off;
            end
            if ~isempty(dirIndex2)
                hold on;
                quiver(0+xshift2(t+1),-2+yshift2(t+1),Xsm(dirIndex2),Ysm(dirIndex2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5)
                hold off;
            end

            % Plot signal
            if signal==1
                [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
                [Xsig,Ysig] = pol2cart(th,rad);
                if t<=500
                    hold on;
                    scatter(Xsig(sigBound2),Ysig(sigBound2)-2,'black','.')
                    hold off;
                else
                    hold on;
                    scatter(Xsig(sigBound1),Ysig(sigBound1),'black','.')
                    hold off;
                end
            end

            ohf = findobj(gcf);
            figaxes = findobj(ohf(1), 'Type', 'axes');
            set(figaxes(1),'Fontsize',15)
            set(figaxes(2),'Fontsize',14)
            camroll(90)


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


            %save(strcat('./vid_matfiles/sigswitch_vid_files/vars_t',int2str(t)),'a1','b1','a2','b2','Xa','s1','s2','xC1','yC1','xC2','yC2','boundC1','boundC2','sigBound1','sigBound2');

        end

        % Check if polarized
        a1New = a1;
        a1New(a1New<0.0001)=0;
        if (a1New(1)~=0 && a1New(length(a1New))~=0)
            zeroInd1=find(a1New==0,1,'first');
            zeroInd2=find(a1New==0,1,'last');
            dirIndexa1=ceil((zeroInd1+zeroInd2)/2) - 50;
        else
            ind1=find(a1New~=0,1,'first');
            ind2=find(a1New~=0,1,'last');
            dirIndexa1=ceil((ind1+ind2)/2);
        end
        b1New = b1;
        b1New(b1New<0.0001)=0;
        if (b1New(1)~=0 && b1New(length(b1New))~=0)
            zeroInd1=find(b1New==0,1,'first');
            zeroInd2=find(b1New==0,1,'last');
            dirIndexb1=ceil((zeroInd1+zeroInd2)/2) - 50;
        else
            ind1=find(b1New~=0,1,'first');
            ind2=find(b1New~=0,1,'last');
            dirIndexb1=ceil((ind1+ind2)/2);
        end
        a2New = a2;
        a2New(a2New<0.0001)=0;
        if (a2New(1)~=0 && a2New(length(a2New))~=0)
            zeroInd1=find(a2New==0,1,'first');
            zeroInd2=find(a2New==0,1,'last');
            dirIndexa2=ceil((zeroInd1+zeroInd2)/2) - 50;
        else
            ind1=find(a2New~=0,1,'first');
            ind2=find(a2New~=0,1,'last');
            dirIndexa2=ceil((ind1+ind2)/2);
        end
        b2New = b2;
        b2New(b2New<0.0001)=0;
        if (b2New(1)~=0 && b2New(length(b2New))~=0)
            zeroInd1=find(b2New==0,1,'first');
            zeroInd2=find(b2New==0,1,'last');
            dirIndexb2=ceil((zeroInd1+zeroInd2)/2) - 50;
        else
            ind1=find(b2New~=0,1,'first');
            ind2=find(b2New~=0,1,'last');
            dirIndexb2=ceil((ind1+ind2)/2);
        end
        if dirIndexa1<1
                dirIndexa1=dirIndexa1+101;
        end
        if dirIndexb1<1
                dirIndexb1=dirIndexb1+101;
        end
        if dirIndexa2<1
                dirIndexa2=dirIndexa2+101;
        end
        if dirIndexb2<1
                dirIndexb2=dirIndexb2+101;
        end
        if countpol==1
            if polarizedc1==0 && ~isempty(dirIndexa1) && ~isempty(dirIndexb1)
                sprintf('Cell 1 polarized after %d steps', t)
                polarize_time_c1=polarize_time_c1+t;
                num_pol_c1=num_pol_c1+1;
                polarizedc1=1;
            end
            if polarizedc2==0 && ~isempty(dirIndexa2) && ~isempty(dirIndexb2)
                sprintf('Cell 2 polarized after %d steps', t)
                polarize_time_c2=polarize_time_c2+t;
                num_pol_c2=num_pol_c2+1;
                polarizedc2=1;
            end
            % if ~isempty(dirIndexa1) && ~isempty(dirIndexb1) && ~isempty(dirIndexa2) && ~isempty(dirIndexb2)
            %     sprintf('Both polarized after %d steps', t)
            %     polarize_time=polarize_time+t;
            %     num_polarized=num_polarized+1;
            %     break
            % end
            if ~isempty(dirIndexa1) && ~isempty(dirIndexb1) && ~isempty(dirIndexa2) && ~isempty(dirIndexb2)
                [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
                medang1 = th(1,dirIndexa1);
                medang2 = th(1,dirIndexa2);
                angdiff = min(abs(medang1-medang2),abs(2*pi-abs(medang1-medang2)));
                angTolerance=pi/4;
                if angdiff < angTolerance
                    samedir_counter=samedir_counter+1;
                else
                    samedir_counter=0;
                end
            else
                samedir_counter=0;
            end
            if samedir_counter>=10
                sprintf('Both polarized in same direction after %d steps', t)
                polarize_time=polarize_time+t;
                num_polarized=num_polarized+1;
                break
            end
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
        if savefigs==1
            % savefig(figc1,filenameC1);
            % savefig(figc2,filenameC2);
            savefig(figcells,filenameCells);
            savefig(scatplot,filenameScatter);
        end
        save(strcat('vid_matfiles/moving_cells/racupc1_rhoupc2_forcedependent/1000bRacOn_1000aRhoOn',int2str(ppp),'.mat'),...
            'boundC1','boundC2','posx1','posx2','posy1','posy2','NNx1','NNx2',...
            'NNy1','NNy2','a1all','a2all','b1all','b2all','Xa','Xb','s1','s2',...
            'xC1','xC2','yC1','yC2','xshift1','yshift1','xshift2','yshift2')
        ppp = ppp + 1;
        
        if writem==1
            if strcmp(samedirection, 'yes')
                res_counters(1)=res_counters(1)+1;
            elseif strcmp(samedirection, 'strong no; collision')
                res_counters(2)=res_counters(2)+1;
            elseif strcmp(samedirection, '1NP')
                res_counters(3)=res_counters(3)+1;
            elseif strcmp(samedirection, '2NP')
                res_counters(4)=res_counters(4)+1;
            else
                res_counters(5)=res_counters(5)+1;
            end

            anglelf=pi/4;
            if ~isempty(dirIndex1) && ~isempty(dirIndex2)
                if xor(abs(medang1-3*pi/2)<anglelf, abs(medang2-pi/2)<anglelf)
                    res_counters(6)=res_counters(6)+1;
                end
            end

            angledist=pi/4;
            if isempty(dirIndex1) && ~isempty(dirIndex2) && max(b1)>1
                medang2 = th(1,dirIndex2);
                if abs(medang2-pi/2)>angledist
                    res_counters(7)=res_counters(7)+1;
                end
            end
            if isempty(dirIndex2) && ~isempty(dirIndex1) && max(b2)>1
                medang1 = th(1,dirIndex1);
                if abs(medang1-3*pi/2)>angledist
                    res_counters(7)=res_counters(7)+1;
                end
            end
            if isempty(dirIndex1) && isempty(dirIndex2) && ((max(b1)>1 && max(a2)>1) || (max(b2)>1 && max(a1)>1))
                res_counters(7)=res_counters(7)+1;
            end
        end
    end

     if writem==1
         writematrix(res_counters,strcat('./allparamsresults/branchedbundled/',...
             string(ka_vals(ka_ind)),'ka_',string(kb_vals(kb_ind)),'kb_',...
             string(kc_vals(kc_ind)),'kc_',string(kd_vals(kd_ind)),'kd.xls'))
         % options=["Bkonx","Akony","Akoffx","Bkoffy"];
         % writematrix(res_counters,strcat('./allparamsresults/forcedependent/',...
         %     '1000',options(c1_ind),'C1_','1000',options(c2_ind),'C2.xls'))
         %options=["koffx","koffy","konx","kony"];
         %writematrix(res_counters,strcat('./allparamsresults/allracrho/',...
         %       string(coeff_vals(c1coeff_ind)), options(c1_ind), 'C1_',...
         %       string(coeff_vals(c2coeff_ind)), options(c2_ind), 'C2.xls'))
         sprintf(int2str(res_counters))
     end
end

if countpol==1
    % if num_polarized>0
    %     sprintf('Average polarization time for both cells with %d runs: %d timesteps',num_polarized,polarize_time/num_polarized)
    % end
    % if num_pol_c1>0
    %     sprintf('Average polarization time for cell 1 with %d runs: %d timesteps',num_pol_c1,polarize_time_c1/num_pol_c1)
    % end
    % if num_pol_c2>0
    %     sprintf('Average polarization time for cell 2 with %d runs: %d timesteps',num_pol_c2,polarize_time_c2/num_pol_c2)
    % end
    % if num_pol_c1+num_pol_c2>0
    %     sprintf('Average polarization time for one cell with %d runs: %d timesteps',num_pol_c1+num_pol_c2,(polarize_time_c1+polarize_time_c2)/(num_pol_c1+num_pol_c2))
    % end
    if num_polarized>0
        sprintf('Average polarization time for both cells in same direction with %d runs: %d timesteps',num_polarized,polarize_time/num_polarized)
    else
        sprintf('None polarized in the same direction')
    end
end

% all_results_matrix((c1_ind-1)*length(c1_vals)+c2_ind,:) = res_counters;


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
%end


 %            end
 %        end
 %    end
 % end

