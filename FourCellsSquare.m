% Simulate competition for resources for actin filaments
% and fully coupled to stochastic biochemistry for polarity proteins
% on a circular cell membrane (periodic BC, 1.5D) with two cells
%
% Mechanics -- two actin networks: branched (a) and contractile (b)
% Polarity proteins : Rac (X) and Rho (Y)
%
% Last updated: 4/30/2024
% Katie Levandosky
% Calina Copos
addpath('./freeze_colors')

clear;
close all;
clc;


save_matfile=0;
mat_location = '';

move_cells=0;

res_counters = [0,0,0,0,0,0,0]; %[yes, strong no, 1NP, 2NP, no, LF, dist. effort]

counter_ppp = 1;
ppp = 1;

while (ppp<=1)
    close all;
    savefigs=0;
    setnum=int2str(ppp);
    savelocation='';
    if savefigs==1
        filenameCells=strcat(savelocation,'Cells_',setnum);
        filenameScatter=strcat(savelocation,'Scatter_',setnum);
    end

    vid = 0;
    vidObjScatter = VideoWriter(strcat(savelocation,'ScatterVid_',setnum,'.mp4'),'MPEG-4');
    vidObjCells = VideoWriter(strcat(savelocation,'CellsVid_',setnum,'.mp4'),'MPEG-4');

    counter_ppp = counter_ppp+1;

    clearvars -except counter_ppp vid vidObjScatter ppp vidObjCells ...
        vidObj2 vidObjCol2 vidObjRR2 savefigs filenameC1 filenameC2 ...
        filenameScatter filenameCells res_counters c1_vals c2_vals c1_ind ...
        c2_ind all_results_matrix polarize_time polarize_time_c1 ...
        polarize_time_c2 num_polarized num_pol_c1 num_pol_c2 countpol writem ...
        polarize_time_c3 polarize_time_c4 num_pol_c3 num_pol_c4 save_matfile ...
        mat_location move_cells

    rng('shuffle');
    set(0,'DefaultFigureVisible','on')

    polarizedc1=0; %has cell1 polarized yet
    polarizedc2=0;
    polarizedc3=0;
    polarizedc4=0;


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

    posx3 = zeros(N,Nt);              % array of positions of X(t) cell 3
    posy3 = zeros(N,Nt);              % array of positions of Y(t) cell 3

    posx4 = zeros(N,Nt);              % array of positions of X(t) cell 4
    posy4 = zeros(N,Nt);              % array of positions of Y(t) cell 4

    

    % Boundary between cells
    blen = floor(bper * L); % length of overlap between cell membranes
    boundC1_4 = (floor((Na-1)*1/2 - floor((Na-1)*bper/2)))+1:(floor((Na-1)*1/2 + floor((Na-1)*bper/2)))+1;
    boundC1_2 = (floor((Na-1)*3/4 - floor((Na-1)*bper/2)))+1:(floor((Na-1)*3/4 + floor((Na-1)*bper/2)))+1;
    boundC2_1 = (floor((Na-1)*1/4 - floor((Na-1)*bper/2)))+1:(floor((Na-1)*1/4 + floor((Na-1)*bper/2)))+1;
    boundC2_3 = (floor((Na-1)*1/2 - floor((Na-1)*bper/2)))+1:(floor((Na-1)*1/2 + floor((Na-1)*bper/2)))+1;
    boundC3_2 = [(floor((Na-1)*1 - floor((Na-1)*bper/2)))+1:(floor((Na-1)*1)), 1:(floor((Na-1)*0 + floor((Na-1)*bper/2)))+1];
    boundC3_4 = (floor((Na-1)*1/4 - floor((Na-1)*bper/2)))+1:(floor((Na-1)*1/4 + floor((Na-1)*bper/2)))+1;
    boundC4_3 = (floor((Na-1)*3/4 - floor((Na-1)*bper/2)))+1:(floor((Na-1)*3/4 + floor((Na-1)*bper/2)))+1;
    boundC4_1 = [(floor((Na-1)*1 - floor((Na-1)*bper/2)))+1:(floor((Na-1)*1)), 1:(floor((Na-1)*0 + floor((Na-1)*bper/2)))+1];

    % Signal
    signal=0;
    sigper=0.40;
    sigBound1 = (floor((Na-1)*3/8 - floor((Na-1)*sigper/2)))+1:(floor((Na-1)*3/8 + floor((Na-1)*sigper/2)))+1;
    sigBound3 = (floor((Na-1)*5/8 - floor((Na-1)*sigper/2)))+1:(floor((Na-1)*5/8 + floor((Na-1)*sigper/2)))+1;

    % options: uncoupled, alt-racup-rhoup, alt-racup-rhoup-forcedependent,
    % branched-bundled-promotion
    pathway='branched-bundled-promotion';

    epsilon=0.5; % distance to detect other molecules (finding nearby rac/rho to remove)
    numToRemove=0;
    counter1=0;
    counter2=0;


    ka_vals=0.9*[0,1,-1];
    kb_vals=0.9*[0,1,-1];
    kc_vals=0.9*[0,1,-1];
    kd_vals=0.9*[0,1,-1];
    ka_ind=1;
    kb_ind=1;
    kc_ind=1;
    kd_ind=1;


    % Competition for limited resource (actin monomers) term
    %
    %F = @(U,V) -U.*U - m0*U.*V;
    F = @(U,V) -m0*U.*V;

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

    a3       = zeros(N,1);
    anew3    = zeros(N,1);
    b3       = zeros(N,1);
    bnew3    = zeros(N,1);

    a4       = zeros(N,1);
    anew4    = zeros(N,1);
    b4       = zeros(N,1);
    bnew4    = zeros(N,1);

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
        a3   = ones(N,1);
        anew3= ones(N,1);
        b3   = 0.5*ones(N,1);
        bnew3= 0.5*ones(N,1);
        a4   = ones(N,1);
        anew4= ones(N,1);
        b4   = 0.5*ones(N,1);
        bnew4= 0.5*ones(N,1);

        % (2) random
    elseif (ictype==2)
        a1 = 0.1 + 0.9.*rand(length(Xa),1);
        b1 = 0.1 + 0.9.*rand(length(Xa),1);
        a2 = 0.1 + 0.9.*rand(length(Xb),1);
        b2 = 0.1 + 0.9.*rand(length(Xb),1);
        a3 = 0.1 + 0.9.*rand(length(Xa),1);
        b3 = 0.1 + 0.9.*rand(length(Xa),1);
        a4 = 0.1 + 0.9.*rand(length(Xb),1);
        b4 = 0.1 + 0.9.*rand(length(Xb),1);

        % (3) arctangent
    elseif (ictype==3)
        steepness = 20;
        a1 = (tanh(steepness*(X1-0.375)) - tanh(steepness*(X1-1.125)) + 0.2)/2.2;
        b1 = (2 - tanh(steepness*(X1-0.375)) + tanh(steepness*(X1-1.125)) + 0.2)/2.2;
        a2 = (tanh(steepness*(X2-0.375)) - tanh(steepness*(X2-1.125)) + 0.2)/2.2;
        b2 = (2 - tanh(steepness*(X2-0.375)) + tanh(steepness*(X2-1.125)) + 0.2)/2.2;
        a3 = (tanh(steepness*(X3-0.375)) - tanh(steepness*(X3-1.125)) + 0.2)/2.2;
        b3 = (2 - tanh(steepness*(X3-0.375)) + tanh(steepness*(X3-1.125)) + 0.2)/2.2;
        a4 = (tanh(steepness*(X4-0.375)) - tanh(steepness*(X4-1.125)) + 0.2)/2.2;
        b4 = (2 - tanh(steepness*(X4-0.375)) + tanh(steepness*(X4-1.125)) + 0.2)/2.2;

        %a = (tanh(steepness*(X-0.5)) - tanh(steepness*(X-1.5)) + 0.2)/2.2;
        %b = (2 - tanh(steepness*(X-0.5)) + tanh(steepness*(X-1.5)) +0.2)/2.2;
    elseif (ictype==4)
        % (4) odd condition #1 (multiple peaks)
        steepness = 10;
        a1 = (1-cos(3*Xa*pi/5))/2; a1=a1';
        b1 = (tanh(steepness*(Xa-7.5))+1)/2 + (1-tanh(steepness*(Xa-2.5)))/2; b1=b1';
        a2 = (1-cos(3*Xb*pi/5))/2; a2=a2';
        b2 = (tanh(steepness*(Xb-7.5))+1)/2 + (1-tanh(steepness*(Xb-2.5)))/2; b2=b2';
        a3 = (1-cos(3*Xa*pi/5))/2; a3=a3';
        b3 = (tanh(steepness*(Xa-7.5))+1)/2 + (1-tanh(steepness*(Xa-2.5)))/2; b3=b3';
        a4 = (1-cos(3*Xb*pi/5))/2; a4=a4';
        b4 = (tanh(steepness*(Xb-7.5))+1)/2 + (1-tanh(steepness*(Xb-2.5)))/2; b4=b4';

    elseif (ictype==5)
        % (4) odd condition #2
        steepness = 10;
        b1 = (1-cos(Xa*pi/5))/2; b1=b1';
        a1 = (tanh(steepness*(Xa-7.5))+1)/2 + (1-tanh(steepness*(Xa-2.5)))/2; a1=a1';
        b2 = (1-cos(Xb*pi/5))/2; b2=b2';
        a2 = (tanh(steepness*(Xb-7.5))+1)/2 + (1-tanh(steepness*(Xb-2.5)))/2; a2=a2';
        b3 = (1-cos(Xa*pi/5))/2; b3=b3';
        a3 = (tanh(steepness*(Xa-7.5))+1)/2 + (1-tanh(steepness*(Xa-2.5)))/2; a3=a3';
        b4 = (1-cos(Xb*pi/5))/2; b4=b4';
        a4 = (tanh(steepness*(Xb-7.5))+1)/2 + (1-tanh(steepness*(Xb-2.5)))/2; a4=a4';

    elseif (ictype==6)
        % (5) odd condition #3
        mu = 1.8; sigma = 0.1;
        a1 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); a1=a1';
        a2 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); a2=a2';
        a3 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); a3=a3';
        a4 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); a4=a4';
        mu = 1.9; sigma = 0.1;
        b1 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); b1=b1';
        b2 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); b2=b2';
        b3 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); b3=b3';
        b4 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); b4=b4';
    end

    a1all=zeros(length(Xa),Nt);
    a2all=zeros(length(Xb),Nt);
    b1all=zeros(length(Xa),Nt);
    b2all=zeros(length(Xb),Nt);
    a3all=zeros(length(Xa),Nt);
    a4all=zeros(length(Xb),Nt);
    b3all=zeros(length(Xa),Nt);
    b4all=zeros(length(Xb),Nt);
    a1all(:,1)=a1;
    a2all(:,1)=a2;
    b1all(:,1)=b1;
    b2all(:,1)=b2;
    a3all(:,1)=a3;
    a4all(:,1)=a4;
    b3all(:,1)=b3;
    b4all(:,1)=b4;

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
    II3 = speye(Na,Na);
    Lapdiff3 = spdiags(ones(Na,1)*[0 1 -2 1 0], [-Na+1, -1, 0 , 1, Na-1], Na, Na);
    Lapdiff3(1,1) = -2; Lapdiff3(1,2) = 1; Lapdiff3(1,Na) = 1;
    Lapdiff3(Na,1) = 1; Lapdiff3(Na,Na-1) = 1; Lapdiff3(Na,Na) = -2;
    Hm3 = II3+(pa/2)*Lapdiff3;
    Hs3 = II3-(pa/2)*Lapdiff3;
    II4 = speye(Na,Na);
    Lapdiff4 = spdiags(ones(Na,1)*[0 1 -2 1 0], [-Na+1, -1, 0 , 1, Na-1], Na, Na);
    Lapdiff4(1,1) = -2; Lapdiff4(1,2) = 1; Lapdiff4(1,Na) = 1;
    Lapdiff4(Na,1) = 1; Lapdiff4(Na,Na-1) = 1; Lapdiff4(Na,Na) = -2;
    Hm4 = II4+(pa/2)*Lapdiff4;
    Hs4 = II4-(pa/2)*Lapdiff4;

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

    nx3   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
    ny3   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
    Tx3   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for X(t)
    Ty3   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for Y(t)
    X3    = zeros(MAX_OUTPUT_LENGTH,1);       % number of X(t) molecules on the membrane
    Y3    = zeros(MAX_OUTPUT_LENGTH,1);       % number of Y(t) molecules on the membrane
    NNx3  = zeros(Nt,1);
    NNy3  = zeros(Nt,1);

    nx4   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
    ny4   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
    Tx4   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for X(t)
    Ty4   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for Y(t)
    X4    = zeros(MAX_OUTPUT_LENGTH,1);       % number of X(t) molecules on the membrane
    Y4    = zeros(MAX_OUTPUT_LENGTH,1);       % number of Y(t) molecules on the membrane
    NNx4  = zeros(Nt,1);
    NNy4  = zeros(Nt,1);

    % Set initial conditions for polarity molecules distribution
    %
    rxn_count_x1       = 1;
    rxn_count_y1       = 1;
    rxn_count_x2       = 1;
    rxn_count_y2       = 1;
    rxn_count_x3       = 1;
    rxn_count_y3       = 1;
    rxn_count_x4       = 1;
    rxn_count_y4       = 1;

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

    X3(1)              = 0.1*N;                 % # of particles on membrane
    Y3(1)              = 0.1*N;                 % # of particles on membrane
    NNx3(1)            = X3(1);
    NNy3(1)            = Y3(1);
    Tx3(1)             = 0.0;
    Ty3(1)             = 0.0;
    nx3(1:X3(1),1)      = 1;                     % activate mem-bound particles
    ny3(1:X3(1),1)      = 1;
    r3 = randperm(ceil(L/(0.0102)),X3(1)+Y3(1))*0.0102;
    posx3(1:X3(1),1)=r3(1:X3(1));
    posy3(1:Y3(1),1)=r3(X3(1)+1:end);

    X4(1)              = 0.1*N;                 % # of particles on membrane
    Y4(1)              = 0.1*N;                 % # of particles on membrane
    NNx4(1)            = X4(1);
    NNy4(1)            = Y4(1);
    Tx4(1)             = 0.0;
    Ty4(1)             = 0.0;
    nx4(1:X4(1),1)      = 1;                     % activate mem-bound particles
    ny4(1:X4(1),1)      = 1;
    r4 = randperm(ceil(L/(0.0102)),X4(1)+Y4(1))*0.0102;
    posx4(1:X4(1),1)=r4(1:X4(1));
    posy4(1:Y4(1),1)=r4(X4(1)+1:end);

    % Sample concentration at actin filament spatial scale
    %
    [s1,xC1,yC1] = resamplePolarityMolecules(posx1(1:NNx1(1),1),posy1(1:NNy1(1),1),NNx1(1),NNy1(1),L,Na);
    [s2,xC2,yC2] = resamplePolarityMolecules(posx2(1:NNx2(1),1),posy2(1:NNy2(1),1),NNx2(1),NNy2(1),L,Na);
    [s3,xC3,yC3] = resamplePolarityMolecules(posx3(1:NNx3(1),1),posy3(1:NNy3(1),1),NNx3(1),NNy3(1),L,Na);
    [s4,xC4,yC4] = resamplePolarityMolecules(posx4(1:NNx4(1),1),posy4(1:NNy4(1),1),NNx4(1),NNy4(1),L,Na);

    aic1 = a1;
    bic1 = b1;
    aic2 = a2;
    bic2 = b2;
    aic3 = a3;
    bic3 = b3;
    aic4 = a4;
    bic4 = b4;


    % Setup convergence checks for actin quantities
    %
    conv1_1 = zeros(Nt,2);
    conv2_1 = zeros(Nt,2);
    convsup_1 = zeros(Nt,2);
    conv1_2 = zeros(Nt,2);
    conv2_2 = zeros(Nt,2);
    convsup_2 = zeros(Nt,2);
    conv1_3 = zeros(Nt,2);
    conv2_3 = zeros(Nt,2);
    convsup_3 = zeros(Nt,2);
    conv1_4 = zeros(Nt,2);
    conv2_4 = zeros(Nt,2);
    convsup_4 = zeros(Nt,2);

    %amount cells have moved
    xshift1=zeros(1,Nt);
    yshift1=zeros(1,Nt);
    xshift2=zeros(1,Nt);
    yshift2=zeros(1,Nt);
    xshift3=zeros(1,Nt);
    yshift3=zeros(1,Nt);
    xshift4=zeros(1,Nt);
    yshift4=zeros(1,Nt);

    % Set movie making
    %
    if vid==1
        vidObjScatter.FrameRate = 5;
        vidObjScatter.Quality = 75;
        open(vidObjScatter);

        vidObjCells.FrameRate = 5;
        vidObjCells.Quality = 75;
        open(vidObjCells);
    end

    if vid==1
    %Define colors
    colorLength = 50;
    white = [1,1,1];
    darkyellow = [227/256,180/256,76/256];
    yellow2 = [254/256,254/256,98/256];
    pink = [211/256,95/256,183/256];
    darkpink = [141/256,45/256,113/256];
    whiteyellow2 = [linspace(white(1),yellow2(1),colorLength)',linspace(white(2),yellow2(2),colorLength)',linspace(white(3),yellow2(3),colorLength)'];
    yellow2darkyellow = [linspace(yellow2(1),darkyellow(1),colorLength)',linspace(yellow2(2),darkyellow(2),colorLength)',linspace(yellow2(3),darkyellow(3),colorLength)'];
    whitedarkyellow2 = [whiteyellow2;yellow2darkyellow];
    whitepink = [linspace(white(1),pink(1),colorLength)',linspace(white(2),pink(2),colorLength)',linspace(white(3),pink(3),colorLength)'];
    pinkdarkpink = [linspace(pink(1),darkpink(1),colorLength)',linspace(pink(2),darkpink(2),colorLength)',linspace(pink(3),darkpink(3),colorLength)'];
    whitedarkpink = [whitepink;pinkdarkpink];
    linecolor=[0 0 0];
    linewidth=1;


    branchedColor = whitedarkpink;
    bundledColor = whitedarkyellow2;
    branchedColName = 'Pink';
    bundledColName = 'Yellow';


    % Plot the initial condition
    figscatter = figure(ppp);
    subplot(1,4,1); %Cell 1
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

    subplot(1,4,2); %Cell 2
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

    subplot(1,4,3); %Cell 3
    plot(Xa,a3,'-o','color',branchedColor(end,:),'linewidth',3); hold on;
    plot(Xa,b3,'-ok','color',bundledColor(end,:),'linewidth',3);
    plot(s3,xC3,'-.','color',branchedColor(end,:),'linewidth',1);
    plot(s3,yC3,'-.k','color',bundledColor(end,:),'linewidth',1);
    % xlim([0 10]); ylim([0 2]);
    %title('Time = 0');
    set(gca,'fontname','times','fontsize',20); box on;
    lgd = legend('Branched network','Bundled network','Rac','Rho','Location','northeast');
    lgd.NumColumns = 2;
    set(gcf,'color','w');
    title('Cell 3')
    hold off;

    subplot(1,4,4); %Cell 4
    plot(Xa,a4,'-o','color',branchedColor(end,:),'linewidth',3); hold on;
    plot(Xa,b4,'-ok','color',bundledColor(end,:),'linewidth',3);
    plot(s4,xC4,'-.','color',branchedColor(end,:),'linewidth',1);
    plot(s4,yC4,'-.k','color',bundledColor(end,:),'linewidth',1);
    % xlim([0 10]); ylim([0 2]);
    %title('Time = 0');
    set(gca,'fontname','times','fontsize',20); box on;
    lgd = legend('Branched network','Bundled network','Rac','Rho','Location','northeast');
    lgd.NumColumns = 2;
    set(gcf,'color','w');
    title('Cell 4')
    hold off;

    % Define circles
    gapsize=0;
    [th,rad] = meshgrid((0:3.6:360)*pi/180,0.85:0.01:1);
    [Xcol,Ycol] = pol2cart(th,rad);
    Ycol1=Ycol;
    Ycol2=Ycol;
    Ycol3=Ycol;
    Ycol4=Ycol;
    Xcol1=Xcol;
    Xcol2=Xcol;
    Xcol3=Xcol;
    Xcol4=Xcol;
    Ycol2 = Ycol2 - 2*abs(max(max(Ycol2)))-gapsize;
    Ycol3 = Ycol3 - 2*abs(max(max(Ycol3)))-gapsize;
    Ycol4 = Ycol4 - 0*abs(max(max(Ycol4)))-gapsize;
    Xcol3=Xcol3-2;
    Xcol4=Xcol4-2;
    ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1]';
    ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1]';
    ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2]';
    ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2]';
    ZBranch3 = [a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3]';
    ZBund3 = [b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3]';
    ZBranch4 = [a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4]';
    ZBund4 = [b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4]';
    [th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
    [Xsm,Ysm] = pol2cart(th,rad);
    Ysm1=Ysm;
    Ysm2=Ysm;

    allmax = 12;


    % Circles
    % Cell 1
    t=0;
    figcells=figure(ppp+1);
    clf
    hold on;
    plot3(cos(th(1,:))+xshift1(t+1),sin(th(1,:))+yshift1(t+1),ones(1,length(th))*(allmax+1),'color',linecolor,'LineWidth',linewidth)
    surf(Xcol1+xshift1(t+1),Ycol1+yshift1(t+1),ZBranch1,'AlphaData',ZBranch1,'FaceAlpha','interp','FaceColor','interp');
    view(2)
    colormap(branchedColor)
    clim([0,12])
    freezeColors;
    clim([0,12])
    shading interp
    surf(Xcol1+xshift1(t+1),Ycol1+yshift1(t+1),ZBund1,'AlphaData',ZBund1,'FaceAlpha','interp','FaceColor','interp');
    colormap(bundledColor)
    clim([0,12])
    freezeColors;
    shading interp

    % Cell 2
    plot3(cos(th(1,:))+xshift2(t+1),sin(th(1,:))-2+yshift2(t+1),ones(1,length(th))*(allmax+1),'color',linecolor,'LineWidth',linewidth)
    surf(Xcol2+xshift2(t+1),Ycol2+yshift2(t+1),ZBranch2,'AlphaData',ZBranch2,'FaceAlpha','interp','FaceColor','interp');
    view(2)
    colormap(branchedColor)
    freezeColors;
    clim([0,12])
    shading interp
    surf(Xcol2+xshift2(t+1),Ycol2+yshift2(t+1),ZBund2,'AlphaData',ZBund2,'FaceAlpha','interp','FaceColor','interp');
    colormap(bundledColor)
    freezeColors;
    clim([0,12])
    shading interp

    % Cell 3
    plot3(cos(th(1,:))-2+xshift3(t+1),sin(th(1,:))-2+yshift3(t+1),ones(1,length(th))*(allmax+1),'color',linecolor,'LineWidth',linewidth)
    surf(Xcol3+xshift3(t+1),Ycol3+yshift3(t+1),ZBranch3,'AlphaData',ZBranch3,'FaceAlpha','interp','FaceColor','interp');
    view(2)
    colormap(branchedColor)
    freezeColors;
    clim([0,12])
    shading interp
    surf(Xcol3+xshift3(t+1),Ycol3+yshift3(t+1),ZBund3,'AlphaData',ZBund3,'FaceAlpha','interp','FaceColor','interp');
    colormap(bundledColor)
    freezeColors;
    clim([0,12])
    shading interp

    % Cell 4
    plot3(cos(th(1,:))-2+xshift4(t+1),sin(th(1,:))+yshift4(t+1),ones(1,length(th))*(allmax+1),'color',linecolor,'LineWidth',linewidth)
    surf(Xcol4+xshift4(t+1),Ycol4+yshift4(t+1),ZBranch4,'AlphaData',ZBranch4,'FaceAlpha','interp','FaceColor','interp');
    view(2)
    colormap(branchedColor)
    clim([0,12])
    freezeColors;
    cb=colorbar('Location','eastoutside');
    freezeColors(cb);
    cbpos=cb.Position;
    set(cb,'Position',[cbpos(1)+2*cbpos(3),cbpos(2),cbpos(3),cbpos(4)/2])
    set(cb,'TickLabels',[])
    cbpos=cb.Position;
    shading interp
    surf(Xcol4+xshift4(t+1),Ycol4+yshift4(t+1),ZBund4,'AlphaData',ZBund4,'FaceAlpha','interp','FaceColor','interp');
    colormap(bundledColor)
    clim([0,12])
    freezeColors;
    jcb=jicolorbar;
    freezeColors(jcb);
    set(jcb,'Position',[cbpos(1)+cbpos(3),cbpos(2),cbpos(3),cbpos(4)])
    shading interp
    grid off
    xlim([-5,3])
    ylim([-5,3])
    pbaspect([8 8 1])

    hold off;
    box off;
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gcf,'color','w');


    % Plot signal
    if signal==1
        [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
        [Xsig,Ysig] = pol2cart(th,rad);
        hold on;
        scatter(Xsig(sigBound3),Ysig(sigBound3)-2,'black','.')
        hold off;
    end

    ohf = findobj(gcf);
    figaxes = findobj(ohf(1), 'Type', 'axes');
    set(figaxes(1),'Fontsize',15)
    set(figaxes(2),'Fontsize',14)
    camroll(90)

    if vid==1
        scatterframe = getframe(figscatter);
        writeVideo(vidObjScatter,scatterframe);
        cellsframe = getframe(figcells);
        writeVideo(vidObjCells,cellsframe);
    end
    end
    

    %% Run simulation
    %
    tic
    quit_cond = 0;
    cond = 0;
    for t=1:(Nt-1)

        %% Run biochemistry
        [Konx1,Kony1,Kfbx1,Kfby1,Koffx1,Koffy1] = spatialrates(ron,rfb,roff,a1,b1,s1,beta,cond,boundC1_2); % set rates
        [Konx2,Kony2,Kfbx2,Kfby2,Koffx2,Koffy2] = spatialrates(ron,rfb,roff,a2,b2,s2,beta,cond,[boundC2_1,boundC2_3]);
        [Konx3,Kony3,Kfbx3,Kfby3,Koffx3,Koffy3] = spatialrates(ron,rfb,roff,a3,b3,s3,beta,cond,[boundC3_2,boundC3_4]);
        [Konx4,Kony4,Kfbx4,Kfby4,Koffx4,Koffy4] = spatialrates(ron,rfb,roff,a4,b4,s4,beta,cond,boundC4_3);


        % this works
        if signal==1
            steepness = 20;
            Konx3 = (ron*(tanh(steepness*(s3-s3(sigBound3(1)))) - tanh(steepness*(s3-s3(sigBound3(end)))) + 0.2)/2.2)';
            Kony3 = (ron*(2 - tanh(steepness*(s3-s3(sigBound3(1)))) + tanh(steepness*(s3-s3(sigBound3(end)))) + 0.2)/2.2)';
            Kfbx3 = (rfb*(tanh(steepness*(s3-s3(sigBound3(1)))) - tanh(steepness*(s3-s3(sigBound3(end)))) + 0.2)/2.2)';
            Kfby3 = (rfb*(2 - tanh(steepness*(s3-s3(sigBound3(1)))) + tanh(steepness*(s3-s3(sigBound3(end)))) + 0.2)/2.2)';
            Koffx3 = (roff*(2 - tanh(steepness*(s3-s3(sigBound3(1)))) + tanh(steepness*(s3-s3(sigBound3(end)))) + 0.2)/2.2)';
            Koffy3 = (roff*(tanh(steepness*(s3-s3(sigBound3(1)))) - tanh(steepness*(s3-s3(sigBound3(end)))) + 0.2)/2.2)';
        end


        switch pathway
            case 'uncoupled'
                
            case 'alt-racup-rhoup'
                Konx1(boundC1_2)=Konx1(boundC1_2)*1000;
                Kony2(boundC2_1)=Kony2(boundC2_1)*1000;
                Konx2(boundC2_3)=Konx2(boundC2_3)*1000;
                Kony3(boundC3_2)=Kony3(boundC3_2)*1000;
                Konx3(boundC3_4)=Konx3(boundC3_4)*1000;
                Kony4(boundC4_3)=Kony4(boundC4_3)*1000;
                Konx4(boundC4_1)=Konx4(boundC4_1)*1000;
                Kony1(boundC1_4)=Kony1(boundC1_4)*1000;

            case 'alt-racrup-rhoup-forcedependent'
                Konx1(boundC1_2) = Konx1(boundC1_2).*flip(b2(boundC2_1))*1000;
                Kony2(boundC2_1) = Kony2(boundC2_1).*flip(a1(boundC1_2))*1000;
                Konx2(boundC2_3) = Konx2(boundC2_3).*flip(b3(boundC3_2))*1000;
                Kony3(boundC3_2) = Kony3(boundC3_2).*flip(a2(boundC2_3))*1000;
                Konx3(boundC3_4) = Konx3(boundC3_4).*flip(b4(boundC4_3))*1000;
                Kony4(boundC4_3) = Kony4(boundC4_3).*flip(a3(boundC3_4))*1000;
                Konx4(boundC4_1) = Konx4(boundC4_1).*flip(b1(boundC1_4))*1000;
                Kony1(boundC1_4) = Kony1(boundC1_4).*flip(a4(boundC4_1))*1000;

            case 'branched-bundled-promotion'
                kb_ind=2;
                kc_ind=2;
        end


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
                dn2(j) = (rr<((konx2+kfbx2*nnx2/N)*(N/nnx2-1)/a0))*1.0 + (rr>=((konx2+kfbx2*nnx2/N)*(N/nnx2-1)/a0))*(-1.0);
            end

            [mintaux2,minidx2] = min(taux2(1:j));       % find first chemical rxn
            Tx2(rxn_count_x2+1) = Tx2(rxn_count_x2) + mintaux2;
            X2(rxn_count_x2+1) = nnx2 + dn2(minidx2);
            rxn_count_x2 = rxn_count_x2 + 1;
            NNx2(t+1) = X2(rxn_count_x2-1);
        end

        %Cell 3
        if((t-1)*dt<Tx3(rxn_count_x3))
            NNx3(t+1) = X3(rxn_count_x3-1);
        else
            nnx3 = X3(rxn_count_x3);
            taux3 = zeros(nnx3,1);
            dn3 = zeros(nnx3,1);
            r3 = rand(nnx3,1);

            if(nnx3==0)
                sprintf('here 3rac')
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nnx3          % all agents
                konx3 = interp1(s3,Konx3,posx3(j,t));
                koffx3 = interp1(s3,Koffx3,posx3(j,t));
                kfbx3 = interp1(s3,Kfbx3,posx3(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffx3 + (konx3+kfbx3*nnx3/N)*(N/nnx3-1);
                taux3(j) = -log(r3(j))/a0;
                rr = rand(1,1);
                dn3(j) = (rr<((konx3+kfbx3*nnx3/N)*(N/nnx3-1)/a0))*1.0 + (rr>=((konx3+kfbx3*nnx3/N)*(N/nnx3-1)/a0))*(-1.0);
            end

            [mintaux3,minidx3] = min(taux3(1:j));       % find first chemical rxn
            Tx3(rxn_count_x3+1) = Tx3(rxn_count_x3) + mintaux3;
            X3(rxn_count_x3+1) = nnx3 + dn3(minidx3);
            rxn_count_x3 = rxn_count_x3 + 1;
            NNx3(t+1) = X3(rxn_count_x3-1);
        end

        %Cell 4
        if((t-1)*dt<Tx4(rxn_count_x4))
            NNx4(t+1) = X4(rxn_count_x4-1);
        else
            nnx4 = X4(rxn_count_x4);
            taux4 = zeros(nnx4,1);
            dn4 = zeros(nnx4,1);
            r4 = rand(nnx4,1);

            if(nnx4==0)
                sprintf('here 4rac')
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nnx4          % all agents
                konx4 = interp1(s4,Konx4,posx4(j,t));
                koffx4 = interp1(s4,Koffx4,posx4(j,t));
                kfbx4 = interp1(s4,Kfbx4,posx4(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffx4 + (konx4+kfbx4*nnx4/N)*(N/nnx4-1);
                taux4(j) = -log(r4(j))/a0;
                rr = rand(1,1);
                dn4(j) = (rr<((konx4+kfbx4*nnx4/N)*(N/nnx4-1)/a0))*1.0 + (rr>=((konx4+kfbx4*nnx4/N)*(N/nnx4-1)/a0))*(-1.0);
            end

            [mintaux4,minidx4] = min(taux4(1:j));       % find first chemical rxn
            Tx4(rxn_count_x4+1) = Tx4(rxn_count_x4) + mintaux4;
            X4(rxn_count_x4+1) = nnx4 + dn4(minidx4);
            rxn_count_x4 = rxn_count_x4 + 1;
            NNx4(t+1) = X4(rxn_count_x4-1);
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

        %Cell 3
        if((t-1)*dt<Ty3(rxn_count_y3))
            NNy3(t+1) = Y3(rxn_count_y3-1);
        else
            nny3 = Y3(rxn_count_y3);
            tauy3 = zeros(nny3,1);
            dn3 = zeros(nny3,1);
            r3 = rand(nny3,1);

            if(nny3==0)
                sprintf('here 3rho')
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nny3          % all agents
                kony3 = interp1(s3,Kony3,posy3(j,t));
                koffy3 = interp1(s3,Koffy3,posy3(j,t));
                kfby3 = interp1(s3,Kfby3,posy3(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffy3 + (kony3+kfby3*nny3/N)*(N/nny3-1);
                tauy3(j) = -log(r3(j))/a0;
                rr = rand(1,1);
                dn3(j) = (rr<((kony3+kfby3*nny3/N)*(N/nny3-1)/a0))*1.0 + (rr>=((kony3+kfby3*nny3/N)*(N/nny3-1)/a0))*(-1.0);
            end

            [mintauy3,minidy3] = min(tauy3(1:j));       % find first chemical rxn
            Ty3(rxn_count_y3+1) = Ty3(rxn_count_y3) + mintauy3;
            Y3(rxn_count_y3+1) = nny3 + dn3(minidy3);
            rxn_count_y3 = rxn_count_y3 + 1;
            NNy3(t+1) = Y3(rxn_count_y3-1);
        end

        %Cell 4
        if((t-1)*dt<Ty4(rxn_count_y4))
            NNy4(t+1) = Y4(rxn_count_y4-1);
        else
            nny4 = Y4(rxn_count_y4);
            tauy4 = zeros(nny4,1);
            dn4 = zeros(nny4,1);
            r4 = rand(nny4,1);

            if(nny4==0)
                sprintf('here 4rho')
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nny4          % all agents
                kony4 = interp1(s4,Kony4,posy4(j,t));
                koffy4 = interp1(s4,Koffy4,posy4(j,t));
                kfby4 = interp1(s4,Kfby4,posy4(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffy4 + (kony4+kfby4*nny4/N)*(N/nny4-1);
                tauy4(j) = -log(r4(j))/a0;
                rr = rand(1,1);
                dn4(j) = (rr<((kony4+kfby4*nny4/N)*(N/nny4-1)/a0))*1.0 + (rr>=((kony4+kfby4*nny4/N)*(N/nny4-1)/a0))*(-1.0);
            end

            [mintauy4,minidy4] = min(tauy4(1:j));       % find first chemical rxn
            Ty4(rxn_count_y4+1) = Ty4(rxn_count_y4) + mintauy4;
            Y4(rxn_count_y4+1) = nny4 + dn4(minidy4);
            rxn_count_y4 = rxn_count_y4 + 1;
            NNy4(t+1) = Y4(rxn_count_y4-1);
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
        K1_3 = NNx3(t+1);
        K2_3 = NNy3(t+1);
        K1_4 = NNx4(t+1);
        K2_4 = NNy4(t+1);

        % Between reactions, perform Brownian motion with periodic BC
        r1_1 = rand(K1_1,1);    % coin flip
        nx1(1:K1_1,t+1) = 1;
        posx1(1:K1_1,t+1) = posx1(1:K1_1,t) + dx*((r1_1<p)*1.0 + (r1_1>(1-p))*(-1.0));
        r1_2 = rand(K1_2,1);    % coin flip
        nx2(1:K1_2,t+1) = 1;
        posx2(1:K1_2,t+1) = posx2(1:K1_2,t) + dx*((r1_2<p)*1.0 + (r1_2>(1-p))*(-1.0));
        r1_3 = rand(K1_3,1);    % coin flip
        nx3(1:K1_3,t+1) = 1;
        posx3(1:K1_3,t+1) = posx3(1:K1_3,t) + dx*((r1_3<p)*1.0 + (r1_3>(1-p))*(-1.0));
        r1_4 = rand(K1_4,1);    % coin flip
        nx4(1:K1_4,t+1) = 1;
        posx4(1:K1_4,t+1) = posx4(1:K1_4,t) + dx*((r1_4<p)*1.0 + (r1_4>(1-p))*(-1.0));
        r2_1 = rand(K2_1,1);    % coin flip
        ny1(1:K2_1,t+1) = 1;
        posy1(1:K2_1,t+1) = posy1(1:K2_1,t) + dx*((r2_1<p)*1.0 + (r2_1>(1-p))*(-1.0));
        r2_2 = rand(K2_2,1);    % coin flip
        ny2(1:K2_2,t+1) = 1;
        posy2(1:K2_2,t+1) = posy2(1:K2_2,t) + dx*((r2_2<p)*1.0 + (r2_2>(1-p))*(-1.0));
        r2_3 = rand(K2_3,1);    % coin flip
        ny3(1:K2_3,t+1) = 1;
        posy3(1:K2_3,t+1) = posy3(1:K2_3,t) + dx*((r2_3<p)*1.0 + (r2_3>(1-p))*(-1.0));
        r2_4 = rand(K2_4,1);    % coin flip
        ny4(1:K2_4,t+1) = 1;
        posy4(1:K2_4,t+1) = posy4(1:K2_4,t) + dx*((r2_4<p)*1.0 + (r2_4>(1-p))*(-1.0));


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

        % Cell 3
        firstcoll = sum(ismembertol(posx3(1:K1_3,t+1),posy3(1:K2_3,t+1),0.005,'DataScale',1));
        if firstcoll~=0
            % Get indices of collisions
            aa3 = ismembertol(posx3(1:K1_3,t+1),posy3(1:K2_3,t+1),0.005,'DataScale',1);
            list_idx3 = find(aa3~=0);
            bb3 = ismembertol(posy3(1:K2_3,t+1),posx3(1:K1_3,t+1),0.005,'DataScale',1);
            list_idy3 = find(bb3~=0);

            posx3(list_idx3,t+1) = posx3(list_idx3,t);
            posy3(list_idy3,t+1) = posy3(list_idy3,t);
        end

        % Cell 4
        firstcoll = sum(ismembertol(posx4(1:K1_4,t+1),posy4(1:K2_4,t+1),0.005,'DataScale',1));
        if firstcoll~=0
            % Get indices of collisions
            aa4 = ismembertol(posx4(1:K1_4,t+1),posy4(1:K2_4,t+1),0.005,'DataScale',1);
            list_idx4 = find(aa4~=0);
            bb4 = ismembertol(posy4(1:K2_4,t+1),posx4(1:K1_4,t+1),0.005,'DataScale',1);
            list_idy4 = find(bb4~=0);

            posx4(list_idx4,t+1) = posx4(list_idx4,t);
            posy4(list_idy4,t+1) = posy4(list_idy4,t);
        end

        % Enforce periodic boundary conditions
        posx1(1:K1_1,t+1) = posx1(1:K1_1,t+1) + (-L).*(posx1(1:K1_1,t+1)>L) + (L).*(posx1(1:K1_1,t+1)<0.0);
        posy1(1:K2_1,t+1) = posy1(1:K2_1,t+1) + (-L).*(posy1(1:K2_1,t+1)>L) + (L).*(posy1(1:K2_1,t+1)<0.0);
        posx2(1:K1_2,t+1) = posx2(1:K1_2,t+1) + (-L).*(posx2(1:K1_2,t+1)>L) + (L).*(posx2(1:K1_2,t+1)<0.0);
        posy2(1:K2_2,t+1) = posy2(1:K2_2,t+1) + (-L).*(posy2(1:K2_2,t+1)>L) + (L).*(posy2(1:K2_2,t+1)<0.0);
        posx3(1:K1_3,t+1) = posx3(1:K1_3,t+1) + (-L).*(posx3(1:K1_3,t+1)>L) + (L).*(posx3(1:K1_3,t+1)<0.0);
        posy3(1:K2_3,t+1) = posy3(1:K2_3,t+1) + (-L).*(posy3(1:K2_3,t+1)>L) + (L).*(posy3(1:K2_3,t+1)<0.0);
        posx4(1:K1_4,t+1) = posx4(1:K1_4,t+1) + (-L).*(posx4(1:K1_4,t+1)>L) + (L).*(posx4(1:K1_4,t+1)<0.0);
        posy4(1:K2_4,t+1) = posy4(1:K2_4,t+1) + (-L).*(posy4(1:K2_4,t+1)>L) + (L).*(posy4(1:K2_4,t+1)<0.0);

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

        % Find spontaneous association location cell 3
        ss3 = sort(posx3(1:K1_3,t));
        [ijk3] = find(ss3==posx3(minidx3,t),1);
        prevind3 = (ijk3-1)*(ijk3>1) + (K1_3)*(ijk3==1);
        nextind3 = (ijk3+1)*(ijk3<K1_3) + 1*(ijk3==K1_3);
        x2_3 = posx3(minidx3,t)+(ss3(prevind3)-posx3(minidx3,t))/2;
        x1_3 = posx3(minidx3,t)+(ss3(nextind3)-posx3(minidx3,t))/2;
        locx3 = (x2_3-x1_3).*rand(1,1) + x1_3; % random location halfway between the closest left/right particles
        ss3 = sort(posy3(1:K2_3,t));
        [ijk3] = find(ss3==posy3(minidy3,t),1);
        prevind3 = (ijk3-1)*(ijk3>1) + (K2_3)*(ijk3==1);
        nextind3 = (ijk3+1)*(ijk3<K2_3) + 1*(ijk3==K2_3);
        y2_3 = posy3(minidy3,t)+(ss3(prevind3)-posy3(minidy3,t))/2;
        y1_3 = posy3(minidy3,t)+(ss3(nextind3)-posy3(minidy3,t))/2;
        locy3 = (y2_3-y1_3).*rand(1,1) + y1_3; % random location halfway between the closest left/right particles

        ponx3 = ron/(ron+rfb*(N-K1_3));
        pony3 = ron/(ron+rfb*(N-K2_3));

        % Find spontaneous association location cell 4
        ss4 = sort(posx4(1:K1_4,t));
        [ijk4] = find(ss4==posx4(minidx4,t),1);
        prevind4 = (ijk4-1)*(ijk4>1) + (K1_4)*(ijk4==1);
        nextind4 = (ijk4+1)*(ijk4<K1_4) + 1*(ijk4==K1_4);
        x2_4 = posx4(minidx4,t)+(ss4(prevind4)-posx4(minidx4,t))/2;
        x1_4 = posx4(minidx4,t)+(ss4(nextind4)-posx4(minidx4,t))/2;
        locx4 = (x2_4-x1_4).*rand(1,1) + x1_4; % random location halfway between the closest left/right particles
        ss4 = sort(posy4(1:K2_4,t));
        [ijk4] = find(ss4==posy4(minidy4,t),1);
        prevind4 = (ijk4-1)*(ijk4>1) + (K2_4)*(ijk4==1);
        nextind4 = (ijk4+1)*(ijk4<K2_4) + 1*(ijk4==K2_4);
        y2_4 = posy4(minidy4,t)+(ss4(prevind4)-posy4(minidy4,t))/2;
        y1_4 = posy4(minidy4,t)+(ss4(nextind4)-posy4(minidy4,t))/2;
        locy4 = (y2_4-y1_4).*rand(1,1) + y1_4; % random location halfway between the closest left/right particles

        ponx4 = ron/(ron+rfb*(N-K1_4));
        pony4 = ron/(ron+rfb*(N-K2_4));

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
            if numToRemove>0
                boundC1Scaled=(L*boundC1_2/Na);
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
            if numToRemove>0
                boundC2Scaled1=(L*boundC2_1/Na);
                boundC2Scaled2=(L*boundC2_3/Na);
                locRemovey2 = find(abs(posy2(:,t+1)-posx2(K1_2,t+1))<epsilon,numToRemove);
                numFound = length(locRemovey2);
                if ~isempty(locRemovey2) && ((boundC2Scaled1(1)<=posx2(K1_2,t+1) && boundC2Scaled1(end)>=posx2(K1_2,t+1)) ...
                        || (boundC2Scaled2(1)<=posx2(K1_2,t+1) && boundC2Scaled2(end)>=posx2(K1_2,t+1)))
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

        %Cell 3
        if(NNx3(t+1) < NNx3(t))                % diassociation event (particle off)
            oldcol = posx3(minidx3,1:end);
            othercols = posx3([1:minidx3-1,minidx3+1:K1_3],1:end);
            otherothercols = posx3(K1_3+1:end,1:end);
            newpos = [othercols;oldcol;otherothercols];
            posx3 = newpos;
            nx3(K1_3,t+1) = 0;
        elseif(NNx3(t+1) > NNx3(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posx3(K1_3,t+1) = posx3(K1_3,t)+(rr<ponx3)*locx3;              % on event
            posx3(K1_3,t+1) = posx3(K1_3,t)+(rr>=ponx3)*posx3(minidx3,t);   % recruitment event
            nx3(K1_3,t+1) = 1;
            % Look for nearby rho (posy3), take them off
            % locx3=location of rac binding
            if numToRemove>0
                boundC3Scaled1=(L*boundC3_2/Na);
                boundC3Scaled2=(L*boundC3_4/Na);
                locRemovey3 = find(abs(posy3(:,t+1)-posx3(K1_3,t+1))<epsilon,numToRemove);
                numFound = length(locRemovey3);
                if ~isempty(locRemovey3) && ((boundC3Scaled1(1)<=posx3(K1_3,t+1) && boundC3Scaled1(end)>=posx3(K1_3,t+1)) ...
                        || (boundC3Scaled2(1)<=posx3(K1_3,t+1) && boundC3Scaled2(end)>=posx2(K1_3,t+1)))
                    % posy2(locRemovey2,t+1)=0;
                    oldcol = posy3(locRemovey3,1:end); % Find the particle to be removed
                    othercols = posy3(setdiff(1:K2_3,locRemovey3),1:end); % Gather other "on" particles
                    otherothercols = posy3(K2_3+1:end,1:end); % Gather "off" particles
                    newpos = [othercols;oldcol;otherothercols]; % Put removed particle at the end of "on" particles
                    posy3 = newpos;
                    ny3(K2_3-numFound+1:K2_3,t+1) = 0;
                end
            end
        end

        %Cell 4
        if(NNx4(t+1) < NNx4(t))                % diassociation event (particle off)
            oldcol = posx4(minidx4,1:end); % Find the particle to be removed
            othercols = posx4([1:minidx4-1,minidx4+1:K1_4],1:end); % Gather other "on" particles
            otherothercols = posx4(K1_4+1:end,1:end); % Gather "off" particles
            newpos = [othercols;oldcol;otherothercols]; % Put removed particle at the end of "on" particles
            posx4 = newpos;
            nx4(K1_4,t+1) = 0; % Set the removed particle to inactive
        elseif(NNx4(t+1) > NNx4(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posx4(K1_4,t+1) = posx4(K1_4,t)+(rr<ponx4)*locx4; % on event
            posx4(K1_4,t+1) = posx4(K1_4,t)+(rr>=ponx4)*posx4(minidx4,t);   % recruitment event
            nx4(K1_4,t+1) = 1;
            % Look for nearby rho (posy1), take them off
            % posx1(K1_1,t+1)=location of rac binding
            if numToRemove>0
                boundC4Scaled=(L*boundC4_3/Na);
                locRemovey4 = find(abs(posy4(:,t+1)-posx4(K1_4,t+1))<epsilon,numToRemove);
                numFound = length(locRemovey4);
                if ~isempty(locRemovey4) && boundC4Scaled(1)<=posx4(K1_4,t+1) && boundC4Scaled(end)>=posx4(K1_4,t+1)
                    % posy1(locRemovey1,t+1)=0;
                    oldcol = posy4(locRemovey4,1:end); % Find the particle(s) to be removed
                    othercols = posy4(setdiff(1:K2_4,locRemovey4),1:end); % Gather other "on" particles
                    otherothercols = posy4(K2_4+1:end,1:end); % Gather "off" particles
                    newpos = [othercols;oldcol;otherothercols]; % Put removed particle at the end of "on" particles
                    posy4 = newpos;
                    ny4(K2_4-numFound+1:K2_4,t+1) = 0;
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

        %Cell 3
        if (NNy3(t+1) < NNy3(t))                % diassociation event (particle off)
            oldcol = posy3(minidy3,1:end);
            othercols = posy3([1:minidy3-1,minidy3+1:K2_3],1:end);
            otherothercols = posy3(K2_3+1:end,1:end);
            newpos = [othercols;oldcol;otherothercols];
            posy3 = newpos;
            ny3(K2_3,t+1) = 0;
        elseif(NNy3(t+1) > NNy3(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posy3(K2_3,t+1) = posy3(K2_3,t)+(rr<pony3)*locy3;               % on event
            posy3(K2_3,t+1) = posy3(K2_3,t)+(rr>=pony3)*posy3(minidy3,t);    % recruitment event
            ny3(K2_3,t+1) = 1;
        end

        %Cell 4
        if (NNy4(t+1) < NNy4(t))                % diassociation event (particle off)
            oldcol = posy4(minidy4,1:end);
            othercols = posy4([1:minidy4-1,minidy4+1:K2_4],1:end);
            otherothercols = posy4(K2_4+1:end,1:end);
            newpos = [othercols;oldcol;otherothercols];
            posy4 = newpos;
            ny4(K2_4,t+1) = 0;
        elseif(NNy4(t+1) > NNy4(t))             % association event (on or recruitment)
            rr = rand(1,1);
            posy4(K2_4,t+1) = posy4(K2_4,t)+(rr<pony4)*locy4;               % on event
            posy4(K2_4,t+1) = posy4(K2_4,t)+(rr>=pony4)*posy4(minidy4,t);    % recruitment event
            ny4(K2_4,t+1) = 1;
        end

        [s1,xC1,yC1] = resamplePolarityMolecules(posx1(1:K1_1,t+1),posy1(1:K2_1,t+1),K1_1,K2_1,L,Na);
        [s2,xC2,yC2] = resamplePolarityMolecules(posx2(1:K1_2,t+1),posy2(1:K2_2,t+1),K1_2,K2_2,L,Na);
        [s3,xC3,yC3] = resamplePolarityMolecules(posx3(1:K1_3,t+1),posy3(1:K2_3,t+1),K1_3,K2_3,L,Na);
        [s4,xC4,yC4] = resamplePolarityMolecules(posx4(1:K1_4,t+1),posy4(1:K2_4,t+1),K1_4,K2_4,L,Na);

        %% Update actin filaments
        diffRHSa1 = Hm1*a1;
        diffRHSb1 = Hm1*b1;
        diffRHSa2 = Hm2*a2;
        diffRHSb2 = Hm2*b2;
        diffRHSa3 = Hm3*a3;
        diffRHSb3 = Hm3*b3;
        diffRHSa4 = Hm4*a4;
        diffRHSb4 = Hm4*b4;

        %kb = branched pushing on bundled
        %kc = bundled pulling on branched
        cell1_bound2=zeros(length(a1),1);
        cell1_bound2(boundC1_2)=ones(length(boundC1_2),1);
        cell1_bound4=zeros(length(a1),1);
        cell1_bound4(boundC1_4)=ones(length(boundC1_4),1);

        cell2_bound1=zeros(length(a2),1);
        cell2_bound1(boundC2_1)=ones(length(boundC2_1),1);
        cell2_bound3=zeros(length(a2),1);
        cell2_bound3(boundC2_3)=ones(length(boundC2_3),1);

        cell3_bound2=zeros(length(a3),1);
        cell3_bound2(boundC3_2)=ones(length(boundC3_2),1);
        cell3_bound4=zeros(length(a3),1);
        cell3_bound4(boundC3_4)=ones(length(boundC3_4),1);

        cell4_bound3=zeros(length(a4),1);
        cell4_bound3(boundC4_3)=ones(length(boundC4_3),1);
        cell4_bound1=zeros(length(a4),1);
        cell4_bound1(boundC4_1)=ones(length(boundC4_1),1);
        abmax=50;

        gamma=1.5;


        a1_bound4=[a1((Na-1)/2:Na);a1(1:(Na-1)/2-1)];
        b1_bound4=[b1((Na-1)/2:Na);b1(1:(Na-1)/2-1)];
        a2_bound3=[a1((Na-1)/2:Na);a2(1:(Na-1)/2-1)];
        b2_bound3=[b2((Na-1)/2:Na);b2(1:(Na-1)/2-1)];
        a3_bound2=[a3((Na-1)/2:Na);a3(1:(Na-1)/2-1)];
        b3_bound2=[b3((Na-1)/2:Na);b3(1:(Na-1)/2-1)];
        a4_bound1=[a4((Na-1)/2:Na);a4(1:(Na-1)/2-1)];
        b4_bound1=[b4((Na-1)/2:Na);b4(1:(Na-1)/2-1)];

        % Growth term maxes out version
        rxna1 = dt*( F(a1,b1) + (a1.*(1+alpha(1)*xC1 ...
            + ka_vals(ka_ind)*cell1_bound2.* (flip(a2).*(flip(a2)<=abmax) + abmax*(flip(a2)>abmax)) ...
            + kb_vals(kb_ind)*cell1_bound2.* (flip(b2).*(flip(b2)<=abmax) + abmax*(flip(b2)>abmax)) ...
            + ka_vals(ka_ind)*cell1_bound4.* (flip(a4_bound1).*(flip(a4_bound1)<=abmax) + abmax*(flip(a4_bound1)>abmax)) ...
            + kb_vals(kb_ind)*cell1_bound4.* (flip(b4_bound1).*(flip(b4_bound1)<=abmax) + abmax*(flip(b4_bound1)>abmax)) )) - a1.*a1); %Cell 1 branched
        rxnb1 = dt*( F(b1,a1) + (b1.*(1+alpha(1)*yC1 ...
            + kc_vals(kc_ind)*cell1_bound2.* (flip(a2).*(flip(a2)<=abmax) + abmax*(flip(a2)>abmax)) ...
            + kd_vals(kd_ind)*cell1_bound2.* (flip(b2).*(flip(b2)<=abmax) + abmax*(flip(b2)>abmax)) ...
            + kc_vals(kc_ind)*cell1_bound4.* (flip(a4_bound1).*(flip(a4_bound1)<=abmax) + abmax*(flip(a4_bound1)>abmax)) ...
            + kd_vals(kd_ind)*cell1_bound4.* (flip(b4_bound1).*(flip(b4_bound1)<=abmax) + abmax*(flip(b4_bound1)>abmax)) )) - b1.*b1); %Cell 1 bundled
        rxna2 = dt*( F(a2,b2) + (a2.*(1+alpha(1)*xC2 ...
            + ka_vals(ka_ind)*cell2_bound1.* (flip(a1).*(flip(a1)<=abmax) + abmax*(flip(a1)>abmax)) ...
            + kb_vals(kb_ind)*cell2_bound1.* (flip(b1).*(flip(b1)<=abmax) + abmax*(flip(b1)>abmax)) ...
            + ka_vals(ka_ind)*cell2_bound3.* (flip(a3_bound2).*(flip(a3_bound2)<=abmax) + abmax*(flip(a3_bound2)>abmax)) ...
            + kb_vals(kb_ind)*cell2_bound3.* (flip(b3_bound2).*(flip(b3_bound2)<=abmax) + abmax*(flip(b3_bound2)>abmax)) )) - a2.*a2); %Cell 2 branched
        rxnb2 = dt*( F(b2,a2) + (b2.*(1+alpha(1)*yC2 ...
            + kc_vals(kc_ind)*cell2_bound1.* (flip(a1).*(flip(a1)<=abmax) + abmax*(flip(a1)>abmax)) ...
            + kd_vals(kd_ind)*cell2_bound1.* (flip(b1).*(flip(b1)<=abmax) + abmax*(flip(b1)>abmax)) ...
            + kc_vals(kc_ind)*cell2_bound3.* (flip(a3_bound2).*(flip(a3_bound2)<=abmax) + abmax*(flip(a3_bound2)>abmax)) ...
            + kd_vals(kd_ind)*cell2_bound3.* (flip(b3_bound2).*(flip(b3_bound2)<=abmax) + abmax*(flip(b3_bound2)>abmax)) )) - b2.*b2); %Cell 2 bundled
        rxna3 = dt*( F(a3,b3) + (a3.*(1+alpha(1)*xC3 ...
            + ka_vals(ka_ind)*cell3_bound2.* (flip(a2_bound3).*(flip(a2_bound3)<=abmax) + abmax*(flip(a2_bound3)>abmax)) ...
            + kb_vals(kb_ind)*cell3_bound2.* (flip(b2_bound3).*(flip(b2_bound3)<=abmax) + abmax*(flip(b2_bound3)>abmax)) ...
            + ka_vals(ka_ind)*cell3_bound4.* (flip(a4).*(flip(a4)<=abmax) + abmax*(flip(a4)>abmax)) ...
            + kb_vals(kb_ind)*cell3_bound4.* (flip(b4).*(flip(b4)<=abmax) + abmax*(flip(b4)>abmax)) )) - a3.*a3); %Cell 3 branched
        rxnb3 = dt*( F(b3,a3) + (b3.*(1+alpha(1)*yC3 ...
            + kc_vals(kc_ind)*cell3_bound2.* (flip(a2_bound3).*(flip(a2_bound3)<=abmax) + abmax*(flip(a2_bound3)>abmax)) ...
            + kd_vals(kd_ind)*cell3_bound2.* (flip(b2_bound3).*(flip(b2_bound3)<=abmax) + abmax*(flip(b2_bound3)>abmax)) ...
            + kc_vals(kc_ind)*cell3_bound4.* (flip(a4).*(flip(a4)<=abmax) + abmax*(flip(a4)>abmax)) ...
            + kd_vals(kd_ind)*cell3_bound4.* (flip(b4).*(flip(b4)<=abmax) + abmax*(flip(b4)>abmax)) )) - b3.*b3); %Cell 3 bundled
        rxna4 = dt*( F(a4,b4) + (a4.*(1+alpha(1)*xC4 ...
            + ka_vals(ka_ind)*cell4_bound3.* (flip(a3).*(flip(a3)<=abmax) + abmax*(flip(a3)>abmax)) ...
            + kb_vals(kb_ind)*cell4_bound3.* (flip(b3).*(flip(b3)<=abmax) + abmax*(flip(b3)>abmax)) ...
            + ka_vals(ka_ind)*cell4_bound1.* (flip(a1_bound4).*(flip(a1_bound4)<=abmax) + abmax*(flip(a1_bound4)>abmax)) ...
            + kb_vals(kb_ind)*cell4_bound1.* (flip(b1_bound4).*(flip(b1_bound4)<=abmax) + abmax*(flip(b1_bound4)>abmax)) )) - a4.*a4); %Cell 4 branched
        rxnb4 = dt*( F(b4,a4) + (b4.*(1+alpha(1)*yC4 ...
            + kc_vals(kc_ind)*cell4_bound3.* (flip(a3).*(flip(a3)<=abmax) + abmax*(flip(a3)>abmax)) ...
            + kd_vals(kd_ind)*cell4_bound3.* (flip(b3).*(flip(b3)<=abmax) + abmax*(flip(b3)>abmax)) ...
            + kc_vals(kc_ind)*cell4_bound1.* (flip(a1_bound4).*(flip(a1_bound4)<=abmax) + abmax*(flip(a1_bound4)>abmax)) ...
            + kd_vals(kd_ind)*cell4_bound1.* (flip(b1_bound4).*(flip(b1_bound4)<=abmax) + abmax*(flip(b1_bound4)>abmax)) )) - b4.*b4); %Cell 4 bundled

        a1 = Hs1\(diffRHSa1+rxna1);
        b1 = Hs1\(diffRHSb1+rxnb1);
        a2 = Hs2\(diffRHSa2+rxna2);
        b2 = Hs2\(diffRHSb2+rxnb2);
        a3 = Hs3\(diffRHSa3+rxna3);
        b3 = Hs3\(diffRHSb3+rxnb3);
        a4 = Hs4\(diffRHSa4+rxna4);
        b4 = Hs4\(diffRHSb4+rxnb4);

        a1all(:,t)=a1;
        a2all(:,t)=a2;
        b1all(:,t)=b1;
        b2all(:,t)=b2;
        a3all(:,t)=a3;
        a4all(:,t)=a4;
        b3all(:,t)=b3;
        b4all(:,t)=b4;

        %Calculate direction angles
        % cell 1
        a1New = a1;
        a1New(a1New<1)=0;
        if (a1New(1)~=0 && a1New(end)~=0)
            zeroInda1_1=find(a1New==0,1,'first');
            zeroInda2_1=find(a1New==0,1,'last');
            dirIndexa1=ceil((zeroInda1_1+zeroInda2_1)/2) - 50;
        else
            inda1_1=find(a1New~=0,1,'first');
            inda2_1=find(a1New~=0,1,'last');
            dirIndexa1=ceil((inda1_1+inda2_1)/2);
        end
        b1New = b1;
        b1New(b1New<1)=0;
        if (b1New(1)~=0 && b1New(end)~=0)
            zeroIndb1_1=find(b1New==0,1,'first');
            zeroIndb2_1=find(b1New==0,1,'last');
            dirIndexb1=ceil((zeroIndb1_1+zeroIndb2_1)/2) - 50;
        else
            indb1_1=find(b1New~=0,1,'first');
            indb2_1=find(b1New~=0,1,'last');
            dirIndexb1=ceil((indb1_1+indb2_1)/2);
        end
        % cell 2
        a2New = a2;
        a2New(a2New<1)=0;
        if (a2New(1)~=0 && a2New(end)~=0)
            zeroInda1_2=find(a2New==0,1,'first');
            zeroInda2_2=find(a2New==0,1,'last');
            dirIndexa2=ceil((zeroInda1_2+zeroInda2_2)/2) - 50;
        else
            inda1_2=find(a2New~=0,1,'first');
            inda2_2=find(a2New~=0,1,'last');
            dirIndexa2=ceil((inda1_2+inda2_2)/2);
        end
        b2New = b2;
        b2New(b2New<1)=0;
        if (b2New(1)~=0 && b2New(end)~=0)
            zeroIndb1_2=find(b2New==0,1,'first');
            zeroIndb2_2=find(b2New==0,1,'last');
            dirIndexb2=ceil((zeroIndb1_2+zeroIndb2_2)/2) - 50;
        else
            indb1_2=find(b2New~=0,1,'first');
            indb2_2=find(b2New~=0,1,'last');
            dirIndexb2=ceil((indb1_2+indb2_2)/2);
        end
        % cell 3
        a3New = a3;
        a3New(a3New<1)=0;
        if (a3New(1)~=0 && a3New(end)~=0)
            zeroInda1_3=find(a3New==0,1,'first');
            zeroInda2_3=find(a3New==0,1,'last');
            dirIndexa3=ceil((zeroInda1_3+zeroInda2_3)/2) - 50;
        else
            inda1_3=find(a3New~=0,1,'first');
            inda2_3=find(a3New~=0,1,'last');
            dirIndexa3=ceil((inda1_3+inda2_3)/2);
        end
        b3New = b3;
        b3New(b3New<1)=0;
        if (b3New(1)~=0 && b3New(end)~=0)
            zeroIndb1_3=find(b3New==0,1,'first');
            zeroIndb2_3=find(b3New==0,1,'last');
            dirIndexb3=ceil((zeroIndb1_3+zeroIndb2_3)/2) - 50;
        else
            indb1_3=find(b3New~=0,1,'first');
            indb2_3=find(b3New~=0,1,'last');
            dirIndexb3=ceil((indb1_3+indb2_3)/2);
        end
        % cell 4
        a4New = a4;
        a4New(a4New<1)=0;
        if (a4New(1)~=0 && a4New(end)~=0)
            zeroInda1_4=find(a4New==0,1,'first');
            zeroInda2_4=find(a4New==0,1,'last');
            dirIndexa4=ceil((zeroInda1_4+zeroInda2_4)/2) - 50;
        else
            inda1_4=find(a4New~=0,1,'first');
            inda2_4=find(a4New~=0,1,'last');
            dirIndexa4=ceil((inda1_4+inda2_4)/2);
        end
        b4New = b4;
        b4New(b4New<1)=0;
        if (b4New(1)~=0 && b4New(end)~=0)
            zeroIndb1_4=find(b4New==0,1,'first');
            zeroIndb2_4=find(b4New==0,1,'last');
            dirIndexb4=ceil((zeroIndb1_4+zeroIndb2_4)/2) - 50;
        else
            indb1_4=find(b4New~=0,1,'first');
            indb2_4=find(b4New~=0,1,'last');
            dirIndexb4=ceil((indb1_4+indb2_4)/2);
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
        if dirIndexa3<1
            dirIndexa3=dirIndexa3+101;
        end
        if dirIndexb3<1
            dirIndexb3=dirIndexb3+101;
        end
        if dirIndexa4<1
            dirIndexa4=dirIndexa4+101;
        end
        if dirIndexb4<1
            dirIndexb4=dirIndexb4+101;
        end

        %move cells
        [th,rad] = meshgrid((0:3.6:360)*pi/180,1);
        if move_cells==1
            if ~isempty(dirIndexa1) && ~isempty(dirIndexb1)
                xshift1(t+1)=xshift1(t)+cos(th(dirIndexa1))*0.0005;
                yshift1(t+1)=yshift1(t)+sin(th(dirIndexa1))*0.0005;
            end
            if ~isempty(dirIndexa2) && ~isempty(dirIndexb2)
                xshift2(t+1)=xshift2(t)+cos(th(dirIndexa2))*0.0005;
                yshift2(t+1)=yshift2(t)+sin(th(dirIndexa2))*0.0005;
            end
            if ~isempty(dirIndexa3) && ~isempty(dirIndexb3)
                xshift3(t+1)=xshift3(t)+cos(th(dirIndexa3))*0.0005;
                yshift3(t+1)=yshift3(t)+sin(th(dirIndexa3))*0.0005;
            end
            if ~isempty(dirIndexa4) && ~isempty(dirIndexb4)
                xshift4(t+1)=xshift4(t)+cos(th(dirIndexa4))*0.0005;
                yshift4(t+1)=yshift4(t)+sin(th(dirIndexa4))*0.0005;
            end
        end

        %% Plot the solution(s)
        if vid==1
            make_plot=(mod(t,tplot) == 0 || t==Nt-1);
        else
            make_plot=(t==Nt-1);
        end
        if make_plot
            % if t==(Nt-1)

            %Define colors
            colorLength = 50;
            white = [1,1,1];
            darkyellow = [227/256,180/256,76/256];
            yellow2 = [254/256,254/256,98/256];
            pink = [211/256,95/256,183/256];
            darkpink = [141/256,45/256,113/256];

            whiteyellow2 = [linspace(white(1),yellow2(1),colorLength)',linspace(white(2),yellow2(2),colorLength)',linspace(white(3),yellow2(3),colorLength)'];
            yellow2darkyellow = [linspace(yellow2(1),darkyellow(1),colorLength)',linspace(yellow2(2),darkyellow(2),colorLength)',linspace(yellow2(3),darkyellow(3),colorLength)'];
            whitedarkyellow2 = [whiteyellow2;yellow2darkyellow];
            whitepink = [linspace(white(1),pink(1),colorLength)',linspace(white(2),pink(2),colorLength)',linspace(white(3),pink(3),colorLength)'];
            pinkdarkpink = [linspace(pink(1),darkpink(1),colorLength)',linspace(pink(2),darkpink(2),colorLength)',linspace(pink(3),darkpink(3),colorLength)'];
            whitedarkpink = [whitepink;pinkdarkpink];


            branchedColor = whitedarkpink;
            bundledColor = whitedarkyellow2;
            branchedColName = 'Pink';
            bundledColName = 'Yellow';
            linecolor=[0 0 0];
            linewidth=1;


            % Make scatterplots
            scatplot=figure(ppp);
            subplot(1,4,1); %Cell 1
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

            subplot(1,4,2); %Cell 2
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

            subplot(1,4,3); %Cell 3
            plot(Xa,a3,'-o','color',branchedColor(end,:),'linewidth',3); hold on;
            plot(Xa,b3,'-ok','color',bundledColor(end,:),'linewidth',3);
            plot(s3,xC3,'-.','color',branchedColor(end,:),'linewidth',1);
            plot(s3,yC3,'-.k','color',bundledColor(end,:),'linewidth',1);
            % xlim([0 10]); ylim([0 2]);
            %title('Time = 0');
            set(gca,'fontname','times','fontsize',20); box on;
            lgd = legend('Branched network','Bundled network','Rac','Rho','Location','northeast');
            lgd.NumColumns = 2;
            set(gcf,'color','w');
            title('Cell 3')
            hold off;

            subplot(1,4,4); %Cell 4
            plot(Xa,a4,'-o','color',branchedColor(end,:),'linewidth',3); hold on;
            plot(Xa,b4,'-ok','color',bundledColor(end,:),'linewidth',3);
            plot(s4,xC4,'-.','color',branchedColor(end,:),'linewidth',1);
            plot(s4,yC4,'-.k','color',bundledColor(end,:),'linewidth',1);
            % xlim([0 10]); ylim([0 2]);
            %title('Time = 0');
            set(gca,'fontname','times','fontsize',20); box on;
            lgd = legend('Branched network','Bundled network','Rac','Rho','Location','northeast');
            lgd.NumColumns = 2;
            set(gcf,'color','w');
            title('Cell 4')
            hold off;




            % Define circles
            gapsize=0;
            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.85:0.01:1);
            [Xcol,Ycol] = pol2cart(th,rad);
            Ycol1=Ycol;
            Ycol2=Ycol;
            Ycol3=Ycol;
            Ycol4=Ycol;
            Xcol1=Xcol;
            Xcol2=Xcol;
            Xcol3=Xcol;
            Xcol4=Xcol;
            Ycol2 = Ycol2 - 2*abs(max(max(Ycol2)))-gapsize;
            Ycol3 = Ycol3 - 2*abs(max(max(Ycol3)))-gapsize;
            Ycol4 = Ycol4 - 0*abs(max(max(Ycol4)))-gapsize;
            Xcol3=Xcol3-2;
            Xcol4=Xcol4-2;
            ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1]';
            ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1]';
            ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2]';
            ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2]';
            ZBranch3 = [a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3 a3]';
            ZBund3 = [b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3 b3]';
            ZBranch4 = [a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4 a4]';
            ZBund4 = [b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4 b4]';
            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
            [Xsm,Ysm] = pol2cart(th,rad);
            Ysm1=Ysm;
            Ysm2=Ysm;

            allmax = 12;

            % Circles
            % Cell 1
            showtime=1;
            figcells=figure(ppp+1);
            clf
            hold on;
            plot3(cos(th(1,:))+xshift1(t+1),sin(th(1,:))+yshift1(t+1),ones(1,length(th))*(allmax+1),'color',linecolor,'LineWidth',linewidth)
            surf(Xcol1+xshift1(t+1),Ycol1+yshift1(t+1),ZBranch1,'AlphaData',ZBranch1,'FaceAlpha','interp','FaceColor','interp');
            view(2)
            colormap(branchedColor)
            clim([0,12])
            freezeColors;
            clim([0,12])
            shading interp
            surf(Xcol1+xshift1(t+1),Ycol1+yshift1(t+1),ZBund1,'AlphaData',ZBund1,'FaceAlpha','interp','FaceColor','interp');
            colormap(bundledColor)
            clim([0,12])
            freezeColors;
            shading interp

            % Cell 2
            plot3(cos(th(1,:))+xshift2(t+1),sin(th(1,:))-2+yshift2(t+1),ones(1,length(th))*(allmax+1),'color',linecolor,'LineWidth',linewidth)
            surf(Xcol2+xshift2(t+1),Ycol2+yshift2(t+1),ZBranch2,'AlphaData',ZBranch2,'FaceAlpha','interp','FaceColor','interp');
            view(2)
            colormap(branchedColor)
            freezeColors;
            clim([0,12])
            shading interp
            surf(Xcol2+xshift2(t+1),Ycol2+yshift2(t+1),ZBund2,'AlphaData',ZBund2,'FaceAlpha','interp','FaceColor','interp');
            colormap(bundledColor)
            freezeColors;
            clim([0,12])
            shading interp

            % Cell 3
            plot3(cos(th(1,:))-2+xshift3(t+1),sin(th(1,:))-2+yshift3(t+1),ones(1,length(th))*(allmax+1),'color',linecolor,'LineWidth',linewidth)
            surf(Xcol3+xshift3(t+1),Ycol3+yshift3(t+1),ZBranch3,'AlphaData',ZBranch3,'FaceAlpha','interp','FaceColor','interp');
            view(2)
            colormap(branchedColor)
            freezeColors;
            clim([0,12])
            shading interp
            surf(Xcol3+xshift3(t+1),Ycol3+yshift3(t+1),ZBund3,'AlphaData',ZBund3,'FaceAlpha','interp','FaceColor','interp');
            colormap(bundledColor)
            freezeColors;
            clim([0,12])
            shading interp

            % Cell 4
            plot3(cos(th(1,:))-2+xshift4(t+1),sin(th(1,:))+yshift4(t+1),ones(1,length(th))*(allmax+1),'color',linecolor,'LineWidth',linewidth)
            surf(Xcol4+xshift4(t+1),Ycol4+yshift4(t+1),ZBranch4,'AlphaData',ZBranch4,'FaceAlpha','interp','FaceColor','interp');
            view(2)
            colormap(branchedColor)
            clim([0,12])
            freezeColors;
            cb=colorbar('Location','eastoutside');
            freezeColors(cb);
            cbpos=cb.Position;
            set(cb,'Position',[cbpos(1)+2*cbpos(3),cbpos(2),cbpos(3),cbpos(4)/2])
            set(cb,'TickLabels',[])
            cbpos=cb.Position;
            shading interp
            surf(Xcol4+xshift4(t+1),Ycol4+yshift4(t+1),ZBund4,'AlphaData',ZBund4,'FaceAlpha','interp','FaceColor','interp');
            colormap(bundledColor)
            clim([0,12])
            freezeColors;
            jcb=jicolorbar;
            freezeColors(jcb);
            set(jcb,'Position',[cbpos(1)+cbpos(3),cbpos(2),cbpos(3),cbpos(4)])
            shading interp
            grid off
            xlim([-5,3])
            ylim([-5,3])
            pbaspect([8 8 1])

            hold off;
            box off;
            set(gca,'XColor','w')
            set(gca,'YColor','w')
            set(gcf,'color','w');



            % Plot arrows
            figure(ppp+1)
            if ~isempty(dirIndexa1) && ~isempty(dirIndexb1)
                hold on;
                quiver(0+xshift1(t+1),0+yshift1(t+1),Xsm(dirIndexa1),Ysm(dirIndexa1),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',2);
                hold off;
            end
            if ~isempty(dirIndexa2) && ~isempty(dirIndexb2)
                hold on;
                quiver(0+xshift2(t+1),-2+yshift2(t+1),Xsm(dirIndexa2),Ysm(dirIndexa2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',2)
                hold off;
            end
            if ~isempty(dirIndexa3) && ~isempty(dirIndexb3)
                hold on;
                quiver(-2+xshift3(t+1),-2+yshift3(t+1),Xsm(dirIndexa3),Ysm(dirIndexa3),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',2)
                hold off;
            end
            if ~isempty(dirIndexa4) && ~isempty(dirIndexb4)
                hold on;
                quiver(-2+xshift4(t+1),0+yshift4(t+1),Xsm(dirIndexa4),Ysm(dirIndexa4),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',2)
                hold off;
            end

            % Plot signal
            if signal==1
                [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
                [Xsig,Ysig] = pol2cart(th,rad);
                hold on;
                scatter(Xsig(sigBound2),Ysig(sigBound2)-2*abs(max(max(Ycol2)))-gapsize,'black','.')
                hold off;
            end

            if showtime==1
                if t==Nt-1
                    timebox=annotation('textbox', [0.75, cbpos(2), 0.1, 0.05], 'String', "t = " + (Nt)*0.01,'FitBoxToText','on','EdgeColor','none','FontSize',20);
                else
                    timebox=annotation('textbox', [0.75, cbpos(2), 0.1, 0.05], 'String', "t = " + (t)*0.01,'FitBoxToText','on','EdgeColor','none','FontSize',20);
                end
            end

            ohf = findobj(gcf);
            figaxes = findobj(ohf(1), 'Type', 'axes');
            set(figaxes(1),'Fontsize',15)
            set(figaxes(2),'Fontsize',14)
            camroll(90)


            % Add frame to video
            if vid==1
                scatterframe = getframe(scatplot);
                writeVideo(vidObjScatter,scatterframe);
                cellsframe = getframe(figcells);
                writeVideo(vidObjCells,cellsframe);
            end

            if savefigs==1 && t==Nt-1
                savefig(figcells,filenameCells);
                savefig(scatplot,filenameScatter);
            end


        end
    end

    % measure of polarized state (1 if polarized and 0 otherwise)
    %st = 1*( (abs(a(1)-b(end))<1e-3 || abs(a(end)-b(1))<1e-3 ) && abs(a(1)-a(end))>1e-3 && abs(b(1)-b(end))>1e-3 );

    if vid==1
        close(vidObjScatter);
        close(vidObjCells);
    end
    sprintf('Simulation %d done',ppp)
    toc
    if(quit_cond==0)
        

        if save_matfile==1
            save(strcat(mat_location,int2str(ppp),'.mat'),...
                'boundC1_4','boundC1_2','boundC2_1','boundC2_3','boundC3_2','boundC3_4','boundC4_3','boundC4_1','posx1','posx2','posx3','posx4','posy1','posy2','posy3','posy4','NNx1','NNx2','NNx3','NNx4',...
                'NNy1','NNy2','NNy3','NNy4','a1all','a2all','a3all','a4all','b1all','b2all','b3all','b4all','Xa','Xb','s1','s2','s3','s4',...
                'xC1','xC2','xC3','xC4','yC1','yC2','yC3','yC4','xshift1','yshift1','xshift2','yshift2','xshift3','yshift3','xshift4','yshift4','Nt')
        end
        ppp = ppp + 1;
    end
end

