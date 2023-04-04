% Simulate competition for resources for actin filaments 
% and fully coupled to stochastic biochemistry for polarity proteins
%
% Mechanics -- two actin networks: branched (a) and contractile (b)
% Polarity proteins : Rac (X) and Rho (Y)
%
% Last updated: 7/12/2019

clear;
close all;

vid = 0;
%vidObj = VideoWriter('sim2cells_uncoupled4.avi');

%for ppp = 1:5
ppp=1;

% Set actin filament parameters
% 
Da      = 0.1;                  % diffusion coefficient for actin
m0      = 2;                  % competition for actin monomers
%kon     = 0.0001;
%koff    = 0.11;
%kfb     = 0.1;

% Set polarity protein parameters
% 
N       = 200;                  % total number of molecules in the cell (conserved)
ron     = 0.001;                % units: 1/min
rfb     = 1;                    % units: 1/min
roff    = 0.9;                  % units: 1/min
D       = 0.012;                % units: um^2/min

% Set feedback
%
alpha = 0.5;
beta = 1.0;

% Set discretization 
%
L      = 2.0;  
dt     = 0.001;                  % temporal discretization
Na     = 101;                    % number of space steps
dxa    = 1.0/((Na-1)/2);         % spatial discretization
Xa     = dxa*((1:Na)'-(Na+1)/2); % grid points
Xa     = Xa + 1.0;
Ya     = Xa;
dya    = dxa;
pa     = dt*Da/(dxa^2);  

Tend   = 10.0;
Nt     = Tend/dt;
dx     = sqrt(2*D*dt);
tplot  = 500; 

posx1 = zeros(N,Nt);              % array of positions of x1(t)
posy1 = zeros(N,Nt);              % array of positions of y1(t)
posx2 = zeros(N,Nt);              % array of positions of x2(t)
posy2 = zeros(N,Nt);              % array of positions of y2(t)

% Actin reaction term
%
%F = @(U,V) U-U.*U - m0*U.*V; 
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
a1new    = zeros(N,1);
b1       = zeros(N,1);
b1new    = zeros(N,1);

a2       = zeros(N,1);
a2new    = zeros(N,1);
b2       = zeros(N,1);
b2new    = zeros(N,1);

% (1) pulse in middle
if (ictype==1)
    a1   = ones(N,1);
    a1new= ones(N,1);
  
    b1   = 0.5*ones(N,1);        
    b1new= 0.5*ones(N,1);
    
    a2   = ones(N,1);
    a2new= ones(N,1);
  
    b2   = 0.5*ones(N,1);        
    b2new= 0.5*ones(N,1);

% (2) random
elseif (ictype==2)
    a1 = 0.1 + 0.9.*rand(length(Xa),1);
    b1 = 0.1 + 0.9.*rand(length(Xa),1);
    
    a2 = 0.1 + 0.9.*rand(length(Ya),1);
    b2 = 0.1 + 0.9.*rand(length(Ya),1);

% (3) arctangent 
elseif (ictype==3)
    steepness = 25;
    a1 = (tanh(steepness*(Xa-1))+1)/2;
    b1 = (1-tanh(steepness*(Xa-1)))/2;
    
    a2 = (tanh(steepness*(Ya-1))+1)/2;
    b2 = (1-tanh(steepness*(Ya-1)))/2;
   
elseif (ictype==4)
% (4) odd condition #1
    steepness = 25;
    a1 = (1-cos(Xa*pi))/2;
    b1 = (tanh(steepness*(Xa-1.5))+1)/2 + (1-tanh(steepness*(Xa-0.5)))/2;
    
elseif (ictype==5)
% (5) odd condition #2
    mu = 0.8; sigma = 0.1;
    a1 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20);
    mu = 0.8; sigma = 0.1;
    b1 = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20);
end

% Laplacian difference operator with no flux boundary conditions
% Crank-Nicolson operators
II = speye(Na,Na);
Lapdiff = spdiags(ones(Na,1)*[0 1 -2 1 0], [-Na+1, -1, 0 , 1, Na-1], Na, Na);
Lapdiff(1,1) = -1;
Lapdiff(end,end) = -1;
Hm = II+(pa/2)*Lapdiff;
Hs = II-(pa/2)*Lapdiff;

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
X1(1)              = 0.1*N;                 % # of particles on membrane
Y1(1)              = 0.1*N;                 % # of particles on membrane
NNx1(1)            = X1(1);
NNy1(1)            = Y1(1);
Tx1(1)             = 0.0;                   
Ty1(1)             = 0.0;
nx1(1:X1(1),1)      = 1;                    % activate mem-bound particles
ny1(1:X1(1),1)      = 1;
r1 = randperm(ceil(L/(0.0102)),X1(1)+Y1(1))*0.0102;
posx1(1:X1(1),1)=r1(1:X1(1));                         
posy1(1:Y1(1),1)=r1(X1(1)+1:end);

rxn_count_x2       = 1;
rxn_count_y2       = 1;
X2(1)              = 0.1*N;                 % # of particles on membrane
Y2(1)              = 0.1*N;                 % # of particles on membrane
NNx2(1)            = X2(1);
NNy2(1)            = Y2(1);
Tx2(1)             = 0.0;                   
Ty2(1)             = 0.0;
nx2(1:X2(1),1)      = 1;                    % activate mem-bound particles
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

% Set movie making
%
if vid==1
    vidObj.FrameRate = 5;
    vidObj.Quality = 75;
    open(vidObj);
end

% Plot the initial condition
%
if vid == 1
    figure(ppp); 
    plot(Xa,a1,'-o','markerfacecolor',[159 219 229]/255,'linewidth',2); hold on;
    plot(Xa,b1,'-dk','markerfacecolor','k','linewidth',2);
    plot(Ya+2.1,a2,'-o','markerfacecolor',[217 83 25]/255,'linewidth',2); hold on;
    plot(Ya+2.1,b2,'-dk','markerfacecolor','k','linewidth',2);
    plot(s1,xC1,'--','color',[0 0.45 0.75],'linewidth',2);
    plot(s1,yC1,'--k','linewidth',2);
    plot(s2+2.1,xC2,'--','color',[217 83 25]/255,'linewidth',2);
    plot(s2+2.1,yC2,'--k','linewidth',2);
    title('Time = 0.0');
    xlim([0 4]); %ylim([0 35]);
    xticks([0 0.5 1 1.5 2]);
    xticklabels({'0','0.25','0.5','0.75','1'});
    xlabel('Location on cell membrane');
    ylabel('Concentration');
    set(gcf,'color','w');
    set(gca,'fontname','times','fontsize',20); box on;
    hold off;
    pause(1.0);
end

if vid==1
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end

cond = 0; % switch in external gradient cue 

% Run simulation
% 
tic
for t=1:(Nt-1)
    
    %% Run biochemisty
    [Konx1,Kony1,Kfbx1,Kfby1,Koffx1,Koffy1] = spatialrates_2CELLS(ron,rfb,roff,a1,b1,s1,beta,cond); % set rates
    [Konx2,Kony2,Kfbx2,Kfby2,Koffx2,Koffy2] = spatialrates_2CELLS(ron,rfb,roff,a2,b2,s2,beta,cond); % set rates
    
    % Rac cell 1
    if((t-1)*dt<Tx1(rxn_count_x1))
        NNx1(t+1) = X1(rxn_count_x1-1);
    else
       	nx1 = X1(rxn_count_x1);
        taux1 = zeros(nx1,1);
        dn = zeros(nx1,1);
        r = rand(nx1,1);
        a0 = 0.;
        
        for j=1:nx1          % all agents
            konx1 = interp1(s1,Konx1,posx1(j,t));
            koffx1 = interp1(s1,Koffx1,posx1(j,t));
            kfbx1 = interp1(s1,Kfbx1,posx1(j,t));
            % Sample earliest time-to-fire (tau)
            a0 = koffx1 + (konx1+kfbx1*nx1/N)*(N/nx1-1);
            taux1(j) = -log(r(j))/a0;
            rr = rand(1,1);
            dn(j) = (rr<((konx1+kfbx1*nx1/N)*(N/nx1-1)/a0))*1.0 + (rr>=((konx1+kfbx1*nx1/N)*(N/nx1-1)/a0))*(-1.0);
        end
        
        [mintaux1,minidx1] = min(taux1(1:j));       % find first chemical rxn
        Tx1(rxn_count_x1+1) = Tx1(rxn_count_x1) + mintaux1;
        X1(rxn_count_x1+1) = nx1 + dn(minidx1);
        rxn_count_x1 = rxn_count_x1 + 1;
        NNx1(t+1) = X1(rxn_count_x1-1);
    end
    
    % Rac cell 2
    if((t-1)*dt<Tx2(rxn_count_x2))
        NNx2(t+1) = X2(rxn_count_x2-1);
    else
       	nx2 = X2(rxn_count_x2);
        taux2 = zeros(nx2,1);
        dn = zeros(nx2,1);
        r = rand(nx2,1);
        a0 = 0.;
        
        for j=1:nx2          % all agents
            konx2 = interp1(s2,Konx2,posx2(j,t));
            koffx2 = interp1(s2,Koffx2,posx2(j,t));
            kfbx2 = interp1(s2,Kfbx2,posx2(j,t));
            % Sample earliest time-to-fire (tau)
            a0 = koffx2 + (konx2+kfbx2*nx2/N)*(N/nx2-1);
            taux2(j) = -log(r(j))/a0;
            rr = rand(1,1);
            dn(j) = (rr<((konx2+kfbx2*nx2/N)*(N/nx2-1)/a0))*1.0 + (rr>=((konx2+kfbx2*nx2/N)*(N/nx2-1)/a0))*(-1.0);
        end
        
        [mintaux2,minidx2] = min(taux2(1:j));       % find first chemical rxn
        Tx2(rxn_count_x2+1) = Tx2(rxn_count_x2) + mintaux2;
        X2(rxn_count_x2+1) = nx2 + dn(minidx2);
        rxn_count_x2 = rxn_count_x2 + 1;
        NNx2(t+1) = X2(rxn_count_x2-1);
    end
    
    % Rho cell 1
    if((t-1)*dt<Ty1(rxn_count_y1))
        NNy1(t+1) = Y1(rxn_count_y1-1);
    else
        ny1 = Y1(rxn_count_y1);
        tauy1 = zeros(ny1,1);
        dn = zeros(ny1,1);
        r = rand(ny1,1);
        a0 = 0.;
        
        for j=1:ny1          % all agents
            kony1 = interp1(s1,Kony1,posy1(j,t));
            koffy1 = interp1(s1,Koffy1,posy1(j,t));
            kfby1 = interp1(s1,Kfby1,posy1(j,t));
            % Sample earliest time-to-fire (tau)
            a0 = koffy1 + (kony1+kfby1*ny1/N)*(N/ny1-1);
            tauy1(j) = -log(r(j))/a0;
            rr = rand(1,1);
            dn(j) = (rr<((kony1+kfby1*ny1/N)*(N/ny1-1)/a0))*1.0 + (rr>=((kony1+kfby1*ny1/N)*(N/ny1-1)/a0))*(-1.0);
        end
        
        [mintauy1,minidy1] = min(tauy1(1:j));       % find first chemical rxn      
        Ty1(rxn_count_y1+1) = Ty1(rxn_count_y1) + mintauy1;
        Y1(rxn_count_y1+1) = ny1 + dn(minidy1);
        rxn_count_y1 = rxn_count_y1 + 1;
        NNy1(t+1) = Y1(rxn_count_y1-1);
    end
    
    % Rho cell 2
    if((t-1)*dt<Ty2(rxn_count_y2))
        NNy2(t+1) = Y2(rxn_count_y2-1);
    else
        ny2 = Y2(rxn_count_y2);
        tauy2 = zeros(ny2,1);
        dn = zeros(ny2,1);
        r = rand(ny2,1);
        a0 = 0.;
        
        for j=1:ny2          % all agents
            kony2 = interp1(s2,Kony2,posy2(j,t));
            koffy2 = interp1(s2,Koffy2,posy2(j,t));
            kfby2 = interp1(s2,Kfby2,posy2(j,t));
            % Sample earliest time-to-fire (tau)
            a0 = koffy2 + (kony2+kfby2*ny2/N)*(N/ny2-1);
            tauy2(j) = -log(r(j))/a0;
            rr = rand(1,1);
            dn(j) = (rr<((kony2+kfby2*ny2/N)*(N/ny2-1)/a0))*1.0 + (rr>=((kony2+kfby2*ny2/N)*(N/ny2-1)/a0))*(-1.0);
        end
        
        [mintauy2,minidy2] = min(tauy2(1:j));       % find first chemical rxn      
        Ty2(rxn_count_y2+1) = Ty2(rxn_count_y2) + mintauy2;
        Y2(rxn_count_y2+1) = ny2 + dn(minidy2);
        rxn_count_y2 = rxn_count_y2 + 1;
        NNy2(t+1) = Y2(rxn_count_y2-1);
    end

    %% Run diffusion of membrane-bound polarity proteins
    p  = 0.5;                  % probability of hoping left or right

    % Fetch the number of particles at this time
    Kx1 = NNx1(t+1);
    Ky1 = NNy1(t+1);
    Kx2 = NNx2(t+1);
    Ky2 = NNy2(t+1);

    % Between reactions, perform Brownian motion with periodic BC
    % cell 1
    r = rand(Kx1,1);    % coin flip
    nx1(1:Kx1,t+1) = 1;   
    posx1(1:Kx1,t+1) = posx1(1:Kx1,t) + dx*((r<p)*1.0 + (r>(1-p))*(-1.0));

    r = rand(Ky1,1);    % coin flip
    ny1(1:Ky1,t+1) = 1;
    posy1(1:Ky1,t+1) = posy1(1:Ky1,t) + dx*((r<p)*1.0 + (r>(1-p))*(-1.0));
    
    % cell 2
    r = rand(Kx2,1);    % coin flip
    nx2(1:Kx2,t+1) = 1;   
    posx2(1:Kx2,t+1) = posx2(1:Kx2,t) + dx*((r<p)*1.0 + (r>(1-p))*(-1.0));

    r = rand(Ky2,1);    % coin flip
    ny2(1:Ky2,t+1) = 1;
    posy2(1:Ky2,t+1) = posy2(1:Ky2,t) + dx*((r<p)*1.0 + (r>(1-p))*(-1.0));

    % Check for collision(s) and resolve any collisions
    %
    % Resolution: winner moves, loser stays goes back to previous timestep
    % (did not work)
    % New resolution: no one advances
    % (works)
    % cell 1
    firstcoll1 = sum(ismembertol(posx1(1:Kx1,t+1),posy1(1:Ky1,t+1),0.005,'DataScale',1));  
    if firstcoll1~=0           
        %sprintf('collision occured @ t=%.2f\n',t+1)
        % Get indices of collisions
        aa = ismembertol(posx1(1:Kx1,t+1),posy1(1:Ky1,t+1),0.005,'DataScale',1);
        list_idx = find(aa~=0);
        bb = ismembertol(posy1(1:Ky1,t+1),posx1(1:Kx1,t+1),0.005,'DataScale',1);
        list_idy = find(bb~=0);

        posx1(list_idx,t+1) = posx1(list_idx,t);
        posy1(list_idy,t+1) = posy1(list_idy,t);

    %         i = firstcoll;
    %         while (i>0) 
    %             r = rand(1,1); 
    %             posx(list_idx(i),t+1) = (r<0.5)*posx(list_idx(i),t+1) + (r>0.5)*(posx(list_idx(i),t)-dx);
    %             posy(list_idy(i),t+1) = (r<0.5)*(posy(list_idy(i),t)-dx) + (r>0.5)*posy(list_idy(i),t+1);  
    %             i = i-1;
    %         end
    end

    % cell 2
    clear aa bb list_idx list_idy
    firstcoll2 = sum(ismembertol(posx2(1:Kx2,t+1),posy2(1:Ky2,t+1),0.005,'DataScale',1));  
    if firstcoll2~=0           
        %sprintf('collision occured @ t=%.2f\n',t+1)
        % Get indices of collisions
        aa = ismembertol(posx2(1:Kx2,t+1),posy2(1:Ky2,t+1),0.005,'DataScale',1);
        list_idx = find(aa~=0);
        bb = ismembertol(posy2(1:Ky2,t+1),posx2(1:Kx2,t+1),0.005,'DataScale',1);
        list_idy = find(bb~=0);

        posx2(list_idx,t+1) = posx2(list_idx,t);
        posy2(list_idy,t+1) = posy2(list_idy,t);
    end
    
    % Enforce periodic boundary conditions
    %posx(1:K1,t+1) = posx(1:K1,t+1) + (-L).*(posx(1:K1,t+1)>L) + (L).*(posx(1:K1,t+1)<0.0);
    %posy(1:K2,t+1) = posy(1:K2,t+1) + (-L).*(posy(1:K2,t+1)>L) + (L).*(posy(1:K2,t+1)<0.0);

    % Enforce no-flux boundary conditions
    posx1(1:Kx1,t+1) = posx1(1:Kx1,t+1) + (posx1(1:Kx1,t)-posx1(1:Kx1,t+1)).*(posx1(1:Kx1,t+1)>=L) + (posx1(1:Kx1,t)-posx1(1:Kx1,t+1)).*(posx1(1:Kx1,t+1)<=0);
    posy1(1:Ky1,t+1) = posy1(1:Ky1,t+1) + (posy1(1:Ky1,t)-posy1(1:Ky1,t+1)).*(posy1(1:Ky1,t+1)>=L) + (posy1(1:Ky1,t)-posy1(1:Ky1,t+1)).*(posy1(1:Ky1,t+1)<=0);
    posx2(1:Kx2,t+1) = posx2(1:Kx2,t+1) + (posx2(1:Kx2,t)-posx2(1:Kx2,t+1)).*(posx2(1:Kx2,t+1)>=L) + (posx2(1:Kx2,t)-posx2(1:Kx2,t+1)).*(posx2(1:Kx2,t+1)<=0);
    posy2(1:Ky2,t+1) = posy2(1:Ky2,t+1) + (posy2(1:Ky2,t)-posy2(1:Ky2,t+1)).*(posy2(1:Ky2,t+1)>=L) + (posy2(1:Ky2,t)-posy2(1:Ky2,t+1)).*(posy2(1:Ky2,t+1)<=0);
    
    % Determine if a biochemical rxn has occured - update positions
    locx1 = posy1(1:Ky1,t)+(1+4*rand(1,1))*dx; % randomly chosen loc not Y1(t)
    locy1 = posx1(1:Kx1,t)+(1+4*rand(1,1))*dx; % randomly chosen loc not X1(t)
    %ponx1 = konx1/(konx1+kfbx1*(N-Kx1));
    %pony1 = kony1/(kony1+kfby1*(N-Ky1));
    ponx1 = ron/(ron+rfb*(N-Kx1));
    pony1 = ron/(ron+rfb*(N-Ky1));
    
    locx2 = posy2(1:Ky2,t)+(1+4*rand(1,1))*dx; % randomly chosen loc not Y2(t)
    locy2 = posx2(1:Kx2,t)+(1+4*rand(1,1))*dx; % randomly chosen loc not X2(t)
    %ponx2 = konx2/(konx2+kfbx2*(N-Kx2));
    %pony2 = kony2/(kony2+kfby2*(N-Ky2));
    ponx2 = ron/(ron+rfb*(N-Kx2));
    pony2 = ron/(ron+rfb*(N-Ky2));
    
    % Rac cell 1
    if(NNx1(t+1) < NNx1(t))                % diassociation event (particle off)
        oldcol = posx1(minidx1,1:end);
        othercols = posx1([1:minidx1-1,minidx1+1:Kx1],1:end);
        otherothercols = posx1(Kx1+1:end,1:end);
        newpos = [othercols;oldcol;otherothercols];
        posx1 = newpos;
        nx1(Kx1,t+1) = 0;
    elseif(NNx1(t+1) > NNx1(t))             % association event (on or recruitment)
        rr = rand(1,1); 
        posx1(Kx1,t+1) = posx1(Kx1,t)+(rr<ponx1)*locx1(randi(Ky1,1));   % on event
        posx1(Kx1,t+1) = posx1(Kx1,t)+(rr>=ponx1)*posx1(minidx1,t);   % recruitment event
        nx1(Kx1+1,t+1) = 1;
    end

    % Rho cell 1
    if (NNy1(t+1) < NNy1(t))                % diassociation event (particle off)
        oldcol = posy1(minidy1,1:end);
        othercols = posy1([1:minidy1-1,minidy1+1:Ky1],1:end);
        otherothercols = posy1(Ky1+1:end,1:end);
        newpos = [othercols;oldcol;otherothercols];
        posy1 = newpos;
        ny1(Ky1,t+1) = 0;
    elseif(NNy1(t+1) > NNy1(t))             % association event (on or recruitment)
        rr = rand(1,1); 
        posy1(Ky1,t+1) = posy1(Ky1,t)+(rr<pony1)*locy1(randi(Kx1,1));   % on event
        posy1(Ky1,t+1) = posy1(Ky1,t)+(rr>=pony1)*posy1(minidy1,t);   % recruitment event
        ny1(Ky1,t+1) = 1;
    end
    
    % Rac cell 2
    if(NNx2(t+1) < NNx2(t))                % diassociation event (particle off)
        oldcol = posx2(minidx2,1:end);
        othercols = posx2([1:minidx2-1,minidx2+1:Kx2],1:end);
        otherothercols = posx1(Kx2+1:end,1:end);
        newpos = [othercols;oldcol;otherothercols];
        posx2 = newpos;
        nx2(Kx2,t+1) = 0;
    elseif(NNx2(t+1) > NNx2(t))             % association event (on or recruitment)
        rr = rand(1,1); 
        posx2(Kx2,t+1) = posx2(Kx2,t)+(rr<ponx2)*locx2(randi(Ky2,1));   % on event
        posx2(Kx2,t+1) = posx2(Kx2,t)+(rr>=ponx2)*posx2(minidx2,t);   % recruitment event
        nx2(Kx2+1,t+1) = 1;
    end

    % Rho cell 2
    if (NNy2(t+1) < NNy2(t))                % diassociation event (particle off)
        oldcol = posy2(minidy2,1:end);
        othercols = posy2([1:minidy2-1,minidy2+1:Ky2],1:end);
        otherothercols = posy2(Ky2+1:end,1:end);
        newpos = [othercols;oldcol;otherothercols];
        posy2 = newpos;
        ny2(Ky2,t+1) = 0;
    elseif(NNy2(t+1) > NNy2(t))             % association event (on or recruitment)
        rr = rand(1,1); 
        posy2(Ky2,t+1) = posy2(Ky2,t)+(rr<pony2)*locy2(randi(Kx2,1));   % on event
        posy2(Ky2,t+1) = posy2(Ky2,t)+(rr>=pony2)*posy2(minidy2,t);   % recruitment event
        ny2(Ky2,t+1) = 1;
    end
    
    [s1,xC1,yC1] = resamplePolarityMolecules(posx1(1:Kx1,t+1),posy1(1:Ky1,t+1),Kx1,Ky1,L,Na);
    [s2,xC2,yC2] = resamplePolarityMolecules(posx2(1:Kx2,t+1),posy2(1:Ky2,t+1),Kx2,Ky2,L,Na);

    
    %% Update actin filaments
    diffRHSa1 = Hm*a1;
    diffRHSb1 = Hm*b1;
    diffRHSa2 = Hm*a2;
    diffRHSb2 = Hm*b2;
  
    rxna1 = dt*( F(a1,b1) + a1.*(1+alpha*xC1));
    rxnb1 = dt*( F(b1,a1) + b1.*(1+alpha*yC1));
    rxna2 = dt*( F(a2,b2) + a2.*(1+alpha*xC2));
    rxnb2 = dt*( F(b2,a2) + b2.*(1+alpha*yC2));

    a1 = Hs\(diffRHSa1+rxna1);
    b1 = Hs\(diffRHSb1+rxnb1);
    a2 = Hs\(diffRHSa2+rxna2);
    b2 = Hs\(diffRHSb2+rxnb2);

%   % Do a local perturbation analysis of polarized solution
%     if(t==5000) 
%         a1 = a;
%         b1 = b;
%         
%         s = zeros(size(a)); shift = 10;
%         %s(1:end+shift) = a(-shift+1:end);
%         %s((end+shift+1):end) = 1;
%         s(shift+1:end) = a(1:end-shift);
%         s(1:shift) = s(shift+1);
%         a = s;
%         s = zeros(size(b)); shift = 10;
%         s(shift+1:end) = b(1:end-shift);
%         s(1:shift) = s(shift+1);
%         %s(1:end+shift) = b(-shift+1:end);
%         b = s;
%     end
  
    %% Plot the solution(s)
    %if mod(t,tplot) == 0
    if t==(Nt-1)
        figure(ppp);
        plot(Xa,a1,'-o','markerfacecolor',[159 219 229]/255,'linewidth',2); hold on;
        plot(Xa,b1,'-dk','markerfacecolor','k','linewidth',2);
        plot(Ya+2.1,a2,'-o','markerfacecolor',[217 83 25]/255,'linewidth',2); hold on;
        plot(Ya+2.1,b2,'-dk','markerfacecolor','k','linewidth',2);
        plot(s1,xC1,'--','color',[0 0.45 0.75],'linewidth',2);
        plot(s1,yC1,'--k','linewidth',2);
        plot(s2+2.1,xC2,'--','color',[217 83 25]/255,'linewidth',2);
        plot(s2+2.1,yC2,'--k','linewidth',2);
    
        title(['Time = ',num2str(t*dt)]);
        xlim([0 4]); %ylim([0 1.0]);
        set(gca,'fontname','times','fontsize',20); box on;
        %legend('Branched network','Contractile network','Rac','Rho','Location','northeast');
        %legend('Branched network (-- Rac)','Bundled network (-- Rho)','Location','northeast');
        xticks([0 0.5 1 1.5 2]);
        xticklabels({'0','0.25','0.5','0.75','1'});
        xlabel('Location on cell membrane');
        ylabel('Concentration');
        set(gcf,'color','w');
        drawnow;
        hold off;
        
        if vid==1
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame);
        end
    end 
end

% measure of polarized state (1 if polarized and 0 otherwise)
%st = 1*( (abs(a(1)-b(end))<1e-3 || abs(a(end)-b(1))<1e-3 ) && abs(a(1)-a(end))>1e-3 && abs(b(1)-b(end))>1e-3 );     

if vid==1    
    close(vidObj);
end
%sprintf('Simulation %d done',ppp)
toc
%end

%% Plot all particle trajectories

% ccx = [255 219 88]/256.*ones(Nt,1);     % mustard yellow
% ccy = [0 0 255]/256.*ones(Nt,1);        % blue
% figure(11); hold on;
% subplot(2,1,2);
% for j=1:10:max(max([NNx1,NNy1]))
%     hold on;
%     %plot(linspace(0,Tend,Nt),pos(j,:));
%     scatter(linspace(0,Tend,Nt),posx1(j,:),1,ccx);
%     scatter(linspace(0,Tend,Nt),posy1(j,:),1,ccy);
%     box on;
%     set(gca,'Color','k','fontsize',30);
%     set(gcf,'color','w');
%     %xlabel('Time');
%     %ylabel('Cell 1 membrane location');
% end
% 
% subplot(2,1,1);
% for j=1:10:max(max([NNx2,NNy2]))
%     hold on;
%     %plot(linspace(0,Tend,Nt),pos(j,:));
%     scatter(linspace(0,Tend,Nt),posx2(j,:),1,ccx);
%     scatter(linspace(0,Tend,Nt),posy2(j,:),1,ccy);
%     box on;
%     set(gca,'Color','k','fontsize',30);
%     set(gcf,'color','w');
%     %xlabel('Time');
%     %ylabel('Cell 2 membrane location');
% end

% Check convergence
% 
% figure(3)
% plot(linspace(0,Tend,Nt),conv1(:,1))
% hold on;
% plot(linspace(0,Tend,Nt),conv2(:,1))
% plot(linspace(0,Tend,Nt),convsup(:,1))