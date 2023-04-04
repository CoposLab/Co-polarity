% Stochastic Altschuler model for polarization along cell membrane
% with 2 molecule types (i.e. Rac and Rho)
% Chemical reaction rates are spatially-dependent
%
% Diffusion of X(t) and Y(t) particles on the cell membrane with exclusion
%
% Three types of chemical reactions can happen on the membrane:
% (1) association: a particle will randomly appear on the membrane 
% (2) disassociation: a particle will randomly dissapear from the membrane
% (3) association thru recruitment: a particle will recruit another
% particle to the membrane
% 
% Last updated: 5/09/2019

clear
close all;

% Set parameters
% 
N       = 200;                  % total number of molecules in the cell (conserved)
kon     = 0.001;                % units: 1/min
kfb     = 1.0;                  % units: 1/min
koff    = 0.9;                  % units: 1/min
D       = 0.012;                % units: um^2/min

% Set discretization
% 
Tend = 10.0;                    % end simulation time
p    = 0.5;                     % probability of hoping left or right
dt   = 0.001;                   % temporal discretization (units: min)
Nt   = Tend/dt;                 % total time steps
dx   = sqrt(2*D*dt);            % spatial discretization (units: um)
L    = 2.0;                     % length of membrane box (units: um)

% Initialize book-keeping vectors
%
MAX_OUTPUT_LENGTH = 10000;
posx1 = zeros(N,Nt);                      % array of positions of X(t)
posy1 = zeros(N,Nt);                      % array of positions of Y(t)
posx2 = zeros(N,Nt);                      % array of positions of X(t)
posy2 = zeros(N,Nt);                      % array of positions of Y(t)
nx1   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
ny1   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
nx2   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
ny2   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
NNx1  = zeros(Nt,1);
NNy1  = zeros(Nt,1);
NNx2  = zeros(Nt,1);
NNy2  = zeros(Nt,1);
Tx1 = zeros(MAX_OUTPUT_LENGTH,1);         % times of chemical reactions for X(t)
Ty1 = zeros(MAX_OUTPUT_LENGTH,1);         % times of chemical reactions for Y(t)
X1  = zeros(MAX_OUTPUT_LENGTH,1);         % number of X(t) molecules on the membrane    
Y1  = zeros(MAX_OUTPUT_LENGTH,1);         % number of Y(t) molecules on the membrane
Tx2 = zeros(MAX_OUTPUT_LENGTH,1);         % times of chemical reactions for X(t)
Ty2 = zeros(MAX_OUTPUT_LENGTH,1);         % times of chemical reactions for Y(t)
X2  = zeros(MAX_OUTPUT_LENGTH,1);         % number of X(t) molecules on the membrane    
Y2  = zeros(MAX_OUTPUT_LENGTH,1);         % number of Y(t) molecules on the membrane

% Set initial conditions (T=0)
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
r1 = randperm(ceil(L/(0.0102)),X1(1)+Y1(1))*0.01;%*0.0102;
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
r2 = randperm(ceil(L/(0.0102)),X2(1)+Y2(1))*0.01;%*0.0102;
posx2(1:X2(1),1)=r2(1:X2(1));                         
posy2(1:Y2(1),1)=r2(X2(1)+1:end);

Ns = L/dx;
ds = dx;
s1  = ds*((1:Ns)'-(Ns+1)/2) + 1.0;
s2  = ds*((1:Ns)'-(Ns+1)/2) + 1.0;
steepness = 10;

%% Spatial rates
assoc = 'coupled';
fb = 'off';
disassoc = 'off';

%% Loop in time
%
tic
for t=1:(Nt-1)
    %% Run biochemistry
    
    % association rate
    switch assoc
        case {'on'}
            Konx1 = 2*kon*((tanh(steepness*(s1-1))+1)/2 + 0.1)/1.1;
            Kony1 = 2*kon*((1-tanh(steepness*(s1-1)))/2 + 0.1)/1.1;
            Konx2 = 2*kon*((tanh(steepness*(s2-1))+1)/2 + 0.1)/1.1;
            Kony2 = 2*kon*((1-tanh(steepness*(s2-1)))/2 + 0.1)/1.1;
        case {'off'}
            Konx1 = kon*ones(length(s1),1);
            Kony1 = kon*ones(length(s1),1);
            Konx2 = kon*ones(length(s2),1);
            Kony2 = kon*ones(length(s2),1);
        case {'coupled'}
            Konx1 = kon*ones(length(s1),1);
            Kony1 = kon*ones(length(s1),1);
            Konx2 = kon*ones(length(s2),1);
            Kony2 = kon*ones(length(s2),1);

            % Rac in cell 2 inhibits Rho in cell 1
            numx2 = nnz(posx2(:,t));        % # of Rac at rear
            ratiox = sum(posx2(1:numx2,t)<=s2(20))/sum(posx2(1:numx2,t));    
            Kony1(end-20:end) = kon/ratiox;

            % Rho in cell 1 inhibits Rac in cell 2
            numy1 = nnz(posy1(:,t));        % # of Rac at rear
            ratioy = sum(posy1(1:numy1,t)>=s1(end-20))/sum(posy1(1:numy1,t));    
            Konx2(1:20) = kon/ratioy;
    end

    % feedback rate
    switch fb
        case {'on'}
            Kfbx1 = 2*kfb*((tanh(steepness*(s1-1))+1)/2 + 0.1)/1.1;
            Kfby1 = 2*kfb*((1-tanh(steepness*(s1-1)))/2 + 0.1)/1.1;
            Kfbx2 = 2*kfb*((tanh(steepness*(s2-1))+1)/2 + 0.1)/1.1;
            Kfby2 = 2*kfb*((1-tanh(steepness*(s2-1)))/2 + 0.1)/1.1;
        case {'off'}
            Kfbx1 = kfb*ones(length(s1),1);
            Kfby1 = kfb*ones(length(s1),1);
            Kfbx2 = kfb*ones(length(s2),1);
            Kfby2 = kfb*ones(length(s2),1);
    end

    % disassociation rate
    switch disassoc
        case {'on'}
            Koffx1 = 2*koff*((1-tanh(steepness*(s1-1)))/2 + 0.1)/1.1;
            Koffy1 = 2*koff*((tanh(steepness*(s1-1))+1)/2 + 0.1)/1.1;
            Koffx2 = 2*koff*((1-tanh(steepness*(s2-1)))/2 + 0.1)/1.1;
            Koffy2 = 2*koff*((tanh(steepness*(s2-1))+1)/2 + 0.1)/1.1;
        case {'off'}
            Koffx1 = koff*ones(length(s1),1);
            Koffy1 = koff*ones(length(s1),1);
            Koffx2 = koff*ones(length(s2),1);
            Koffy2 = koff*ones(length(s2),1);
    end
    
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
    
    % Between reactions, perform Brownian motion
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
     
    % Enforce no-flux boundary conditions
    posx1(1:Kx1,t+1) = posx1(1:Kx1,t+1) + (posx1(1:Kx1,t)-posx1(1:Kx1,t+1)).*(posx1(1:Kx1,t+1)>=L) + (posx1(1:Kx1,t)-posx1(1:Kx1,t+1)).*(posx1(1:Kx1,t+1)<=0);
    posy1(1:Ky1,t+1) = posy1(1:Ky1,t+1) + (posy1(1:Ky1,t)-posy1(1:Ky1,t+1)).*(posy1(1:Ky1,t+1)>=L) + (posy1(1:Ky1,t)-posy1(1:Ky1,t+1)).*(posy1(1:Ky1,t+1)<=0);
    posx2(1:Kx2,t+1) = posx2(1:Kx2,t+1) + (posx2(1:Kx2,t)-posx2(1:Kx2,t+1)).*(posx2(1:Kx2,t+1)>=L) + (posx2(1:Kx2,t)-posx2(1:Kx2,t+1)).*(posx2(1:Kx2,t+1)<=0);
    posy2(1:Ky2,t+1) = posy2(1:Ky2,t+1) + (posy2(1:Ky2,t)-posy2(1:Ky2,t+1)).*(posy2(1:Ky2,t+1)>=L) + (posy2(1:Ky2,t)-posy2(1:Ky2,t+1)).*(posy2(1:Ky2,t+1)<=0);
    
    locx1 = posy1(1:Ky1,t)+(1+4*rand(1,1))*dx; % randomly chosen loc not Y1(t)
    locy1 = posx1(1:Kx1,t)+(1+4*rand(1,1))*dx; % randomly chosen loc not X1(t)
    %ponx1 = konx1/(konx1+kfbx1*(N-Kx1));
    %pony1 = kony1/(kony1+kfby1*(N-Ky1));
    ponx1 = kon/(kon+kfb*(N-Kx1));
    pony1 = kon/(kon+kfb*(N-Ky1));
    
    locx2 = posy2(1:Ky2,t)+(1+4*rand(1,1))*dx; % randomly chosen loc not Y2(t)
    locy2 = posx2(1:Kx2,t)+(1+4*rand(1,1))*dx; % randomly chosen loc not X2(t)
    %ponx2 = konx2/(konx2+kfbx2*(N-Kx2));
    %pony2 = kony2/(kony2+kfby2*(N-Ky2));
    ponx2 = kon/(kon+kfb*(N-Kx2));
    pony2 = kon/(kon+kfb*(N-Ky2));
    
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
    
end
toc

%% Plot all particle trajectories
%
ccx = [255 219 88]/256.*ones(Nt,1);     % mustard yellow
ccy = [0 0 255]/256.*ones(Nt,1);        % blue
figure(11);
subplot(1,2,1);
for j=1:5:max(max([NNx1,NNy1]))
    hold on;
    %plot(linspace(0,Tend,Nt),pos(j,:));
    scatter(linspace(0,Tend,Nt),posx1(j,:),1,ccx);
    scatter(linspace(0,Tend,Nt),posy1(j,:),1,ccy);
    pbaspect([10 4 1]); 
    set(gca,'color','k','fontname','times','fontsize',20); box on;
    set(gcf,'color','w');
    xlabel('Time');
    ylabel('Cell 1');
end

subplot(1,2,2);
for j=1:5:max(max([NNx2,NNy2]))
    hold on;
    %plot(linspace(0,Tend,Nt),pos(j,:));
    scatter(linspace(0,Tend,Nt),posx2(j,:),1,ccx);
    scatter(linspace(0,Tend,Nt),posy2(j,:),1,ccy);
    pbaspect([10 4 1]);
    set(gca,'color','k','fontname','times','fontsize',20); box on;
    set(gcf,'color','w');
    xlabel('Time');
    ylabel('Cell 2');
end


