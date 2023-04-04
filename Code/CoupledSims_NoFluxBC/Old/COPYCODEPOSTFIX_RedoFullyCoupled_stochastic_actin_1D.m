% Simulate competition for resources for actin filaments 
% and fully coupled to stochastic biochemistry for polarity proteins
%
% Mechanics -- two actin networks: branched (a) and contractile (b)
% Polarity proteins : Rac (X) and Rho (Y)
%
% Last updated: 5/14/2019

clear;
close all;

vid = 1;
vidObj = VideoWriter('sim_fullycoupled05_odd3ic.avi');

% Set actin filament parameters
% 
Da      = 0.1;                  % diffusion coefficient for actin
m0      = 2;                    % competition for actin monomers
%kon     = 0.0001;
%koff    = 0.11;
%kfb     = 0.1;
% Whaaaaat are these parameters for?!?!

% Set polarity protein parameters
% 
N       = 100;                  % total number of molecules in the cell (conserved)
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
pa     = dt*Da/(dxa^2);  
Tend   = 30.0;
Nt     = Tend/dt;
dx     = sqrt(2*D*dt);
tplot  = 500; 

posx = zeros(N,Nt);              % array of positions of X(t)
posy = zeros(N,Nt);              % array of positions of Y(t)

% Actin reaction term
%
%F = @(U,V) U-U.*U - m0*U.*V; 
F = @(U,V) -U.*U - m0*U.*V; 

% Set initial conditions for actin distribution
%
ictype = 6;
%    1 = step in the middle
%    2 = random
%    3 = sigmoidal curve
%    4 = odd condition #1: branched peak in the middle, contractile peaks at
%    the front and rear cell
%    5 = odd condition #2: both peaks in the middle (w/ noise)
%    6 = odd condition #1: bundled peak in the middle, branched peaks at
%    the front and rear cell

a       = zeros(N,1);
anew    = zeros(N,1);
b       = zeros(N,1);
bnew    = zeros(N,1);

% (1) pulse in middle
if (ictype==1)
    a   = ones(N,1);
    anew= ones(N,1);
  
    b   = 0.5*ones(N,1);        
    bnew= 0.5*ones(N,1);

% (2) random
elseif (ictype==2)
    a = 0.1 + 0.9.*rand(length(Xa),1);
    b = 0.1 + 0.9.*rand(length(Xa),1);

% (3) arctangent 
elseif (ictype==3)
    steepness = 25;
    a = (tanh(steepness*(Xa-1))+1)/2;
    b = (1-tanh(steepness*(Xa-1)))/2;
   
elseif (ictype==4)
% (4) odd condition #1
    steepness = 25;
    a = (1-cos(Xa*pi))/2;
    b = (tanh(steepness*(Xa-1.5))+1)/2 + (1-tanh(steepness*(Xa-0.5)))/2;
    
elseif (ictype==5)
% (5) odd condition #2
    mu = 0.8; sigma = 0.1;
    a = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20);
    mu = 0.8; sigma = 0.1;
    b = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20);
    
elseif (ictype==6)
% (5) odd condition #3
    steepness = 25;
    b = (1-cos(Xa*pi))/2;
    a = (tanh(steepness*(Xa-1.5))+1)/2 + (1-tanh(steepness*(Xa-0.5)))/2;
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
nx   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
ny   = zeros(N,Nt);                      % state of the particle (0 inactive, 1 active)
Tx   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for X(t)
Ty   = zeros(MAX_OUTPUT_LENGTH,1);       % times of chemical reactions for Y(t)
X    = zeros(MAX_OUTPUT_LENGTH,1);       % number of X(t) molecules on the membrane    
Y    = zeros(MAX_OUTPUT_LENGTH,1);       % number of Y(t) molecules on the membrane
NNx  = zeros(Nt,1);
NNy  = zeros(Nt,1);

% Set initial conditions for polarity molecules distribution
%
rxn_count_x       = 1;
rxn_count_y       = 1;
X(1)              = 0.1*N;                 % # of particles on membrane
Y(1)              = 0.1*N;                 % # of particles on membrane
NNx(1)            = X(1);
NNy(1)            = Y(1);
Tx(1)             = 0.0;                   
Ty(1)             = 0.0;
nx(1:X(1),1)      = 1;                    % activate mem-bound particles
ny(1:X(1),1)      = 1;
r = randperm(ceil(L/(0.0102)),X(1)+Y(1))*0.0102;
posx(1:X(1),1)=r(1:X(1));                         
posy(1:Y(1),1)=r(X(1)+1:end);

% Sample concentration at actin filament spatial scale
%
[s,xC,yC] = resamplePolarityMolecules(posx(1:NNx(1),1),posy(1:NNy(1),1),NNx(1),NNy(1),L,Na);

aic = a;
bic = b;

% Setup convergence checks for actin quantities
%
conv1 = zeros(Nt,2);
conv2 = zeros(Nt,2);
convsup = zeros(Nt,2);

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
    figure(1); 
    plot(Xa,a,'-bo','markerfacecolor',[65 105 225]/255,'linewidth',2); hold on;
    %plot(Xa,a,'-o','markerfacecolor',[159 219 229]/255,'linewidth',2); hold on;
    plot(Xa,b,'-dk','markerfacecolor','k','linewidth',2);
    %plot(s,xC,'--','color',[0 0.45 0.75],'linewidth',2);
    plot(s,xC,'--b','linewidth',2);
    plot(s,yC,'--k','linewidth',2);
    title('Time = 0.0');
    xlim([0 2]); %ylim([0 35]);
    xticks([0 0.5 1 1.5 2]);
    %xticklabels({'0','0.25','0.5','0.75','1'});
    %xlabel('Location on cell membrane');
    %ylabel('Concentration');
    %set(gcf,'color','w');
    set(gca,'fontname','times','fontsize',30); box on;
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
    %if t>=10000
    %    cond=1;
    %end
    
    %% Run biochemisty
    [Konx,Kony,Kfbx,Kfby,Koffx,Koffy] = spatialrates(ron,rfb,roff,a,b,s,beta,cond); % set rates
    
    if((t-1)*dt<Tx(rxn_count_x))
        NNx(t+1) = X(rxn_count_x-1);
    else
       	nx = X(rxn_count_x);
        taux = zeros(nx,1);
        dn = zeros(nx,1);
        r = rand(nx,1);
        
        for j=1:nx          % all agents
            konx = interp1(s,Konx,posx(j,t));
            koffx = interp1(s,Koffx,posx(j,t));
            kfbx = interp1(s,Kfbx,posx(j,t));
            % Sample earliest time-to-fire (tau)
            %a0 = koff + (konx+kfbx*nx/N)*(N/nx-1);
            % fixed 7/15
            a0 = koffx + (konx+kfbx*nx/N)*(N/nx-1);
            taux(j) = -log(r(j))/a0;
            rr = rand(1,1);
            dn(j) = (rr<((konx+kfbx*nx/N)*(N/nx-1)/a0))*1.0 + (rr>=((konx+kfbx*nx/N)*(N/nx-1)/a0))*(-1.0);
        end
        
        [mintaux,minidx] = min(taux(1:j));       % find first chemical rxn
        Tx(rxn_count_x+1) = Tx(rxn_count_x) + mintaux;
        X(rxn_count_x+1) = nx + dn(minidx);
        rxn_count_x = rxn_count_x + 1;
        NNx(t+1) = X(rxn_count_x-1);
    end
    
    if((t-1)*dt<Ty(rxn_count_y))
        NNy(t+1) = Y(rxn_count_y-1);
    else
        ny = Y(rxn_count_y);
        tauy = zeros(ny,1);
        dn = zeros(ny,1);
        r = rand(ny,1);
        
        for j=1:ny          % all agents
            kony = interp1(s,Kony,posy(j,t));
            koffy = interp1(s,Koffy,posy(j,t));
            kfby = interp1(s,Kfby,posy(j,t));
            % Sample earliest time-to-fire (tau)
            %a0 = koff + (kony+kfby*ny/N)*(N/ny-1);
            % fixed 7/15
            a0 = koffy + (kony+kfby*ny/N)*(N/ny-1);
            tauy(j) = -log(r(j))/a0;
            rr = rand(1,1);
            dn(j) = (rr<((kony+kfby*ny/N)*(N/ny-1)/a0))*1.0 + (rr>=((kony+kfby*ny/N)*(N/ny-1)/a0))*(-1.0);
        end
        
        [mintauy,minidy] = min(tauy(1:j));       % find first chemical rxn      
        Ty(rxn_count_y+1) = Ty(rxn_count_y) + mintauy;
        Y(rxn_count_y+1) = ny + dn(minidy);
        rxn_count_y = rxn_count_y + 1;
        NNy(t+1) = Y(rxn_count_y-1);
    end

    %% Run diffusion of membrane-bound polarity proteins
    p  = 0.5;                  % probability of hoping left or right

    % Fetch the number of particles at this time
    K1 = NNx(t+1);
    K2 = NNy(t+1);

    % Between reactions, perform Brownian motion with periodic BC
    r = rand(K1,1);    % coin flip
    nx(1:K1,t+1) = 1;   
    posx(1:K1,t+1) = posx(1:K1,t) + dx*((r<p)*1.0 + (r>(1-p))*(-1.0));

    r = rand(K2,1);    % coin flip
    ny(1:K2,t+1) = 1;
    posy(1:K2,t+1) = posy(1:K2,t) + dx*((r<p)*1.0 + (r>(1-p))*(-1.0));

    % Check for collision(s) and resolve any collisions
    %
    % Resolution: winner moves, loser stays goes back to previous timestep
    % (did not work)
    % New resolution: no one advances
    % (works)
    firstcoll = sum(ismembertol(posx(1:K1,t+1),posy(1:K2,t+1),0.005,'DataScale',1));  
    if firstcoll~=0           
        %sprintf('collision occured @ t=%.2f\n',t+1)
        % Get indices of collisions
        aa = ismembertol(posx(1:K1,t+1),posy(1:K2,t+1),0.005,'DataScale',1);
        list_idx = find(aa~=0);
        bb = ismembertol(posy(1:K2,t+1),posx(1:K1,t+1),0.005,'DataScale',1);
        list_idy = find(bb~=0);

        posx(list_idx,t+1) = posx(list_idx,t);
        posy(list_idy,t+1) = posy(list_idy,t);

    %         i = firstcoll;
    %         while (i>0) 
    %             r = rand(1,1); 
    %             posx(list_idx(i),t+1) = (r<0.5)*posx(list_idx(i),t+1) + (r>0.5)*(posx(list_idx(i),t)-dx);
    %             posy(list_idy(i),t+1) = (r<0.5)*(posy(list_idy(i),t)-dx) + (r>0.5)*posy(list_idy(i),t+1);  
    %             i = i-1;
    %         end
    end

    % Enforce periodic boundary conditions
    %posx(1:K1,t+1) = posx(1:K1,t+1) + (-L).*(posx(1:K1,t+1)>L) + (L).*(posx(1:K1,t+1)<0.0);
    %posy(1:K2,t+1) = posy(1:K2,t+1) + (-L).*(posy(1:K2,t+1)>L) + (L).*(posy(1:K2,t+1)<0.0);

    % Enforce no-flux boundary conditions
    posx(1:K1,t+1) = posx(1:K1,t+1) + (posx(1:K1,t)-posx(1:K1,t+1)).*(posx(1:K1,t+1)>L) + (posx(1:K1,t)-posx(1:K1,t+1)).*(posx(1:K1,t+1)<0.0);
    posy(1:K2,t+1) = posy(1:K2,t+1) + (posy(1:K2,t)-posy(1:K2,t+1)).*(posy(1:K2,t+1)>L) + (posy(1:K2,t)-posy(1:K2,t+1)).*(posy(1:K2,t+1)<0.0);
    
    % Determine if a biochemical rxn has occured - update positions
    locx = posy(1:K2,t)+(1+4*rand(1,1))*dx; % randomly chosen loc not Y(t)
    locy = posx(1:K1,t)+(1+4*rand(1,1))*dx; % randomly chosen loc not X(t)
    %ponx = konx/(konx+kfbx*(N-K1));
    %pony = kony/(kony+kfby*(N-K2));
    ponx = ron/(ron+rfb*(N-K1));
    pony = ron/(ron+rfb*(N-K2));
    % is this a typo? should it be konx, kony instead of kon, kon 
    % and similarly for kfb. found on 7/14. fixed. not checked if all is
    % good... maybe not but definitely should be ron, rfb. (7/15)
    
    
    if(NNx(t+1) < NNx(t))                % diassociation event (particle off)
        oldcol = posx(minidx,1:end);
        othercols = posx([1:minidx-1,minidx+1:K1],1:end);
        otherothercols = posx(K1+1:end,1:end);
        newpos = [othercols;oldcol;otherothercols];
        posx = newpos;
        nx(K1,t+1) = 0;
    elseif(NNx(t+1) > NNx(t))             % association event (on or recruitment)
        rr = rand(1,1); 
        posx(K1,t+1) = posx(K1,t)+(rr<ponx)*locx(randi(K2,1));   % on event
        posx(K1,t+1) = posx(K1,t)+(rr>=ponx)*posx(minidx,t);   % recruitment event
        nx(K1+1,t+1) = 1;
    end

    if (NNy(t+1) < NNy(t))                % diassociation event (particle off)
        oldcol = posy(minidy,1:end);
        othercols = posy([1:minidy-1,minidy+1:K2],1:end);
        otherothercols = posy(K2+1:end,1:end);
        newpos = [othercols;oldcol;otherothercols];
        posy = newpos;
        ny(K2,t+1) = 0;
    elseif(NNy(t+1) > NNy(t))             % association event (on or recruitment)
        rr = rand(1,1); 
        posy(K2,t+1) = posy(K2,t)+(rr<pony)*locy(randi(K1,1));   % on event
        posy(K2,t+1) = posy(K2,t)+(rr>=pony)*posy(minidy,t);   % recruitment event
        ny(K2,t+1) = 1;
    end
    
    [s,xC,yC] = resamplePolarityMolecules(posx(1:K1,t+1),posy(1:K2,t+1),K1,K2,L,Na);
  
    %% Update actin filaments
    diffRHSa = Hm*a;
    diffRHSb = Hm*b;
  
    rxna = dt*( F(a,b) + a.*(1+alpha*xC));
    rxnb = dt*( F(b,a) + b.*(1+alpha*yC));

    a = Hs\(diffRHSa+rxna);
    b = Hs\(diffRHSb+rxnb);

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
    if mod(t,tplot) == 0
    %if t==(Nt-1)
        plot(Xa,a,'-bo','markerfacecolor',[65 105 225]/255,'linewidth',2); hold on;
        plot(Xa,b,'-dk','markerfacecolor','k','linewidth',2);
        plot(s,xC,'--b','linewidth',2);
        plot(s,yC,'--k','linewidth',2);
    
        %title(['Time = ',num2str(t*dt)]);
        title(['Time = ',num2str(t)]);
        xlim([0 2]); %ylim([0 1.0]);
        set(gca,'fontname','times','fontsize',30); box on;
        legend('Branched network','Contractile network');%,'Rac','Rho','Location','northeast');
        %legend('Branched network (-- Rac)','Bundled network (-- Rho)','Location','northeast');
        xticks([0 0.5 1 1.5 2]);
        %xticklabels({'0','0.25','0.5','0.75','1'});
        %xlabel('Location on cell membrane');
        %ylabel('Concentration');
        %set(gcf,'color','w');
        drawnow;
        hold off;
        
        if vid==1
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame);
        end
    end 
    
    %% Check convergence
    conv2(t,1) = norm(a-aic,2);
    conv2(t,2) = norm(b-bic,2);
    conv1(t,1) = norm(a-aic,1);
    conv1(t,2) = norm(b-bic,1);
    convsup(t,1) = max(abs(a-aic));
    convsup(t,2) = max(abs(b-bic));
end

% measure of polarized state (1 if polarized and 0 otherwise)
%st = 1*( (abs(a(1)-b(end))<1e-3 || abs(a(end)-b(1))<1e-3 ) && abs(a(1)-a(end))>1e-3 && abs(b(1)-b(end))>1e-3 );     

if vid==1    
    close(vidObj);
end
toc

% Plot all particle trajectories
% 
% ccx = [255 219 88]/256.*ones(Nt,1);     % mustard yellow
% ccy = [0 0 255]/256.*ones(Nt,1);        % blue
% time = linspace(0,Tend,Nt);
% figure(2);
% for j=1:max(max([NNx,NNy]))
%     hold on;
%     %plot(linspace(0,Tend,Nt),pos(j,:));
%     scatter(linspace(0,Tend,Nt),posx(j,:),1,ccx);
%     scatter(linspace(0,Tend,Nt),posy(j,:),1,ccy);
%     box on;
%     set(gca,'Color','k','fontsize',30); mjyk=
%     xlabel('Time');
%     ylabel('Location on cell membrane');
%     yticks([0 0.5 1 1.5 2]);
%     yticklabels({'0','0.25','0.5','0.75','1'});
% end

% Check convergence
% 
% figure(3)
% plot(linspace(0,Tend,Nt),conv1(:,1))
% hold on;
% plot(linspace(0,Tend,Nt),conv2(:,1))
% plot(linspace(0,Tend,Nt),convsup(:,1))