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
%vidObj = VideoWriter('stronglycoupled.avi');

counter_ppp = 1;
ppp = 1;
       
while (ppp<=1)
    counter_ppp = counter_ppp+1;

    clearvars -except counter_ppp vid vidObj ppp
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
    alpha = 2.0;
    beta = 2.0;

    % Set discretization 
    %
    L      = 10.0;                  % cell length
    dt     = 0.01;                  % temporal discretization
    Na     = 101;                   % number of space steps
    dxa    = 5.0/((Na-1)/2);        % spatial discretization
    Xa     = 0:dxa:L;
    pa     = dt*Da/(dxa^2);  
    Tend   = 25.0;                  % total simulation time
    Nt     = Tend/dt;
    dx     = sqrt(2*D*dt);
    tplot  = 100; 

    posx = zeros(N,Nt);              % array of positions of X(t)
    posy = zeros(N,Nt);              % array of positions of Y(t)

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
        steepness = 20;
        a = (tanh(steepness*(X-0.375)) - tanh(steepness*(X-1.125)) + 0.2)/2.2;
        b = (2 - tanh(steepness*(X-0.375)) + tanh(steepness*(X-1.125)) + 0.2)/2.2;

        %a = (tanh(steepness*(X-0.5)) - tanh(steepness*(X-1.5)) + 0.2)/2.2;
        %b = (2 - tanh(steepness*(X-0.5)) + tanh(steepness*(X-1.5)) +0.2)/2.2;
    elseif (ictype==4)
    % (4) odd condition #1 (multiple peaks)
        steepness = 10;
        a = (1-cos(3*Xa*pi/5))/2; a=a';
        b = (tanh(steepness*(Xa-7.5))+1)/2 + (1-tanh(steepness*(Xa-2.5)))/2; b=b';

    elseif (ictype==5)
    % (4) odd condition #2
        steepness = 10;
        b = (1-cos(Xa*pi/5))/2; b=b';
        a = (tanh(steepness*(Xa-7.5))+1)/2 + (1-tanh(steepness*(Xa-2.5)))/2; a=a';

    elseif (ictype==6)
    % (5) odd condition #3
        mu = 1.8; sigma = 0.1;
        a = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); a=a';
        mu = 1.9; sigma = 0.1;
        b = awgn(exp(-0.5*((Xa-mu)/sigma).^2)./(sigma*sqrt(32*pi)),20); b=b';
    end

    % Laplacian difference operator with no flux boundary conditions
    % Crank-Nicolson operators
    II = speye(Na,Na);
    Lapdiff = spdiags(ones(Na,1)*[0 1 -2 1 0], [-Na+1, -1, 0 , 1, Na-1], Na, Na);
    Lapdiff(1,1) = -2; Lapdiff(1,2) = 1; Lapdiff(1,Na) = 1;
    Lapdiff(Na,1) = 1; Lapdiff(Na,Na-1) = 1; Lapdiff(Na,Na) = -2;
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
    nx(1:X(1),1)      = 1;                     % activate mem-bound particles
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
    figure(ppp); 
    plot(Xa,a,'-o','markerfacecolor',[159 219 229]/255,'linewidth',3); hold on;
    plot(Xa,b,'-ok','markerfacecolor','k','linewidth',3);
    plot(s,xC,'-.','color',[0 0.45 0.75],'linewidth',1);
    plot(s,yC,'-.k','linewidth',1);
    xlim([0 10]); ylim([0 2]);
    %title('Time = 0');
    set(gca,'fontname','times','fontsize',20); box on;
    lgd = legend('Branched network','Bundled network','Rac','Rho','Location','northeast');
    lgd.NumColumns = 2;       
    set(gcf,'color','w');
    hold off;
    %keyboard
    pause(1.0);

    if vid==1
        currFrame = getframe(gcf);
        writeVideo(vidObj,currFrame);
    end

    %% Run simulation
    % 
    tic
    quit_cond = 0;
    cond = 0;
    for t=1:(Nt-1)

        %% Run biochemisty
        [Konx,Kony,Kfbx,Kfby,Koffx,Koffy] = spatialrates(ron,rfb,roff,a,b,s,beta,cond); % set rates

        if((t-1)*dt<Tx(rxn_count_x))
            NNx(t+1) = X(rxn_count_x-1);
        else
            nnx = X(rxn_count_x);
            taux = zeros(nnx,1);
            dn = zeros(nnx,1);
            r = rand(nnx,1);

            if(nnx==0)
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nnx          % all agents
                konx = interp1(s,Konx,posx(j,t));
                koffx = interp1(s,Koffx,posx(j,t));
                kfbx = interp1(s,Kfbx,posx(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffx + (konx+kfbx*nnx/N)*(N/nnx-1);
                taux(j) = -log(r(j))/a0;
                rr = rand(1,1);
                dn(j) = (rr<((konx+kfbx*nnx/N)*(N/nnx-1)/a0))*1.0 + (rr>=((konx+kfbx*nnx/N)*(N/nnx-1)/a0))*(-1.0);
            end

            [mintaux,minidx] = min(taux(1:j));       % find first chemical rxn
            Tx(rxn_count_x+1) = Tx(rxn_count_x) + mintaux;
            X(rxn_count_x+1) = nnx + dn(minidx);
            rxn_count_x = rxn_count_x + 1;
            NNx(t+1) = X(rxn_count_x-1);
        end

        if((t-1)*dt<Ty(rxn_count_y))
            NNy(t+1) = Y(rxn_count_y-1);
        else
            nny = Y(rxn_count_y);
            tauy = zeros(nny,1);
            dn = zeros(nny,1);
            r = rand(nny,1);

            if(nny==0)
                counter_ppp = ppp;
                quit_cond = 1;
                break
            end

            for j=1:nny          % all agents
                kony = interp1(s,Kony,posy(j,t));
                koffy = interp1(s,Koffy,posy(j,t));
                kfby = interp1(s,Kfby,posy(j,t));
                % Sample earliest time-to-fire (tau)
                a0 = koffy + (kony+kfby*nny/N)*(N/nny-1);
                tauy(j) = -log(r(j))/a0;
                rr = rand(1,1);
                dn(j) = (rr<((kony+kfby*nny/N)*(N/nny-1)/a0))*1.0 + (rr>=((kony+kfby*nny/N)*(N/nny-1)/a0))*(-1.0);
            end

            [mintauy,minidy] = min(tauy(1:j));       % find first chemical rxn      
            Ty(rxn_count_y+1) = Ty(rxn_count_y) + mintauy;
            Y(rxn_count_y+1) = nny + dn(minidy);
            rxn_count_y = rxn_count_y + 1;
            NNy(t+1) = Y(rxn_count_y-1);
        end

        if (quit_cond==1)
            sprintf("It's happening at kk = %d, ppp = %d\n",kk,ppp)
            ppp = counter_ppp-1;
            break
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
        % Resolution strategy: No one advances
        %
        firstcoll = sum(ismembertol(posx(1:K1,t+1),posy(1:K2,t+1),0.005,'DataScale',1));
        if firstcoll~=0           
            % Get indices of collisions
            aa = ismembertol(posx(1:K1,t+1),posy(1:K2,t+1),0.005,'DataScale',1);
            list_idx = find(aa~=0);
            bb = ismembertol(posy(1:K2,t+1),posx(1:K1,t+1),0.005,'DataScale',1);
            list_idy = find(bb~=0);

            posx(list_idx,t+1) = posx(list_idx,t);
            posy(list_idy,t+1) = posy(list_idy,t);
        end

        % Enforce periodic boundary conditions
        posx(1:K1,t+1) = posx(1:K1,t+1) + (-L).*(posx(1:K1,t+1)>L) + (L).*(posx(1:K1,t+1)<0.0);
        posy(1:K2,t+1) = posy(1:K2,t+1) + (-L).*(posy(1:K2,t+1)>L) + (L).*(posy(1:K2,t+1)<0.0);

        % Enforce no-flux boundary conditions
        %posx(1:K1,t+1) = posx(1:K1,t+1) + (posx(1:K1,t)-posx(1:K1,t+1)).*(posx(1:K1,t+1)>L) + (posx(1:K1,t)-posx(1:K1,t+1)).*(posx(1:K1,t+1)<0.0);
        %posy(1:K2,t+1) = posy(1:K2,t+1) + (posy(1:K2,t)-posy(1:K2,t+1)).*(posy(1:K2,t+1)>L) + (posy(1:K2,t)-posy(1:K2,t+1)).*(posy(1:K2,t+1)<0.0);

        %% Determine if a biochemical rxn has occured - update positions
        
        % Find spontaneous association location
        ss = sort(posx(1:K1,t));
        [ijk] = find(ss==posx(minidx,t),1);
        prevind = (ijk-1)*(ijk>1) + (K1)*(ijk==1);
        nextind = (ijk+1)*(ijk<K1) + 1*(ijk==K1);
        x2 = posx(minidx,t)+(ss(prevind)-posx(minidx,t))/2;
        x1 = posx(minidx,t)+(ss(nextind)-posx(minidx,t))/2;
        locx = (x2-x1).*rand(1,1) + x1; % random location halfway between the closest left/right particles
        ss = sort(posy(1:K2,t));
        [ijk] = find(ss==posy(minidy,t),1);
        prevind = (ijk-1)*(ijk>1) + (K2)*(ijk==1);
        nextind = (ijk+1)*(ijk<K2) + 1*(ijk==K2);
        y2 = posy(minidy,t)+(ss(prevind)-posy(minidy,t))/2;
        y1 = posy(minidy,t)+(ss(nextind)-posy(minidy,t))/2;
        locy = (y2-y1).*rand(1,1) + y1; % random location halfway between the closest left/right particles

        ponx = ron/(ron+rfb*(N-K1));
        pony = ron/(ron+rfb*(N-K2));

        if(NNx(t+1) < NNx(t))                % diassociation event (particle off)
            oldcol = posx(minidx,1:end);
            othercols = posx([1:minidx-1,minidx+1:K1],1:end);
            otherothercols = posx(K1+1:end,1:end);
            newpos = [othercols;oldcol;otherothercols];
            posx = newpos;
            nx(K1,t+1) = 0;
        elseif(NNx(t+1) > NNx(t))             % association event (on or recruitment)
            rr = rand(1,1); 
            posx(K1,t+1) = posx(K1,t)+(rr<ponx)*locx;              % on event
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
            posy(K2,t+1) = posy(K2,t)+(rr<pony)*locy;               % on event
            posy(K2,t+1) = posy(K2,t)+(rr>=pony)*posy(minidy,t);    % recruitment event
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

        %% Plot the solution(s)
        %if mod(t,tplot) == 0
        if t==(Nt-1)
            figure(ppp); 
            %titfig = sprintf('D = %.2f, m0 = %.2f, beta = %.2f, alpha = %.2f, k_{on} = %.4f, k_{off} = %.2f, k_{fb} = %.4f',D,m0,beta,alpha,ron,roff,rfb); 
            %sgtitle(titfig);
            %subplot(5,4,ppp);
            %figure(ppp); 
            plot(Xa,a,'-o','markerfacecolor',[159 219 229]/255,'linewidth',3); hold on;
            plot(Xa,b,'-ok','markerfacecolor','k','linewidth',3);
            plot(s,xC,'-.','color',[0 0.45 0.75],'linewidth',1);
            plot(s,yC,'-.k','linewidth',1);
            %title(['Time = ',num2str(t*dt,'%.2f')]);
            set(gca,'fontname','times','fontsize',20); box on;
            lgd = legend('Branched network','Bundled network','Rac','Rho','Location','northeast');
            lgd.NumColumns = 2;
            xlim([0 10]);
            set(gcf,'color','w');
            drawnow;
            %keyboard;
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
    sprintf('Simulation %d done',ppp)
    toc
    if(quit_cond==0)
        ppp = ppp + 1;
    end
end
  
%% Plot all particle trajectories
ccx = [0 0 255]/256.*ones(Nt,1);     % blue
ccy = [255 219 88]/256.*ones(Nt,1);  % mustard yellow
time = linspace(0,Tend,Nt);
figure(10);
for j=1:2:max(max([NNx,NNy]))
    hold on;
    %plot(linspace(0,Tend,Nt),pos(j,:));
    scatter(linspace(0,Tend,Nt),posx(j,:),1,ccx);
    scatter(linspace(0,Tend,Nt),posy(j,:),1,ccy);
    box on;
    set(gca,'Color','k','fontsize',20,'fontname','times');
    pbaspect([3 1 1]);
    set(gcf,'color','w');
    %xlabel('Time');
    %ylabel('Location on cell membrane');
    %yticks([0 0.5 1 1.5 2]);
    %yticklabels({'0','0.25','0.5','0.75','1'});
end