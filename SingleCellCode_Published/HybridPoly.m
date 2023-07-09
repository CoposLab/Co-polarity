% Simulate competition for resources for actin filaments
% and fully coupled to stochastic biochemistry for polarity proteins
% on a circular cell membrane (periodic BC, 1.5D)
%
% Mechanics -- two actin networks: branched (a) and contractile (b)
% Polarity proteins : Rac (X) and Rho (Y)
%
% Last updated: 1/31/2020
% Calina Copos
addpath('../TwoCellCode/freeze_colors')

clear;
close all;
clc;

savefigs = 0;
setnum='100';
savelocation='singlecellresults/';
if savefigs==1
    % filenameC1=strcat('savedgraphs/doubleRhoOnCell1_',setnum);
    % filenameC2=strcat('savedgraphs/doubleRhoOnCell2_',setnum);
    filenameCells=strcat(savelocation,'Cell_',setnum);
    filenameScatter=strcat(savelocation,'Scatter_',setnum);
end

vid = 0;
% vidObj = VideoWriter('stronglycoupled.mp4','MPEG-4');
% vidObjCol = VideoWriter('colorplot.mp4','MPEG-4');
% vidObjRR = VideoWriter('colorplotrr.mp4','MPEG-4');

counter_ppp = 1;
ppp = 1;

while (ppp<=1)
    counter_ppp = counter_ppp+1;

    clearvars -except counter_ppp vid vidObj ppp vidObjCol vidObjRR savefigs filenameCells filenameScatter
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
    beta = [2,2];

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

        vidObjCol.FrameRate = 5;
        vidObjCol.Quality = 75;
        open(vidObjCol);

        vidObjRR.FrameRate = 5;
        vidObjRR.Quality = 75;
        open(vidObjRR);
    end

    if vid ==1
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
        Z1 = [a a a a a]';
        Z2 = [b b b b b]';

        % figure(11); %Branched
        % % subplot(1,2,1);
        % contourf(Xcol,Ycol,Z1,100,'LineStyle','none')
        % axis square
        % col=colorbar;
        % colormap(whiteblue);
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
        % colormap(whitered);
        % clim([0,max(max(a),max(b))])
        % title('Bundled Actin')
        % ylabel(col,'Concentration')
        % set(gca,'XTick',[], 'YTick', [])

        %Plot both circles together
        colFrame = figure(13);
        % subplot(1,3,3);
        contourf(Xcol,Ycol,Z1-Z2,100,'LineStyle','none')
        axis square
        col=colorbar;
        colormap(redwhiteblue);
        clim([-max(max(a),max(b)),max(max(a),max(b))])
        title('Combined: Blue=Branched, Red=Bundled')
        ylabel(col,'Concentration Branched - Concentration Bundled')
        set(gca,'XTick',[], 'YTick', [])

        allred = [linspace(red(1),red(1),colorLength)',linspace(red(2),red(2),colorLength)',linspace(red(3),red(3),colorLength)'];
        allblue = [linspace(blue(1),blue(1),colorLength)',linspace(blue(2),blue(2),colorLength)',linspace(blue(3),blue(3),colorLength)'];

        figure(15)
        clf
        h1=surf(Xcol,Ycol,Z1,'FaceAlpha',0.4);
        view(2)
        colormap(whiteblue) % colormap(allblue)
        freezeColors
        % freezeColors(colorbar)
        clim([-max(max(a),max(b)),max(max(a),max(b))])
        shading interp
        % h1.AlphaData=gradient(Z1)./2;
        % h1.FaceAlpha = 'interp';
        hold on;
        h2=surf(Xcol,Ycol,Z2,'FaceAlpha',0.4);
        colormap(whitered) % colormap(allred)
        freezeColors
        % freezeColors(colorbar)
        shading interp
        % h2.AlphaData=gradient(Z2)./2;
        % h2.FaceAlpha = 'interp';
        title('Combined: Blue=Branched, Red=Bundled')
        hold off;

        % [Ma,Ia] = max(a);
        % [Mb,Ib] = max(b);
        %
        % hold on;
        % plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
        % plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
        % plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
        % hold off;

        % Concentric circles
        [th,rad] = meshgrid((0:3.6:360)*pi/180,0.76:0.01:0.8);
        [Xsm,Ysm] = pol2cart(th,rad);

        figure(16)
        clf
        ax1=axes;
        ax2=axes;
        s1=surf(ax2,Xcol,Ycol,Z1);
        view(2)
        colormap(whiteblue);
        % cb1=colorbar(ax1,'Position',[.88 .11 .0675 .815]);
        freezeColors
        freezeColors(colorbar)
        hold on;
        s2=surf(ax2,Xsm,Ysm,Z2);
        colormap(whitered);
        freezeColors
        cb2=colorbar(ax1,'Position',[.05 .11 .0675 .815]);
        clim([0,max(max(a),max(b))])
        linkaxes([ax1,ax2])
        ax1.Visible = 'off';
        ax1.XTick = [];
        ax1.YTick = [];
        ax2.XTick = [];
        ax2.YTick = [];
        shading interp
        axis square
        grid off
        title('Combined: Blue=Branched, Red=Bundled')
        hold off;

        if vid==1
            currFrame = getframe(colFrame);
            writeVideo(vidObjCol,currFrame);
        end

    end

    % Plot Rac and Rho
    % colRRFrame = figure(14);
    % ZRac = [xC xC xC xC xC]';
    % ZRho = [yC yC yC yC yC]';
    % contourf(Xcol,Ycol,ZRac-ZRho,100,'LineStyle','none');
    % axis square
    % col=colorbar;
    % colormap(redwhiteblue);
    % clim([-max(max(xC),max(yC)),max(max(xC),max(yC))])
    % title('Combined: Blue=Rac, Red=Rho')
    % ylabel(col,'Concentration Rac - Concentration Rho')
    % set(gca,'XTick',[], 'YTick', [])

    % if vid==1
    %     currFrame = getframe(colRRFrame);
    %     writeVideo(vidObjRR,currFrame);
    % end

    %% Run simulation
    %
    tic
    quit_cond = 0;
    cond = 0;
    for t=1:(Nt-1)

        %% Run biochemisty
        [Konx,Kony,Kfbx,Kfby,Koffx,Koffy] = spatialrates(ron,rfb,roff,a,b,s,beta,cond,[]); % set rates

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
        % if mod(t,tplot) == 0
        if t==(Nt-1)
            scatplot=figure(ppp);
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

            %Plot on circle
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
            ZBranch1 = [a a a a a a a a]';
            ZBund1 = [b b b b b b b b]';
            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
            [Xsm,Ysm] = pol2cart(th,rad);
            [th,rad] = meshgrid((0:3.6:360)*pi/180,0.83:0.01:0.9);
            [Xmid,Ymid] = pol2cart(th,rad);


            % Concentric circles
            % Cell 1
            figcells=figure(16);
            surf(Xcol,Ycol,ZBranch1);
            view(2)
            colormap(whitebluenavy)
            freezeColors;
            freezeColors(colorbar);
            clim([0,max(b)])
            shading interp
            hold on;
            surf(Xmid,Ymid,ZBund1);
            colormap(whitedarkyellow)
            freezeColors;
            freezeColors(jicolorbar);
            clim([0,max(a)])
            shading interp
            grid off
            set(gca,'XTick',[], 'YTick', [])
            % scatter(Xsm(boundC1),Ysm(boundC1),'black');
            hold off;
            box off;
            set(gca,'XColor','w')
            set(gca,'YColor','w')
            set(gcf,'color','w');
            axis square



            % Find median for cell 1
            a1New = a;
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


        end
    end

    % measure of polarized state (1 if polarized and 0 otherwise)
    %st = 1*( (abs(a(1)-b(end))<1e-3 || abs(a(end)-b(1))<1e-3 ) && abs(a(1)-a(end))>1e-3 && abs(b(1)-b(end))>1e-3 );

    if vid==1
        close(vidObj);
        close(vidObjCol);
        close(vidObjRR);
    end
    sprintf('Simulation %d done',ppp)
    toc
    if(quit_cond==0)
        ppp = ppp + 1;
    end
end

if savefigs==1
    savefig(figcells,filenameCells);
    savefig(scatplot,filenameScatter);
end

%% Plot all particle trajectories
% ccx = [0 0 255]/256.*ones(Nt,1);     % blue
% ccy = [255 219 88]/256.*ones(Nt,1);  % mustard yellow
% time = linspace(0,Tend,Nt);
% figure(10);
% for j=1:2:max(max([NNx,NNy]))
%     hold on;
%     %plot(linspace(0,Tend,Nt),pos(j,:));
%     scatter(linspace(0,Tend,Nt),posx(j,:),1,ccx);
%     scatter(linspace(0,Tend,Nt),posy(j,:),1,ccy);
%     box on;
%     set(gca,'Color','k','fontsize',20,'fontname','times');
%     pbaspect([3 1 1]);
%     set(gcf,'color','w');
%     %xlabel('Time');
%     %ylabel('Location on cell membrane');
%     %yticks([0 0.5 1 1.5 2]);
%     %yticklabels({'0','0.25','0.5','0.75','1'});
% end