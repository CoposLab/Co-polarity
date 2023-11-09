set(0,'DefaultFigureVisible','off')
close all;
clear;
clc;

addpath('./freeze_colors')

signal=0;
% Na=101;
% sigper=0.40;
% sigBound1 = (floor((Na-1)*1/8 - floor((Na-1)*sigper/2)))+1:(floor((Na-1)*1/8 + floor((Na-1)*sigper/2)))+1;
% sigBound1(sigBound1<=0)=sigBound1(sigBound1<=0)+Na;
% sigBound1(sigBound1>Na)=sigBound1(sigBound1>Na)-Na;
% sigBound2 = (floor((Na-1)*5/8 - floor((Na-1)*sigper/2)))+1:(floor((Na-1)*5/8 + floor((Na-1)*sigper/2)))+1;
% sigBound2(sigBound2<=0)=sigBound2(sigBound2<=0)+Na;
% sigBound2(sigBound2>Na)=sigBound2(sigBound2>Na)-Na;

scatvid=0;
branchedbundledvid=0;
racrhovid=0;
circlescatvid=1;
adjacent=1;

for i=1:3

    loadfile='./vid_matfiles/uncoupled_to_coupled/1000stepsuncoupled_1500stepscoupled/racup/1000RacOn';
    
    setnum=int2str(i);
    savelocation='../movies/uncoupled_to_coupled/1000stepsuncoupled_1500stepscoupled/racup/1000RacOn';

    if scatvid==1
        vidObj1 = VideoWriter(strcat(savelocation,'ScatterVid_',setnum,'.mp4'),'MPEG-4');
    end
    if branchedbundledvid==1
        vidObj2 = VideoWriter(strcat(savelocation,'BranchedBundledVid_',setnum,'.mp4'),'MPEG-4');
    end
    if racrhovid==1
        vidObj3 = VideoWriter(strcat(savelocation,'RacRhoVid_',setnum,'.mp4'),'MPEG-4');
    end
    if circlescatvid==1
        vidObj4 = VideoWriter(strcat(savelocation,'CircleScatterVid_',setnum,'.mp4'),'MPEG-4');
    end

    if scatvid==1
        vidObj1.FrameRate = 5;
        vidObj1.Quality = 100;
    end
    if branchedbundledvid==1
        vidObj2.FrameRate = 5;
        vidObj2.Quality = 100;
    end
    if racrhovid==1
        vidObj3.FrameRate = 5;
        vidObj3.Quality = 100;
    end
    if circlescatvid==1
        vidObj4.FrameRate = 5;
        vidObj4.Quality = 100;
    end

    if scatvid==1
        open(vidObj1);
    end
    if branchedbundledvid==1
        open(vidObj2);
    end
    if racrhovid==1
        open(vidObj3);
    end
    if circlescatvid==1
        open(vidObj4);
    end

    % allmax=0;
    % for t=1:49
    %     load(strcat(loadfile,int2str(t*50)));
    %     allmax = max(max(max(max(a1),max(a2)),max(max(b1),max(b2))),allmax);
    % end

    load(strcat(loadfile,setnum));
    set(0,'DefaultFigureVisible','off')

    allmax = max(max(max(max(a1all)),max(max(a2all))),max(max(max(b1all)),max(max(b2all))));

    for t=1:50:2499
    % for t=1:5:150
        a1=a1all(:,t);
        a2=a2all(:,t);
        b1=b1all(:,t);
        b2=b2all(:,t);
        xC1=xC1all(:,t);
        yC1=yC1all(:,t);
        xC2=xC2all(:,t);
        yC2=yC2all(:,t);

        L      = 10.0;
        Na=101;

        % [s1,xC1,yC1] = resamplePolarityMolecules(posx1(1:NNx1(t+1),t+1),posy1(1:NNy1(t+1),t+1),NNx1(t+1),NNy1(t+1),L,Na);
        % [s2,xC2,yC2] = resamplePolarityMolecules(posx2(1:NNx2(t+1),t+1),posy2(1:NNy2(t+1),t+1),NNx2(t+1),NNy2(t+1),L,Na);

        % allmax = max(max(max(a1),max(a2)),max(max(b1),max(b2)));

        % Define circles
        gapsize=0.01;
        [th,rad] = meshgrid((0:3.6:360)*pi/180,0.85:0.01:1);
        [Xcol,Ycol] = pol2cart(th,rad);
        Ycol1=Ycol;
        Ycol2=Ycol;
        Ycol1(:,boundC1)=Ycol1(:,boundC1(1)*ones(1,length(boundC1)));
        Ycol2(:,boundC2)=Ycol2(:,boundC2(1)*ones(1,length(boundC2)));
        Ycol2 = Ycol2 - 2*abs(max(max(Ycol2)))-gapsize;
        ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1]';
        ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1]';
        ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2]';
        ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2]';
        [th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
        [Xsm,Ysm] = pol2cart(th,rad);
        Ysm1=Ysm;
        Ysm2=Ysm;
        Ysm1(:,boundC1)=Ysm1(:,boundC1(1)*ones(1,length(boundC1)));
        Ysm2(:,boundC2)=Ysm2(:,boundC2(1)*ones(1,length(boundC2)));
        [th,rad] = meshgrid((0:3.6:360)*pi/180,0.86:0.01:0.93);
        [Xmid,Ymid] = pol2cart(th,rad);


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

        if scatvid==1
            % Make scatterplots
            scatplot=figure(1);
            clf
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

            scatframe = getframe(scatplot);
            writeVideo(vidObj1,scatframe);
        end

        % allmax = max(max(max(a1),max(a2)),max(max(b1),max(b2)));

        if branchedbundledvid==1
            % Concentric circles
            % Cell 1
            figcells=figure(2);
            clf
            alphaData=ZBranch1+max(0,max(max(ZBranch2))-max(max(ZBranch1)));
            surf(Xcol,Ycol1,ZBranch1,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
            colormap(branchedColor)
            clim([0,allmax/2])
            freezeColors;
            shading interp
            hold on;
            alphaData=ZBund1+max(0,max(max(ZBund2))-max(max(ZBund1)));
            surf(Xcol,Ycol1,ZBund1,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
            colormap(bundledColor)
            clim([0,allmax/2])
            freezeColors;
            shading interp
            view(2)
            grid off
            set(gca,'XTick',[], 'YTick', [])

            % Cell 2
            surf(Xcol,Ycol2,ZBranch2,'AlphaData',ZBranch2+max(0,max(max(ZBranch1))-max(max(ZBranch2))),'FaceAlpha','interp','FaceColor','interp');
            colormap(branchedColor)
            clim([0,allmax/2])
            freezeColors;
            cb=colorbar('Location','eastoutside');
            freezeColors(cb);
            cbpos=cb.Position;
            set(cb,'Position',[cbpos(1)+2*cbpos(3),cbpos(2),cbpos(3),cbpos(4)/2])
            % set(cb,'Position',[0.9062    0.1097    0.0235    0.4077])
            set(cb,'TickLabels',{});
            cbpos=cb.Position;
            shading interp
            % end
            % if max(ZBund2)>0.5
            surf(Xcol,Ycol2,ZBund2,'AlphaData',ZBund2+max(0,max(max(ZBund1))-max(max(ZBund2))),'FaceAlpha','interp','FaceColor','interp');
            colormap(bundledColor)
            clim([0,allmax/2])
            freezeColors;
            jcb=jicolorbar;
            freezeColors(jcb);
            jcbpos=jcb.Position;
            set(jcb,'Position',[cbpos(1)+cbpos(3),cbpos(2),cbpos(3),cbpos(4)])
            shading interp
            % end
            view(2)
            grid off
            axis equal
            set(gca,'XTick',[], 'YTick', [])

            hold off;
            box off;
            set(gca,'XColor','w')
            set(gca,'YColor','w')
            set(gcf,'color','w');
            title(strcat('t=',int2str(t)))
        end


        % Find median for cell 1
        a1New = a1;
        a1New(a1New<1)=0;
        if (a1New(1)~=0 && a1New(end)~=0)
            zeroInd1=find(a1New==0,1,'first');
            zeroInd2=find(a1New==0,1,'last');
            dirIndexa1=ceil((zeroInd1+zeroInd2)/2) - 50;
        else
            ind1=find(a1New~=0,1,'first');
            ind2=find(a1New~=0,1,'last');
            dirIndexa1=ceil((ind1+ind2)/2);
        end
        if dirIndexa1<1
            dirIndexa1=dirIndexa1+101;
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
        if dirIndexb1<1
            dirIndexb1=dirIndexb1+101;
        end
        if ~isempty(dirIndexa1) && ~isempty(dirIndexb1)
            hold on;
            quiver(0,0,Xsm(dirIndexa1),Ysm1(dirIndexa1),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5);
            hold off;
        end



        % Find median for cell 2
        a2New = a2;
        a2New(a2New<1)=0;
        if (a2New(1)~=0 && a2New(end)~=0)
            zeroInd1=find(a2New==0,1,'first');
            zeroInd2=find(a2New==0,1,'last');
            dirIndexa2=ceil((zeroInd1+zeroInd2)/2) - 50;
        else
            ind1=find(a2New~=0,1,'first');
            ind2=find(a2New~=0,1,'last');
            dirIndexa2=ceil((ind1+ind2)/2);
        end
        if dirIndexa2<1
            dirIndexa2=dirIndexa2+101;
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
        if dirIndexb2<1
            dirIndexb2=dirIndexb2+101;
        end
        if ~isempty(dirIndexa2) && ~isempty(dirIndexb2)
            hold on;
            quiver(0,-2*abs(max(max(Ycol2)))-gapsize,Xsm(dirIndexa2),Ysm2(dirIndexa2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5)
            hold off;
        end


        if branchedbundledvid==1
            % Plot signal
            if signal==1
                [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
                [Xsig,Ysig] = pol2cart(th,rad);
                hold on;
                scatter(Xsig(sigBound2),Ysig(sigBound2)-2*abs(max(max(Ycol2)))-gapsize,'black','.')
                hold off;
            end

            % figure(2)
            ohf = findobj(gcf);
            figaxes = findobj(ohf(1), 'Type', 'axes');
            set(figaxes(1),'Fontsize',15)
            set(figaxes(2),'Fontsize',14)
            camroll(90)


            cellsframe = getframe(figcells);
            writeVideo(vidObj2,cellsframe);
        end


        if racrhovid==1
            figracrho=figure(3);
            gapsize=0.05;
            clf
            hold on;
            th = (0:3.6:360)*pi/180;
            Xvals=cos(th);
            Yvals1=sin(th);
            Yvals1(boundC1)=Yvals1(boundC1(1));
            Yvals2=sin(th);
            Yvals2(boundC2)=Yvals2(boundC2(1));
            plot(Xvals,Yvals1,'black')
            plot(Xvals,Yvals2-2*abs(max(Yvals2))-gapsize,'black')
            YRac1=sin(posx1(1:NNx1(t),1)*2*pi/L);
            YRac1(YRac1<Yvals1(boundC1(1)))=Yvals1(boundC1(1));
            YRho1=sin(posy1(1:NNy1(t),1)*2*pi/L);
            YRho1(YRho1<Yvals1(boundC1(1)))=Yvals1(boundC1(1));
            YRac2=sin(posx2(1:NNx2(t),1)*2*pi/L)-gapsize;
            YRac2(YRac2>Yvals2(boundC2(1)))=Yvals2(boundC2(1))-gapsize;
            YRho2=sin(posy2(1:NNy2(t),1)*2*pi/L)-gapsize;
            YRho2(YRho2>Yvals2(boundC2(1)))=Yvals2(boundC2(1))-gapsize;
            scatter(cos(posx1(1:NNx1(t),1)*2*pi/L),YRac1,'MarkerEdgeColor',branchedColor(end,:),'linewidth',2)
            scatter(cos(posx2(1:NNx2(t),1)*2*pi/L),YRac2-2*abs(max(Yvals2)),'MarkerEdgeColor',branchedColor(end,:),'linewidth',2)
            scatter(cos(posy1(1:NNy1(t),1)*2*pi/L),YRho1,'MarkerEdgeColor',bundledColor(end,:),'linewidth',2)
            scatter(cos(posy2(1:NNy2(t),1)*2*pi/L),YRho2-2*abs(max(Yvals2)),'MarkerEdgeColor',bundledColor(end,:),'linewidth',2)
            hold off;
            set(gca,'XColor','w')
            set(gca,'YColor','w')
            set(gcf,'color','w')
            axis([-1 1 -3 1])
            axis equal

            if ~isempty(dirIndexa1)
                figure(3)
                hold on;
                quiver(0,0,Xsm(dirIndexa1),Ysm1(dirIndexa1),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5);
                hold off;
            end
            if ~isempty(dirIndexa2)
                figure(3)
                hold on;
                quiver(0,-2*abs(max(max(Ycol2)))-gapsize,Xsm(dirIndexa2),Ysm2(dirIndexa2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5)
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

            figure(3)
            camroll(90)
            racrhoframe = getframe(figracrho);
            writeVideo(vidObj3,racrhoframe);
        end



        if circlescatvid==1
            figcircle=figure(4);
            clf
            range=3;
            % subplot(1,2,1)
            plot(cos(2*pi*Xa/L),sin(2*pi*Xa/L),'color','black','LineWidth',2)
            hold on
            alphaData=ZBranch1+max(0,max(max(ZBranch2))-max(max(ZBranch1)));
            surf(Xcol,Ycol,ZBranch1,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
            % hold on
            colormap(branchedColor)
            clim([0,allmax/2])
            freezeColors;
            shading interp
            alphaData=ZBund1+max(0,max(max(ZBund2))-max(max(ZBund1)));
            surf(Xcol,Ycol,ZBund1,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
            colormap(bundledColor)
            clim([0,allmax/2])
            freezeColors;
            shading interp
            view(2)
            % plot(cos(2*pi*Xa/L),sin(2*pi*Xa/L),'color','black')
            % plot(2*cos(2*pi*Xa/L),2*sin(2*pi*Xa/L),'color','black')
            if max(xC1)>=(range+2)
                racxvals1=(range-1)*xC1/max(xC1)+1;
                racyvals1=(range-1)*xC1/max(xC1)+1;
            else
                racxvals1=(range-2)*xC1/max(xC1)+1;
                racyvals1=(range-1)*xC1/max(xC1)+1;
            end
            racxvals1=(racxvals1)'.*cos(2*pi*Xa/L);
            racyvals1=(racyvals1)'.*sin(2*pi*Xa/L);
            plot3(racxvals1,racyvals1,(allmax+1)*ones(1,length(racxvals1)),'color',...
                [branchedColor(end,:),1],'LineWidth',3)
            plot3([racxvals1(end),racxvals1(1)],[racyvals1(end),racyvals1(1)],...
                [allmax+1,allmax+1],'color',[branchedColor(end,:),1],'LineWidth',3)
            if max(yC1)>=(range+2)
                rhoxvals1=(range-1)*yC1/max(yC1)+1;
                rhoyvals1=(range-1)*yC1/max(yC1)+1;
            else
                rhoxvals1=(range-2)*yC1/max(yC1)+1;
                rhoyvals1=(range-2)*yC1/max(yC1)+1;
            end
            rhoxvals1=(rhoxvals1)'.*cos(2*pi*Xa/L);
            rhoyvals1=(rhoyvals1)'.*sin(2*pi*Xa/L);
            plot3(rhoxvals1,rhoyvals1,(allmax+1)*ones(1,length(rhoxvals1)),'color',...
                [bundledColor(end,:),1],'LineWidth',3)
            plot3([rhoxvals1(end),rhoxvals1(1)],[rhoyvals1(end),rhoyvals1(1)],...
                [allmax+1,allmax+1],'color',[bundledColor(end,:),1],'LineWidth',3)
            % xlim([-3,3])
            % ylim([-3,3])
            % axis square
            hold off

            %cell 2
            hold on
            if adjacent==0
                plot(cos(2*pi*Xa/L),sin(2*pi*Xa/L)-(2*range),'color','black','LineWidth',2)
            else
                plot(cos(2*pi*Xa/L),sin(2*pi*Xa/L)-2-(range-1),'color','black','LineWidth',2)
            end
            alphaData=ZBranch2+max(0,max(max(ZBranch1))-max(max(ZBranch2)));
            if adjacent==0
                surf(Xcol,Ycol-(2*range),ZBranch2,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
            else
                surf(Xcol,Ycol-2-(range-1),ZBranch2,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
            end
            % hold on
            colormap(branchedColor)
            clim([0,allmax/2])
            cb=colorbar('Location','eastoutside');
            freezeColors;
            freezeColors(cb);
            cbpos=cb.Position;
            % set(cb,'Position',[cbpos(1)+2*cbpos(3),cbpos(2),cbpos(3),cbpos(4)/2])
            set(cb,'Position',[0.9062    0.1097    0.0235    0.4077])
            set(cb,'TickLabels',{});
            cbpos=cb.Position;
            shading interp
            alphaData=ZBund2+max(0,max(max(ZBund1))-max(max(ZBund2)));
            if adjacent==0
                surf(Xcol,Ycol-(2*range),ZBund2,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
            else
                surf(Xcol,Ycol-2-(range-1),ZBund2,'AlphaData',alphaData,'FaceAlpha','interp','FaceColor','interp');
            end
            colormap(bundledColor)
            clim([0,allmax/2])
            freezeColors;
            jcb=jicolorbar;
            freezeColors(jcb);
            jcbpos=jcb.Position;
            set(jcb,'Position',[cbpos(1)+cbpos(3),cbpos(2),cbpos(3),cbpos(4)])
            % set(cb,'Position',[cbpos(1),cbpos(2),cbpos(3),jcbpos(4)])
            shading interp
            view(2)
            % plot(cos(2*pi*Xa/L),sin(2*pi*Xa/L),'color','black')
            % plot(2*cos(2*pi*Xa/L),2*sin(2*pi*Xa/L),'color','black')
            if max(xC2)>=(range+2)
                racxvals2=(range-1)*xC2/max(xC2)+1;
                racyvals2=(range-1)*xC2/max(xC2)+1;
            else
                racxvals2=(range-2)*xC2/max(xC2)+1;
                racyvals2=(range-2)*xC2/max(xC2)+1;
            end
            racxvals2=(racxvals2)'.*cos(2*pi*Xa/L);
            racyvals2=(racyvals2)'.*sin(2*pi*Xa/L);
            if adjacent==0
                plot3(racxvals2,racyvals2-(2*range),(allmax+1)*ones(1,length(racxvals2)),'color',...
                    branchedColor(end,:),'LineWidth',3)
                plot3([racxvals2(end),racxvals2(1)],[racyvals2(end),racyvals2(1)]-(2*range),...
                    [allmax+1,allmax+1],'color',branchedColor(end,:),'LineWidth',3)
            else
                plot3(racxvals2,racyvals2-2-(range-1),(allmax+1)*ones(1,length(racxvals2)),'color',...
                    [branchedColor(end,:),0.5],'LineWidth',3)
                plot3([racxvals2(end),racxvals2(1)],[racyvals2(end),racyvals2(1)]-2-(range-1),...
                    [allmax+1,allmax+1],'color',[branchedColor(end,:),0.5],'LineWidth',3)
            end
            if max(yC2)>=(range+2)
                rhoxvals2=(range-1)*yC2/max(yC2)+1;
                rhoyvals2=(range-1)*yC2/max(yC2)+1;
            else
                rhoxvals2=(range-2)*yC2/max(yC2)+1;
                rhoyvals2=(range-2)*yC2/max(yC2)+1;
            end
            rhoxvals2=(rhoxvals2)'.*cos(2*pi*Xa/L);
            rhoyvals2=(rhoyvals2)'.*sin(2*pi*Xa/L);
            if adjacent==0
                plot3(rhoxvals2,rhoyvals2-(2*range),(allmax+1)*ones(1,length(rhoxvals2)),'color',...
                    bundledColor(end,:),'LineWidth',3)
                plot3([rhoxvals2(end),rhoxvals2(1)],[rhoyvals2(end),rhoyvals2(1)]-2*range,...
                    [allmax+1,allmax+1],'color',bundledColor(end,:),'LineWidth',3)
            else
                plot3(rhoxvals2,rhoyvals2-2-(range-1),(allmax+1)*ones(1,length(rhoxvals1)),'color',...
                    [bundledColor(end,:),0.5],'LineWidth',3)
                plot3([rhoxvals2(end),rhoxvals2(1)],[rhoyvals2(end),rhoyvals2(1)]-2-(range-1),...
                    [allmax+1,allmax+1],'color',[bundledColor(end,:),0.5],'LineWidth',3)
            end
            if adjacent==0
                xlim([-3,3])
                ylim([-9,3])
                set(gca,'plotBoxAspectRatio',[1 2 1]);
            else
                % xlim([-3,3])
                % ylim([-5,3])
                % set(gca,'plotBoxAspectRatio',[6 8 1]);
                xlim([-3,3])
                ylim([-7,3])
                set(gca,'plotBoxAspectRatio',[6 10 1]);
            end
            % axis square
            hold off


            timebox=annotation('textbox', [0.75, 0.1, 0.1, 0.1], 'String', "t = " + t,'FitBoxToText','on');
            tbpos=timebox.Position;
            set(timebox,'Position',[cbpos(1)-tbpos(3), cbpos(2), 0.1, 0.1]);

            if ~isempty(dirIndexa1)
                figure(4)
                hold on;
                quiver(0,0,Xsm(dirIndexa1),Ysm(dirIndexa1),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.7);
                hold off;
            end
            if ~isempty(dirIndexa2)
                figure(4)
                hold on;
                if adjacent==0
                    quiver(0,-2*range,Xsm(dirIndexa2),Ysm(dirIndexa2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.7)
                else
                    quiver(0,-2-(range-1),Xsm(dirIndexa2),Ysm(dirIndexa2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.7)
                end
                hold off;
            end


            if signal==1
                [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
                [Xsig,Ysig] = pol2cart(th,rad);
                if t<=500
                    hold on;
                    scatter(Xsig(sigBound2),Ysig(sigBound2)-2-(range-1),'black','.')
                    hold off;
                else
                    hold on;
                    scatter(Xsig(sigBound1),Ysig(sigBound1),'black','.')
                    hold off;
                end
            end

            grid off
            set(gca,'XTick',[],'YTick',[])
            set(gca,'XColor','w')
            set(gca,'YColor','w')
            set(gcf,'color','w');
            set(gcf,'Position',[209   561   682   474])
            ohf = findobj(gcf);
            figaxes = findobj(ohf(1), 'Type', 'axes');
            set(figaxes(1),'Fontsize',15)
            set(figaxes(2),'Fontsize',14)
            camroll(90)

            circleframe = getframe(figcircle);
            writeVideo(vidObj4,circleframe);

        end

    end

    if scatvid==1
        close(vidObj1);
    end
    if branchedbundledvid==1
        close(vidObj2);
    end
    if racrhovid==1
        close(vidObj3);
    end
    if circlescatvid==1
        close(vidObj4);
    end

end