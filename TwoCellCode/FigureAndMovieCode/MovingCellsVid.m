addpath('./freeze_colors')
set(0,'DefaultFigureVisible','off')
close all;
clear;
clc;

signal=0;
squished=0;
add_shift=1;
timeskip=50;

Nt=2500;

for i=1:1

    loadfile='./vid_matfiles/leader_follower/racupc2/1000RacOn';

    setnum=int2str(i);
    savelocation='../../movies_for_paper/celldoublet_branchedbundledvid_1000RacOnC1_supracellular';

    % vidObj1 = VideoWriter(strcat(savelocation,'ScatterVid_',setnum,'.mp4'),'MPEG-4');
    vidObj2 = VideoWriter(strcat(savelocation,'_BranchedBundledVid_',setnum,'.mp4'),'MPEG-4');
    % vidObj3 = VideoWriter(strcat(savelocation,'RacRhoVid_',setnum,'.mp4'),'MPEG-4');

    % vidObj1.FrameRate = 5;
    % vidObj1.Quality = 100;
    vidObj2.FrameRate = 5;
    vidObj2.Quality = 100;
    % vidObj3.FrameRate = 5;
    % vidObj3.Quality = 100;

    % open(vidObj1);
    open(vidObj2);
    % open(vidObj3);

    load(strcat(loadfile,setnum,'.mat'));
    set(0,'DefaultFigureVisible','off')

    linewidth=1;
    linecolor=[0 0 0];
    showtime=1;

    allmax = max(max(max(max(a1all)),max(max(a2all))),max(max(max(b1all)),max(max(b2all))));

    if add_shift==1
        xshift1=zeros(1,Nt);
        yshift1=zeros(1,Nt);
        xshift2=zeros(1,Nt);
        yshift2=zeros(1,Nt);
    end


    for t=1:timeskip:Nt-1
        a1=a1all(:,t);
        a2=a2all(:,t);
        b1=b1all(:,t);
        b2=b2all(:,t);
        L      = 10.0;
        Na=101;

        % [s1,xC1,yC1] = resamplePolarityMolecules(posx1(1:NNx1(t+1),t+1),posy1(1:NNy1(t+1),t+1),NNx1(t+1),NNy1(t+1),L,Na);
        % [s2,xC2,yC2] = resamplePolarityMolecules(posx2(1:NNx2(t+1),t+1),posy2(1:NNy2(t+1),t+1),NNx2(t+1),NNy2(t+1),L,Na);

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
        
        if add_shift==1
            [th,~] = meshgrid((0:3.6:360)*pi/180,1);
            xshift1(t+1:t+timeskip)=xshift1(t);
            yshift1(t+1:t+timeskip)=yshift1(t);
            xshift2(t+1:t+timeskip)=xshift2(t);
            yshift2(t+1:t+timeskip)=yshift2(t);

            if ~isempty(dirIndexa1) && ~isempty(dirIndexb1)
                xshift1(t+1:t+timeskip)=xshift1(t+1:t+timeskip)+cos(th(dirIndexa1))*0.001*timeskip;
                yshift1(t+1:t+timeskip)=yshift1(t+1:t+timeskip)+sin(th(dirIndexa1))*0.001*timeskip;
            end
            if ~isempty(dirIndexa2) && ~isempty(dirIndexb2)
                xshift2(t+1:t+timeskip)=xshift2(t+1:t+timeskip)+cos(th(dirIndexa2))*0.001*timeskip;
                yshift2(t+1:t+timeskip)=yshift2(t+1:t+timeskip)+sin(th(dirIndexa2))*0.001*timeskip;
            end
        end
        

        % Define circles
        gapsize=0.01 * squished;
        [th,rad] = meshgrid((0:3.6:360)*pi/180,0.85:0.01:1);
        [Xcol,Ycol1] = pol2cart(th,rad);
        [~,Ycol2] = pol2cart(th,rad);
        if squished==1
            Ycol1(:,boundC1)=Ycol1(:,boundC1(1)*ones(1,length(boundC1)));
            Ycol2(:,boundC2)=Ycol2(:,boundC2(1)*ones(1,length(boundC2)));
        end
        Ycol2 = Ycol2 - 2*abs(max(max(Ycol2)))-gapsize;
        ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1 a1]';
        ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1 b1]';
        ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2 a2]';
        ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2 b2]';
        [th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
        [Xsm,Ysm1] = pol2cart(th,rad);
        [~,Ysm2] = pol2cart(th,rad);
        if squished==1
            Ysm1(:,boundC1)=Ysm1(:,boundC1(1)*ones(1,length(boundC1)));
            Ysm2(:,boundC2)=Ysm2(:,boundC2(1)*ones(1,length(boundC2)));
        end


        sigper=0.40;
        sigBound1 = (floor((Na-1)*3/8 - floor((Na-1)*sigper/2)))+1:(floor((Na-1)*3/8 + floor((Na-1)*sigper/2)))+1;
        sigBound2 = (floor((Na-1)*5/8 - floor((Na-1)*sigper/2)))+1:(floor((Na-1)*5/8 + floor((Na-1)*sigper/2)))+1;

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

        % Cell 1
        figcells=figure(1);
        clf
        hold on;
        if squished==1
            plot3(cos(th(1,[1:boundC1(1),boundC1(end):end]))+xshift1(t+1),sin(th(1,[1:boundC1(1),boundC1(end):end]))+yshift1(t+1),ones(1,length(th(1,[1:boundC1(1),boundC1(end):end])))*(allmax+1),'color',linecolor,'LineWidth',linewidth)
        else
            plot3(cos(th(1,:))+xshift1(t+1),sin(th(1,:))+yshift1(t+1),ones(1,length(th))*(allmax+1),'color',linecolor,'LineWidth',linewidth)
        end
        surf(Xcol+xshift1(t+1),Ycol1+yshift1(t+1),ZBranch1,'AlphaData',ZBranch1+max(0,max(max(ZBranch2))-max(max(ZBranch1))),'FaceAlpha','interp','FaceColor','interp');
        view(2)
        clim([0,12])
        colormap(branchedColor)
        freezeColors;
        % freezeColors(colorbar('Location','eastoutside'));
        clim([0,12])
        shading interp
        surf(Xcol+xshift1(t+1),Ycol1+yshift1(t+1),ZBund1,'AlphaData',ZBund1+max(0,max(max(ZBund2))-max(max(ZBund1))),'FaceAlpha','interp','FaceColor','interp');
        colormap(bundledColor)
        clim([0,12])
        freezeColors;
        % freezeColors(jicolorbar);
        shading interp
        grid off
        % set(gca,'XTick',[], 'YTick', [])

        % Cell 2
        if squished==1
            plot3(cos(th(1,[1:boundC2(1),boundC2(end):end]))+xshift2(t+1),sin(th(1,[1:boundC2(1),boundC2(end):end]))-2*abs(max(max(Ycol2)))-gapsize+yshift2(t+1),ones(1,length(th(1,[1:boundC2(1),boundC2(end):end])))*(allmax+1),'color',linecolor,'LineWidth',linewidth)
        else
            plot3(cos(th(1,:))+xshift2(t+1),sin(th(1,:))-2+yshift2(t+1),ones(1,length(th))*(allmax+1),'color',linecolor,'LineWidth',linewidth)
        end
        % plot3(cos(th(1,:))+xshift2(t+1),sin(th(1,:))-2+yshift2(t+1),ones(1,length(th))*(allmax+1),'color',[0.5,0.5,0.5],'LineWidth',1)
        surf(Xcol+xshift2(t+1),Ycol2+yshift2(t+1),ZBranch2,'AlphaData',ZBranch2+max(0,max(max(ZBranch1))-max(max(ZBranch2))),'FaceAlpha','interp','FaceColor','interp');
        view(2)
        colormap(branchedColor)
        clim([0,12])
        cb=colorbar('Location','eastoutside');
        freezeColors;
        freezeColors(cb);
        cbpos=cb.Position;
        % set(cb,'Position',[cbpos(1)+2*cbpos(3),cbpos(2),cbpos(3),cbpos(4)/2])
        set(cb,'Position',[0.9062    0.1097    0.0235    0.4077])
        set(cb,'TickLabels',{});
        cbpos=cb.Position;
        shading interp
        surf(Xcol+xshift2(t+1),Ycol2+yshift2(t+1),ZBund2,'AlphaData',ZBund2+max(0,max(max(ZBund1))-max(max(ZBund2))),'FaceAlpha','interp','FaceColor','interp');
        colormap(bundledColor)
        clim([0,12])
        freezeColors;
        jcb=jicolorbar;
        freezeColors(jcb);
        jcbpos=jcb.Position;
        set(jcb,'Position',[cbpos(1)+cbpos(3),cbpos(2),cbpos(3),cbpos(4)])
        shading interp
        grid off
        xlim([-4,4])
        ylim([-6,2])
        pbaspect([8 8 1])
        % axis square
        set(gca,'XTick',[], 'YTick', [])

        if showtime==1
            % title(strcat('t=',int2str(t)))
            timebox=annotation('textbox', [0.75, cbpos(2), 0.1, 0.05], 'String', "t = " + (t-1)*0.01,'FitBoxToText','on','EdgeColor','none','FontSize',20);
            % tbpos=timebox.Position;
            % set(timebox,'Position',[cbpos(1)-tbpos(3), cbpos(2), 0.1, 0.05]);
        end

        hold off;
        % box off;
        set(gca,'XColor','w')
        set(gca,'YColor','w')
        set(gcf,'color','w');


       
        if ~isempty(dirIndexa1)
            hold on;
            quiver(0+xshift1(t+1),0+yshift1(t+1),Xsm(dirIndexa1),Ysm1(dirIndexa1),0,'color',linecolor,'LineWidth',2,'MaxHeadSize',2);
            hold off;
        end
        if ~isempty(dirIndexa2)
            hold on;
            quiver(0+xshift2(t+1),-2*abs(max(max(Ycol2)))-gapsize+yshift2(t+1),Xsm(dirIndexa2),Ysm2(dirIndexa2),0,'color',linecolor,'LineWidth',2,'MaxHeadSize',2)
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

        % figure(2)
        ohf = findobj(gcf);
        figaxes = findobj(ohf(1), 'Type', 'axes');
        set(figaxes(1),'Fontsize',15)
        set(figaxes(2),'Fontsize',14)
        camroll(90)


        cellsframe = getframe(figcells);
        writeVideo(vidObj2,cellsframe);

    end

    % close(vidObj1);
    close(vidObj2);
    % close(vidObj3);

end