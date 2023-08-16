set(0,'DefaultFigureVisible','on')
close all;
clear;
clc;

signal=1;

loadfile='./sigswitch_vid_files/vars_t';

savelocation='./movies4/sigswitch_uncoupled';
setnum=int2str(1);

vidObj1 = VideoWriter(strcat(savelocation,'ScatterVid_',setnum,'.mp4'),'MPEG-4');
vidObj2 = VideoWriter(strcat(savelocation,'ColorVid_',setnum,'.mp4'),'MPEG-4');

vidObj1.FrameRate = 5;
vidObj1.Quality = 75;
vidObj2.FrameRate = 5;
vidObj2.Quality = 75;

open(vidObj1);
open(vidObj2);

allmax=0;
for t=1:44
    load(strcat(loadfile,int2str(t*50)));
    allmax = max(max(max(max(a1),max(a2)),max(max(b1),max(b2))),allmax);
end

% Define circles
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.93:0.01:1);
[Xcol,Ycol] = pol2cart(th,rad);
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
[Xsm,Ysm] = pol2cart(th,rad);
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.86:0.01:0.93);
[Xmid,Ymid] = pol2cart(th,rad);

for t=1:49
    close all;
    clearvars -except t vidObj1 vidObj2 Xcol Ycol Xsm Ysm Xmid Ymid allmax loadfile savelocation signal sigBound1 sigBound2
    load(strcat(loadfile,int2str(t*50)));

    Na=101;
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

    ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1]';
    ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1]';
    ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2]';
    ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2]';

    % allmax = max(max(max(a1),max(a2)),max(max(b1),max(b2)));

    % Concentric circles
    % Cell 1
    figcells=figure(2);
    clf
    % if max(ZBranch1)>0.5
    surf(Xcol,Ycol,ZBranch1,'AlphaData',ZBranch1+max(0,max(max(ZBranch2))-max(max(ZBranch1))),'FaceAlpha','interp','FaceColor','interp');
    colormap(branchedColor)
    freezeColors;
    freezeColors(colorbar('Location','westoutside'));
    clim([0,allmax])
    shading interp
    % end
    hold on;
    % if max(ZBund1)>0.5
    surf(Xcol,Ycol,ZBund1,'AlphaData',ZBund1+max(0,max(max(ZBund2))-max(max(ZBund1))),'FaceAlpha','interp','FaceColor','interp');
    colormap(bundledColor)
    freezeColors;
    freezeColors(jicolorbar);
    clim([0,allmax])
    shading interp
    % end
    view(2)
    grid off
    set(gca,'XTick',[], 'YTick', [])

    % Cell 2
    % if max(ZBranch2)>0.5
    surf(Xcol,Ycol-2,ZBranch2,'AlphaData',ZBranch2+max(0,max(max(ZBranch1))-max(max(ZBranch2))),'FaceAlpha','interp','FaceColor','interp');
    colormap(branchedColor)
    freezeColors;
    freezeColors(colorbar('Location','westoutside'));
    clim([0,allmax])
    shading interp
    % end
    % if max(ZBund2)>0.5
    surf(Xcol,Ycol-2,ZBund2,'AlphaData',ZBund2+max(0,max(max(ZBund1))-max(max(ZBund2))),'FaceAlpha','interp','FaceColor','interp');
    colormap(bundledColor)
    freezeColors;
    freezeColors(jicolorbar);
    clim([0,allmax])
    shading interp
    % end
    view(2)
    grid off
    axis equal
    xlim([-1.5,1.5]);
    ylim([-3.5,1.5]);
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


    % Find peak for cell 1
    N=length(a1);
    smoothed1=zeros(1,N);
    smoothAmount=20;
    region=-smoothAmount:smoothAmount;
    for i=1:N
        b=i-region;
        smoothed1(i) = sum(a1(b.*(b>0 & b<=N)+(b+N).*(b<=0)+(b-N).*(b>N)))/(2*smoothAmount+1);
    end

    baseline=1;
    peakIndices1=[];
    peakIndex=0;
    peakValue=0;
    for i=1:N
        if smoothed1(i)>baseline
            if peakValue==0 || smoothed1(i)>peakValue
                peakIndex=i;
                peakValue=smoothed1(i);
            end
        elseif smoothed1(i)<baseline && peakIndex~=0
            peakIndices1=[peakIndices1,peakIndex];
            peakIndex=0;
            peakValue=0;
        end
    end
    if peakIndex~=0
        peakIndices1=[peakIndices1,peakIndex];
    end

    [~,maxind] = max(smoothed1(peakIndices1));
    peakInd1 = peakIndices1(maxind);
    if ~isempty(peakInd1) && min(a1)<0.01 && max(b1)>1
        % figure(2)
        hold on
        quiver(0,0,Xsm(peakInd1),Ysm(peakInd1),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5)
        hold off;
    end

    % Find peak for cell 2
    N=length(a2);
    smoothed2=zeros(1,N);
    region=-smoothAmount:smoothAmount;
    for i=1:N
        b=i-region;
        smoothed2(i) = sum(a2(b.*(b>0 & b<=N)+(b+N).*(b<=0)+(b-N).*(b>N)))/(2*smoothAmount+1);
    end

    peakIndices2=[];
    peakIndex=0;
    peakValue=0;
    for i=1:N
        if smoothed2(i)>baseline
            if peakValue==0 || smoothed2(i)>peakValue
                peakIndex=i;
                peakValue=smoothed2(i);
            end
        elseif smoothed2(i)<baseline && peakIndex~=0
            peakIndices2=[peakIndices2,peakIndex];
            peakIndex=0;
            peakValue=0;
        end
    end
    if peakIndex~=0
        peakIndices2=[peakIndices2,peakIndex];
    end
    [maxval,maxind] = max(smoothed2(peakIndices2));
    peakInd2 = peakIndices2(maxind);
    if ~isempty(peakInd2) && min(a2)<0.01 && max(b2)>1
        % figure(2)
        hold on
        quiver(0,-2,Xsm(peakInd2),Ysm(peakInd2),0,'color',[0 0 0],'LineWidth',2,'MaxHeadSize',0.5)
        hold off;
    end

    % Plot signal
        if signal==1
            [th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
            [Xsig,Ysig] = pol2cart(th,rad);
            if t*50<=1500
                hold on;
                scatter(Xsig(sigBound2),Ysig(sigBound2)-2,'black','.')
                hold off;
            else
                hold on;
                scatter(Xsig(sigBound1),Ysig(sigBound1),'black','.')
                hold off;
            end
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

close(vidObj1);
close(vidObj2);