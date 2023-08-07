set(0,'DefaultFigureVisible','on')
close all;

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

% Concentric circles
% Cell 1
figcells=figure(16);
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

% Find peak for cell 1
N=length(a1);
% smoothed=ones(1,N);
% for i=1:N
%     smoothed(i) = (1/7)*(a1((i-3)*(i>3)+(i-3+N)*(i<=3)) ...
%         +a1((i-2)*(i>2)+(i-2+N)*(i<=2)) ...
%         +a1((i-1)*(i>1)+(N)*(i==1)) ...
%         +a1(i) ...
%         +a1((i+1)*(i<N)+(1)*(i==N)) ...
%         +a1((i+2)*(i<N-1)+(i+2-N)*(i>=N-1)) ...
%         +a1((i+3)*(i<N-2)+(i+3-N)*(i>=N-2)));
% end

epsilon=0.01;                     % variance squared (default)
x=linspace(0,L,Na);
a0=sqrt(pi)/(sqrt(2*pi*epsilon)); % area of normal distribution
smoothed=zeros(N,1);
for i=1:N
    smoothed(i) = exp(-(x-a1(i)).^2/(2*epsilon))./(sqrt(2*pi*epsilon));
end
smoothed=smoothed/a0;
% smoothed=sum(smoothed,2);

baseline=0.1;
peakIndices1=[];
peakIndex=0;
peakValue=0;
for i=1:N
    if smoothed(i)>baseline
        if peakValue==0 || smoothed(i)>peakValue
            peakIndex=i;
            peakValue=smoothed(i);
        end
    elseif smoothed(i)<baseline && peakIndex~=0
        peakIndices1=[peakIndices1,peakIndex];
        peakIndex=0;
        peakValue=0;
    end
end
if peakIndex~=0
    peakIndices1=[peakIndices1,peakIndex];
end

[maxval,maxind] = max(smoothed(peakIndices1));
peakInd1 = peakIndices1(maxind);
figure(16)
hold on
quiver(0,0,Xsm(peakInd1),Ysm(peakInd1),0,'color',[1 0 0],'LineWidth',2,'MaxHeadSize',0.5)
hold off;

% Find peak for cell 2
N=length(a2);
smoothed=ones(1,N);
for i=1:N
    smoothed(i) = (1/5)*(a2((i-2)*(i>2)+(i-2+N)*(i<=2))+a2((i-1)*(i>1)+(N)*(i==1))+a2(i)+a2((i+1)*(i<N)+(1)*(i==N))+a2((i+2)*(i<N-1)+(i+2-N)*(i>=N-1)));
end
baseline=0.1;
peakIndices2=[];
peakIndex=0;
peakValue=0;
for i=1:N
    if smoothed(i)>baseline
        if peakValue==0 || smoothed(i)>peakValue
            peakIndex=i;
            peakValue=smoothed(i);
        end
    elseif smoothed(i)<baseline && peakIndex~=0
        peakIndices2=[peakIndices2,peakIndex];
        peakIndex=0;
        peakValue=0;
    end
end
if peakIndex~=0
    peakIndices2=[peakIndices2,peakIndex];
end
[maxval,maxind] = max(smoothed(peakIndices2));
peakInd2 = peakIndices2(maxind);
figure(16)
hold on
quiver(0,-2,Xsm(peakInd2),Ysm(peakInd2),0,'color',[1 0 0],'LineWidth',2,'MaxHeadSize',0.5)
hold off;

% Plot signal
figure(16)
[th,rad] = meshgrid((0:3.6:360)*pi/180,1.1);
[Xsig,Ysig] = pol2cart(th,rad);
hold on;
scatter(Xsig(sigBound),Ysig(sigBound)-2,'black','.')
hold off;

ohf = findobj(gcf);
figaxes = findobj(ohf(1), 'Type', 'axes');
set(figaxes(1),'Fontsize',15)
set(figaxes(2),'Fontsize',14)
camroll(90)


figure(2)
subplot(1,2,1)
hold on
plot(Xa,a1,'markerfacecolor',[159 219 229]/255,'linewidth',2); hold on;
plot(Xa,b1,'markerfacecolor','k','linewidth',2);
scatter(Xa(peakIndices1),a1(peakIndices1),'filled')
lgd = legend('Branched','Bundled','Location','northeast');
lgd.NumColumns = 1;
hold off
subplot(1,2,2)
hold on
plot(Xa,a2,'markerfacecolor',[159 219 229]/255,'linewidth',2); hold on;
plot(Xa,b2,'markerfacecolor','k','linewidth',2);
scatter(Xa(peakIndices2),a2(peakIndices2),'filled')
lgd = legend('Branched','Bundled','Location','northeast');
lgd.NumColumns = 1;
hold off
