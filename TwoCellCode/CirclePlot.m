%Plot on circle
%Define colors
colorLength = 50;
white = [1,1,1];
red = [1,0,0];
blue = [0,0.5,1];
maroon = [0.4,0,0];
navy = [0,0,0.3];
whitered = [linspace(white(1),red(1),colorLength)',linspace(white(2),red(2),colorLength)',linspace(white(3),red(3),colorLength)'];
redmaroon = [linspace(red(1),maroon(1),colorLength)',linspace(red(2),maroon(2),colorLength)',linspace(red(3),maroon(3),colorLength)'];
whiteredmaroon = [whitered;redmaroon];
whiteblue = [linspace(white(1),blue(1),colorLength)',linspace(white(2),blue(2),colorLength)',linspace(white(3),blue(3),colorLength)'];
bluenavy = [linspace(blue(1),navy(1),colorLength)',linspace(blue(2),navy(2),colorLength)',linspace(blue(3),navy(3),colorLength)'];
whitebluenavy = [whiteblue; bluenavy];
myColors = [linspace(red(1),blue(1),colorLength)',linspace(red(2),blue(2),colorLength)',linspace(red(3),blue(3),colorLength)'];
redblue = abs(whiteblue+whitered)./2;
redwhiteblue = [flip(whitered); whiteblue];

[th,rad] = meshgrid((0:3.6:360)*pi/180,0.93:0.01:1);
[Xcol,Ycol] = pol2cart(th,rad);
ZBranch1 = [a1 a1 a1 a1 a1 a1 a1 a1]';
ZBund1 = [b1 b1 b1 b1 b1 b1 b1 b1]';

ZBranch2 = [a2 a2 a2 a2 a2 a2 a2 a2]';
ZBund2 = [b2 b2 b2 b2 b2 b2 b2 b2]';

% Concentric circles
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.8);
[Xsm,Ysm] = pol2cart(th,rad);
[th,rad] = meshgrid((0:3.6:360)*pi/180,0.83:0.01:0.9);
[Xmid,Ymid] = pol2cart(th,rad);

% Cell 1
% figc1=figure(16);
figcells=figure(16);
subplot(2,1,1)
surf(Xcol,Ycol,ZBranch1);
view(2)
colormap(whitebluenavy)
freezeColors;
freezeColors(colorbar);
clim([0,max(b1)])
shading interp
hold on;
surf(Xmid,Ymid,ZBund1);
colormap(whiteredmaroon)
freezeColors;
freezeColors(jicolorbar);
clim([0,max(a1)])
shading interp
grid off
axis square
set(gca,'XTick',[], 'YTick', [])
title('Cell 1 Combined: Blue=Branched, Red=Bundled')
scatter(Xsm(boundC1),Ysm(boundC1),'black');
hold off;
box on;
box off;
set(gca,'XColor','w')
set(gca,'YColor','w')
set(gcf,'color','w');

[Ma1,Ia1] = max(a1);
[Mb1,Ib1] = max(b1);

hold on;
quiver(0,0,Xsm(1,Ia1),Ysm(1,Ia1),0,'color',[0 0 0])
% plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
% plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
% plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
hold off;

% sinsum1=0;
% cossum1=0;
% for i=1:length(a1)
%     sinsum1=sinsum1 + a1(i)*sin(th(i));
%     cossum1=cossum1 + a1(i)*cos(th(i));
% end
% sinsum2=0;
% cossum2=0;
% for i=1:length(a2)
%     sinsum2=sinsum2 + a2(i)*sin(th(i));
%     cossum2=cossum2 + a2(i)*cos(th(i));
% end
% avg1=atan2(sinsum1,cossum1); %atan(sinsum1/cossum1);
% avg2=atan2(sinsum2,cossum2); %atan(sinsum2/cossum2);

% hold on;
% quiver(0,0,0.8*cos(avg1),0.8*sin(avg1),0,'color',[1 0 0])
% hold off;


a1New = a1;
a1New(a1New<1)=0;
% th2=th(1,:);
% if (a1New(1)~=0 && a1New(length(a1New))~=0)
%     th2(ceil(length(th2)/2):length(th2))=-flip(th2(1:ceil(length(th2)/2)));
%     % zeroIndex=find(a1New==0,1,'first');
%     % a1New(zeroIndex:length(a1New))=-a1New(zeroIndex:length(a1New));
% end
% avg1 = sum(a1New.*th(1,:)')/length(a1New);
% avg1=sum(a1New.*th2')/length(a1New);
% hold on;
% quiver(0,0,0.8*cos(avg1),0.8*sin(avg1),0,'color',[1 0 0])
% hold off;

if (a1New(1)~=0 && a1New(length(a1New))~=0)
    zeroInd1=find(a1New==0,1,'first');
    zeroInd2=find(a1New==0,1,'last');
    dirIndex1=ceil((zeroInd1+zeroInd2)/2) - 50;
else
    ind1=find(a1New~=0,1,'first');
    ind2=find(a1New~=0,1,'last');
    dirIndex1=ceil((ind1+ind2)/2);
end
if isempty(dirIndex1)
    dirIndex1=1;
end
if dirIndex1<1
    dirIndex1=dirIndex1+101;
end
hold on;
quiver(0,0,Xsm(dirIndex1),Ysm(dirIndex1),0,'color',[1 0 0])
hold off;

% Cell 2
% figc2=figure(17);
% clf
subplot(2,1,2)
surf(Xcol,Ycol,ZBranch2);
view(2)
colormap(whitebluenavy)
freezeColors;
freezeColors(colorbar);
clim([0,max(b2)])
shading interp
hold on;
surf(Xmid,Ymid,ZBund2);
colormap(whiteredmaroon)
freezeColors;
freezeColors(jicolorbar);
clim([0,max(a2)])
shading interp
grid off
axis square
set(gca,'XTick',[], 'YTick', [])
title('Cell 2 Combined: Blue=Branched, Red=Bundled')
scatter(Xsm(boundC2),Ysm(boundC2),'black');
hold off;
box off;
set(gca,'XColor','w')
set(gca,'YColor','w')
set(gcf,'color','w');

[Ma2,Ia2] = max(a2);
[Mb2,Ib2] = max(b2);
p1=[0 0];
p2=[Xsm(1,Ia2) Ysm(1,Ia2)];
dp=p2-p1;

hold on;
quiver(p1(1),p1(2),dp(1),dp(2),0,'color',[0 0 0])
% plot([Xcol(3,Ia) Xcol(3,Ib)],[Ycol(3,Ia) Ycol(3,Ib)],"black")
% plot([0 Xcol(3,Ia)],[0 Ycol(3,Ia)],"black")
% plot([0 Xcol(3,Ib)],[0 Ycol(3,Ib)],"black")
hold off;

a2New = a2;
a2New(a2New<1)=0;
% th2=th(1,:);
% if (a2New(1)~=0 && a2New(length(a1New))~=0)
%     th2(ceil(length(th2)/2):length(th2))=-flip(th2(1:ceil(length(th2)/2)));
%     % zeroIndex=find(a1New==0,1,'first');
%     % a1New(zeroIndex:length(a1New))=-a1New(zeroIndex:length(a1New));
% end
% avg2 = sum(a2New.*th(1,:)')/length(a2New);
% avg2=sum(a2New.*th2')/length(a2New);
% hold on;
% quiver(0,0,0.8*cos(avg2),0.8*sin(avg2),0,'color',[1 0 0])
% hold off;

if (a2New(1)~=0 && a2New(length(a2New))~=0)
    zeroInd1=find(a2New==0,1,'first');
    zeroInd2=find(a2New==0,1,'last');
    dirIndex2=ceil((zeroInd1+zeroInd2)/2) - 50;
else
    ind1=find(a2New~=0,1,'first');
    ind2=find(a2New~=0,1,'last');
    dirIndex2=ceil((ind1+ind2)/2);
end
if isempty(dirIndex2)
    dirIndex2=1;
end
if dirIndex2<1
    dirIndex2=dirIndex2+101;
end
hold on;
quiver(0,0,Xsm(dirIndex2),Ysm(dirIndex2),0,'color',[1 0 0])
hold off;