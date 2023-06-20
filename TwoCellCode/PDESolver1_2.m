clf;

Da      = 0.5;                  % diffusion coefficient for actin
m0      = 2.0;                  % competition for actin monomers
K       = 1.0;  

Tend    = 25.0;
dt      = 0.1; %0.05
Nt      = Tend/dt; % number of time steps
Na      = 101; % number of space steps
dxa     = 5.0/((Na-1)/2);
L       = 10.0;
Xa      = 0:dxa:L;
Xb      = 0:dxa:L;

Da = 0.5; % diffusion coefficient
pa = dt*Da/(dxa^2);

Ka = 1.5; % branched coeff on overlap
Kb = 1.0; % bundled coeff on overlap
KaElse = 1.0; % branched coeff away from overlap
KbElse = 1.1; % bundled coeff away from overlap
alpha = 2;

F = @(U,V) -m0*U.*V;

rng(11)
a1 = 0.1 + 0.9.*rand(Na,1);
rng(10)
b1 = 0.1 + 0.9.*rand(Na,1);

bper=0.25;
blen=floor(bper*Na);
bound = (floor(Na/2)-floor(blen/2)):(floor(Na/2)+floor(blen/2));

vid = 0;
vidObj = VideoWriter('PDESolver1','MPEG-4');

method = 'crank-nicolsona';

switch method
    case{'crank-nicolson'}
        for i=1:Nt      
            % Laplacian difference operator with no flux boundary conditions
            % Crank-Nicolson operators
            II1 = speye(Na,Na);
            Lapdiff1 = spdiags(ones(Na,1)*[0 1 -2 1 0], [-Na+1, -1, 0 , 1, Na-1], Na, Na);
            Lapdiff1(1,1) = -2; Lapdiff1(1,2) = 1; Lapdiff1(1,Na) = 1;
            Lapdiff1(Na,1) = 1; Lapdiff1(Na,Na-1) = 1; Lapdiff1(Na,Na) = -2;
            Hm1 = II1+(pa/2)*Lapdiff1;
            Hs1 = II1-(pa/2)*Lapdiff1;
        
        
            diffRHSa1 = Hm1*a1;
            diffRHSb1 = Hm1*b1;
        
        
            rxna1 = dt*(F(a1,b1) + KaElse*(a1 - a1.*a1));
            rxnb1 = dt*(F(b1,a1) + KbElse*(b1 - b1.*b1));

            rxna1(bound) = dt*(F(a1(bound),b1(bound)) + Ka*(a1(bound) - a1(bound).*a1(bound)));
            rxnb1(bound) = dt*(F(b1(bound),a1(bound)) + Kb*(b1(bound) - b1(bound).*b1(bound)));


            a1 = Hs1\(diffRHSa1+rxna1);
            b1 = Hs1\(diffRHSb1+rxnb1);
        
            if i==Nt
                figure(1)
                plot(Xa,a1,'-o','markerfacecolor',[159 219 229]/255,'linewidth',3); hold on;
                plot(Xa,b1,'-ok','markerfacecolor','k','linewidth',3);
                set(gca,'fontname','times','fontsize',20); box on;
                set(gcf,'color','w');
                legend('branched','bundled');
                ylim([-1 1])
                hold off

                if vid==1
                    currframe = getframe(gcf);
                    writeVideo(vidObj,currframe);
                end
            end
        
        end
    % method of lines
    otherwise
        reltol=1.0e-04; abstol=1.0e-04;
        options=odeset('RelTol',reltol,'AbsTol',abstol);
        U0 = [a1,b1].';
        
        params = [Da,Ka,Kb,m0,Na,dxa];

        [teval,U] = ode15s(@(t,U) lvs_1D(t,U,params),[0,Tend],U0);
        U = U.';

        %keyboard();
        a1 = U(1:Na,end);
        b1 = U(Na+1:2*Na,end);
        figure(1)
        plot(Xa,a1,'-o','markerfacecolor',[159 219 229]/255,'linewidth',3); hold on;
        plot(Xa,b1,'-ok','markerfacecolor','k','linewidth',3);
end


% Determine if cell is polarized
a1New = a1;
a1New(a1New<0.5)=0;
if (a1New(1)~=0 && a1New(length(a1New))~=0)
    zeroInd1=find(a1New==0,1,'first');
    zeroInd2=find(a1New==0,1,'last');
    dirIndex1=ceil((zeroInd1+zeroInd2)/2) - 50;
else
    ind1=find(a1New~=0,1,'first');
    ind2=find(a1New~=0,1,'last');
    dirIndex1=ceil((ind1+ind2)/2);
end
b1New = b1;
b1New(b1New<0.5)=0;
if (b1New(1)~=0 && b1New(length(b1New))~=0)
    zeroInd1=find(b1New==0,1,'first');
    zeroInd2=find(b1New==0,1,'last');
    dirIndex2=ceil((zeroInd1+zeroInd2)/2) - 50;
else
    ind1=find(b1New~=0,1,'first');
    ind2=find(b1New~=0,1,'last');
    dirIndex2=ceil((ind1+ind2)/2);
end
if dirIndex1<1
    dirIndex1=dirIndex1+101;
end
if dirIndex2<1
    dirIndex2=dirIndex2+101;
end

if isempty(dirIndex1) || isempty(dirIndex2)
    sprintf('Not Polarized')
else
    sprintf('Polarized')
end


function dUdt = lvs_1D(~,U,params)
    Da = params(1); Ka = params(2); Kb = params(3); m0 = params(4);
    Na = params(5); dxa = params(6);
     
    dadt = zeros(Na,1);
    dbdt = zeros(Na,1);

    a(1:Na) = U(1:Na);
    b(1:Na) = U(Na+1:2*Na);

    for i=1:Na
        if i==1 
            dadt(1) = Ka*(a(1)-a(1)^2) - m0*a(1)*b(1) + Da*(a(2)-2.*a(1)+a(Na))./(dxa^2);
            dbdt(1) = Ka*(b(1)-b(1)^2) - m0*a(1)*b(1) + Da*(b(2)-2.*b(1)+b(Na))./(dxa^2);            
        elseif i==Na 
            dadt(Na) = Ka*(a(Na)-a(Na)^2) - m0*a(Na)*b(Na) + Da*(a(1)-2.*a(Na)+a(Na-1))./(dxa^2);
            dbdt(Na) = Ka*(b(Na)-b(Na)^2) - m0*a(Na)*b(Na) + Da*(b(1)-2.*b(Na)+b(Na-1))./(dxa^2);
        else 
            dadt(i) = Ka*(a(i)-a(i)^2) - m0*a(i)*b(i) + Da*(a(i+1)-2.*a(i)+a(i-1))./(dxa^2); 
            dbdt(i) = Ka*(b(i)-b(i)^2) - m0*a(i)*b(i) + Da*(b(i+1)-2.*b(i)+b(i-1))./(dxa^2);  
        end
        
    end

    dUdt = [dadt; dbdt];
end
