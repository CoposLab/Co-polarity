
D=1;

a1 = 0.1 + 0.9.*rand(101,1);

m = 0;

% Function that codes the pde
function [c,f,s] = branched(x,t,u,dudx)
c = 1;
f = D*dudx;
s = u-u.*u;
end

% Initial condition
function u0 = branchedic(x)
u0 = a1;
end

% Boundary condition
function [pl,ql,pr,qr] = branchedbc(xl,ul,xr,ur,t)
pl = 0; %ignored by solver since m=1
ql = 0; %ignored by solver since m=1
pr = ur-besselj(0,n)*exp(-n^2*t);
qr = 0; 
end