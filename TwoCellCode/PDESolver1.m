
dt = 0.01;
Nt = 1/dt; % number of time steps
Na = 101; % number of space steps
dxa = 5.0/((Na-1)/2);

Da = 0.5; % diffusion coefficient
pa = dt*Da/(dxa^2);

Ka = 1; % branched coeff
Kb = 1; % bundled coeff
alpha = 2;

F = @(U,V) -m0*U.*V;

a1 = 0.1 + 0.9.*rand(101,1);
b1 = 0.1 + 0.9.*rand(101,1);


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


    rxna1 = dt*(F(a1,b1) + K*(a1 - a1.*a1));
    rxnb1 = dt*(F(b1,a1) + K*(b1 - b1.*b1));


    a1 = Hs1\(diffRHSa1+rxna1);
    b1 = Hs1\(diffRHSb1+rxnb1);

end

figure(1)
plot(a1)
hold on;
plot(b1)
hold off;
legend('branched','bundled')