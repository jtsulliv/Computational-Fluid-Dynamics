%PROJECT 3 - FINITE VOLUME 
%CONVECTION DIFFUSION EQUATION
%CENTRAL DIFFERENCE VERIFICATION - MESH REFINEMENT STUDY



clear;
clc;

%PARAMETERS
gamma = 0.1;                     %[kg/m*s]
u = 2.5;                         %[m/s]
rho = 1;                         %[kg/m^3]
L = 1;                           %[m]
phi_0 = 1;                       %BOUNDARY CONDITION AT x = 0
phi_L = 0;                       %BOUNDARY CONDITION AT x = L


%ANALYTICAL SOLUTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
N=1000
dx=L/N
x = dx/2;
PHI_A = zeros(N,1);
for i = 1:N;
    phi = ((exp((rho*u*x)/gamma)-1)/(exp((rho*u*L)/gamma)-1))*(phi_L-phi_0)+phi_0;
    PHI_A(i) = phi;
    x = x+dx;
end

X1 = (dx/2):dx:(L-(dx/2));
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%CENTRAL DIFFERENCING N=20%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POPULATING THE [A] MATRIX
tic
N = 5;                          %NUMBER OF CELLS
dx = L/N;                        %DISCRETIZATION IN X
D = gamma/dx;                    %DIFFUSION COEFFICIENT 
F = rho*u;                       %CONVECTION COEFFICIENT



a_E = D - (F/2);                 %EAST NODE COEFFICIENT
a=zeros(N+1,1)
for i=1:N
    a(i)=a_E
end
a(N+1,1)=1

a_W = D + (F/2);                 %WEST NODE COEFFICIENT
a=zeros(N+1,1)
for i=1:N
    a(i)=a_E
end
a(N+1,1)=1

a_P = a_E + a_W;                 %CENTRAL NODE COEFFICIENT

A = zeros(N+2);
A((N/N),(N/N)) = 1;
A((N/N),((N/N)+1)) = 1;
A(N+2,N+2) = 1;
A(N+2,N+1) = 1;

for i = 2:(N+1);
    A(i,i-1) = -a_W;
    A(i,i) = a_P;
    A(i,i+1) = -a_E;
end

%POPULATING THE BOUNDARY CONDITION MATRIX [B]
B = zeros(N+2,1);
B((N/N),1) = 2*phi_0;
B((N+2),1) = 2*phi_L;

%SOLVING FOR phi MATRIX
C = (A^-1)*B;
PHI_cen = zeros(N,1);
 for i=1:N;
     PHI_cen(i) = C(i+1);
 end

X2 = (dx/2):dx:(L-(dx/2));

toc 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









X = (dx/2):dx:(L-(dx/2));

fprintf('THE NUMBER OF NODES IS %d',N)

plot(X1,PHI_A,X2,PHI_cen)

xlabel('x-Location [m]')
ylabel('\phi')
title('Profile of property \phi versus x-Location')
legend('Analytical Solution','Central Differencing')



figureHandle = gcf;
% %# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')


    