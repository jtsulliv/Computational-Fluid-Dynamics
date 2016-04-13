%PROJECT 3 - FINITE VOLUME 
%CONVECTION DIFFUSION EQUATION


tic
clear;
clc;

%PARAMETERS
gamma = 0.1;                     %[kg/m*s]
u = 2.5;                         %[m/s]
rho = 1;                         %[kg/m^3]
L = 1;                           %[m]
phi_0 = 1;                       %BOUNDARY CONDITION AT x = 0
phi_L = 0;                       %BOUNDARY CONDITION AT x = L
N = 20;                          %NUMBER OF CELLS
dx = L/N;                        %DISCRETIZATION IN X
D = gamma/dx;                    %DIFFUSION COEFFICIENT 
F = rho*u;                       %CONVECTION COEFFICIENT
Pe = F/D;                         %PECLET NUMBER


%ANALYTICAL SOLUTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = dx/2;
PHI_A = zeros(N,1);
for i = 1:N;
    phi = ((exp((rho*u*x)/gamma)-1)/(exp((rho*u*L)/gamma)-1))*(phi_L-phi_0)+phi_0;
    PHI_A(i) = phi;
    x = x+dx;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%CENTRAL DIFFERENCING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POPULATING THE [A] MATRIX
a_E = D - (F/2);                 %EAST NODE COEFFICIENT
a_E_cen = a_E;
a_W = D + (F/2);                 %WEST NODE COEFFICIENT
a_W_cen = a_W;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%UPWIND FOR CONVECTION, CENTRAL DIFFERENCE FOR DIFFUSION%%%%%%%%%%%%%%%%%%%
%POPULATING THE [A] MATRIX
a_E = D;                        %EAST NODE COEFFICIENT
a_E_upwind = a_E;
a_W = D + F;                    %WEST NODE COEFFICIENT
a_W_upwind = a_W;
a_P = a_E + a_W;                %CENTRAL NODE COEFFICIENT

AA = zeros(N+2);
AA((N/N),(N/N)) = 1;
AA((N/N),((N/N)+1)) = 1;
AA(N+2,N+2) = 1;
AA(N+2,N+1) = 1;

for i = 2:(N+1);
    AA(i,i-1) = -a_W;
    AA(i,i) = a_P;
    AA(i,i+1) = -a_E;
end


%POPULATING THE BOUNDARY CONDITION MATRIX [B]
BB = zeros(N+2,1);
BB((N/N),1) = 2*phi_0;
BB((N+2),1) = 2*phi_L;


%SOLVING FOR phi MATRIX
CC = (AA^-1)*BB;
PHI_upwind = zeros(N,1);
 for i=1:N;
     PHI_upwind(i) = CC(i+1);
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%HYBRID UPWIND AND CENTRAL DIFFERENCING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POPULATING THE [A] MATRIX
A_W = [F,(D+(F/2)),0];
A_E = [-F,(D-(F/2)),0];
a_W_hybrid = max(A_W);
a_E_hybrid = max(A_E);   


a_P_hybrid = a_E_hybrid + a_W_hybrid;             %CENTRAL NODE COEFFICIENT

AAA = zeros(N+2);
AAA((N/N),(N/N)) = 1;
AAA((N/N),((N/N)+1)) = 1;
AAA(N+2,N+2) = 1;
AAA(N+2,N+1) = 1;

for i = 2:(N+1);
    AAA(i,i-1) = -a_W_hybrid;
    AAA(i,i) = a_P_hybrid;
    AAA(i,i+1) = -a_E_hybrid;
end


%POPULATING THE BOUNDARY CONDITION MATRIX [B]
BBB = zeros(N+2,1);
BBB((N/N),1) = 2*phi_0;
BBB((N+2),1) = 2*phi_L;


%SOLVING FOR phi MATRIX
CCC = (AAA^-1)*BBB;
PHI_hybrid = zeros(N,1);
 for i=1:N;
     PHI_hybrid(i) = CCC(i+1);
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%POWER LAW%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_W = D*max([0,(1-0.1*abs(Pe))^5])+max([F,0]);
a_E = D*max([0,(1-0.1*abs(Pe))^5])+max([-F,0]);
a_P = a_W + a_E;

AAAA = zeros(N+2);
AAAA((N/N),(N/N)) = 1;
AAAA((N/N),((N/N)+1)) = 1;
AAAA(N+2,N+2) = 1;
AAAA(N+2,N+1) = 1;

for i = 2:(N+1);
    AAAA(i,i-1) = -a_W;
    AAAA(i,i) = a_P;
    AAAA(i,i+1) = -a_E;
end


%POPULATING THE BOUNDARY CONDITION MATRIX [B]
BBBB = zeros(N+2,1);
BBBB((N/N),1) = 2*phi_0;
BBBB((N+2),1) = 2*phi_L;

%SOLVING FOR phi MATRIX
CCCC = (AAAA^-1)*BBBB;
PHI_power = zeros(N,1);
 for i=1:N;
     PHI_power(i) = CCCC(i+1);
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%QUICK%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COEFFICIENTS FOR THE FIVE BAND MATRIX 
a_W = D + ((6/8)*F) + ((1/8)*F);
a_WW = -(1/8)*F;
a_E = D - ((3/8)*F); 
a_E_q = a_E;
a_EE = 0;
a_P = a_W + a_E + a_WW + a_EE;


%IMPLEMENTING THE PROPER COEFFICIENTS FOR THE GHOST NODES
AAAAA = zeros(N+4);
AAAAA(1,1) = 1;
AAAAA(1,2) = -2;
AAAAA(1,3) = 1;
AAAAA(2,1) = 0;
AAAAA(2,2) = 1;
AAAAA(2,3) = 1;
AAAAA((N+3),(N+2)) = 1;
AAAAA((N+3),(N+3)) = 1;
AAAAA((N+3),(N+4)) = 0;
AAAAA((N+4),(N+2)) = 1;
AAAAA((N+4),(N+3)) = -2;
AAAAA((N+4),(N+4)) = 1;





%POPULATING THE [A] MATRIX
for i = 3:(N+2);
    AAAAA(i,i-1) = -a_W;
    AAAAA(i,i) = a_P;
    AAAAA(i,i+1) = -a_E;
end



for i = 3:(N+2);
     AAAAA(i,i-2) = -a_WW;
end
 

%POPULATING THE BOUNDARY CONDITION MATRIX [B]
BBBBB = zeros(N+4,1);
BBBBB(2,1) = 2*phi_0;


%SOLVING FOR phi MATRIX
DDDDD = (AAAAA^-1)*(BBBBB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



QUICK = zeros(N,1);
for i=1:N;
    QUICK(i)=DDDDD(i+2);
end






Pe
X = (dx/2):dx:(L-(dx/2));

fprintf('THE NUMBER OF NODES IS %d',N)

plot(X,PHI_A,X,PHI_cen,X,PHI_upwind,X,PHI_hybrid,'o',X,PHI_power,'+',X,QUICK)

xlabel('x-Location [m]')
ylabel('\phi')
title('Profile of property \phi versus x-Location, N=20')
legend('Analytical Solution','Central Differencing','Upwind','Hybrid upwind and central differencing','Power Law','QUICK')



figureHandle = gcf;
% %# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')


toc
