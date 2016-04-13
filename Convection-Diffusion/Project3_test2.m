%QUICK Scheme

gamma = 0.1                     %[kg/m*s]
u = 2.5                         %[m/s]
rho = 1                         %[kg/m^3]
L = 1                           %[m]
phi_0 = 1                       %BOUNDARY CONDITION AT x = 0
phi_L = 0                       %BOUNDARY CONDITION AT x = L
N = 20                           %NUMBER OF CELLS
dx = L/N                        %DISCRETIZATION IN X
D = gamma/dx;                   %DIFFUSION COEFFICIENT 
F = rho*u;  


%ANALYTICAL SOLUTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = dx/2;
PHI_A = zeros(N,1);
for i = 1:N;
    phi = ((exp((rho*u*x)/gamma)-1)/(exp((rho*u*L)/gamma)-1))*(phi_L-phi_0)+phi_0;
    PHI_A(i) = phi;
    x = x+dx;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%QUICK%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COEFFICIENTS FOR THE FIVE BAND MATRIX 
a_W = D + ((6/8)*F) + ((1/8)*F)
a_WW = -(1/8)*F
a_E = D - ((3/8)*F) 
a_E_q = a_E
a_EE = 0
a_P = a_W + a_E + a_WW + a_EE


%IMPLEMENTING THE PROPER COEFFICIENTS FOR THE GHOST NODES
AAAAA = zeros(N+4);
AAAAA(1,1) = 1
AAAAA(1,2) = -2
AAAAA(1,3) = 1
AAAAA(2,1) = 0
AAAAA(2,2) = 1
AAAAA(2,3) = 1
AAAAA((N+3),(N+2)) = 1
AAAAA((N+3),(N+3)) = 1
AAAAA((N+3),(N+4)) = 0
AAAAA((N+4),(N+2)) = 1
AAAAA((N+4),(N+3)) = -2
AAAAA((N+4),(N+4)) = 1





%POPULATING THE [A] MATRIX
for i = 3:(N+2)
    AAAAA(i,i-1) = -a_W
    AAAAA(i,i) = a_P
    AAAAA(i,i+1) = -a_E
end

AAAAA

for i = 3:(N+2)
     AAAAA(i,i-2) = -a_WW
end
 

%POPULATING THE BOUNDARY CONDITION MATRIX [B]
BBBBB = zeros(N+4,1);
BBBBB(2,1) = 2*phi_0


%SOLVING FOR phi MATRIX
CCCCC = linsolve(AAAAA,BBBBB);
DDDDD = (AAAAA^-1)*(BBBBB)
% PHI_QUICK = zeros(N,1);
%  for i=1:N;
%      PHI_QUICK(i) = CCCCC(i+1);
%  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



QUICK = zeros(N,1)
for i=1:N
    QUICK(i)=DDDDD(i+2)
end

X = (dx/2):dx:(L-(dx/2));



plot(X,PHI_A,X,QUICK)

xlabel('x-Location [m]')
ylabel('\phi')
title('Profile of property \phi versus x-Location')
legend('Analytical Solution','QUICK')

figureHandle = gcf;
% %# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')






