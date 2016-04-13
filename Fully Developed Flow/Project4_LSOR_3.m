%PROJECT 4 
%GAUSS-SEIDEL ALGORITHM
%WITH SUCCESSIVE OVER-RELAXATION
%LSOR


clear;
clc;

mu = 0.4;                   %DYNAMIC VISCOSITY [N-s/m^2]
L = 1.5;                    %DUCT LENGTH [m]
h = 1.0;                    %DUCT HEIGHT
y = 1;
z = h/L;
N_y = 81;
N_z = 54;
dy = y/(N_y-1);
dz = z/(N_z-1);
y_o = .25;                  %NON-DIMENSIONAL y TERM
bb = dy/dz;
w = 1

%ANALYTICAL SOLUTION
ii = 1;
U = zeros(1,length(0:dz:z));
Z = zeros(1,length(0:dz:z));

for z_o=0:dz:z;

%INITIALIZING THE ITERATIVE PROCESS FOR THE ANALYTICAL SOLUTION
%n = 1:2:3
n = 3;
i = 1;
for m = 1:2:n;
        a = (cos((m*pi*(1-y_o))/2))/(m^3);
        b = cosh((m*pi*((h/L)-z_o))/2);
        c = cosh((m*pi*h)/(2*L));
        d = (-1)^((m-1)/2);
        r = (d*(1-(b/c))*a);
        A(i) = r;
        i=i+1;
 end
 R1 = sum(A);

 
%n = 1:2:5
n = 5;
i = 1;
for m = 1:2:n;
        a = (cos((m*pi*(1-y_o))/2))/(m^3);
        b = cosh((m*pi*((h/L)-z_o))/2);
        c = cosh((m*pi*h)/(2*L));
        d = (-1)^((m-1)/2);
        r = (d*(1-(b/c))*a);
        A(i) = r;
        i=i+1;
 end
R2 = sum(A);


%ITERATING TO FIND THE ANALYTICAL SOLUTION
i = 1;
tol = 0.0001;
n = 7;

    while abs(R1-R2) > tol;    
        R1 = R2;
        A = zeros(1,((n/2)+.5));
            for m = 1:2:n;
                a = (cos((m*pi*(1-y_o))/2))/(m^3);
                b = cosh((m*pi*((h/L)-z_o))/2);
                c = cosh((m*pi*h)/(2*L));
                d = (-1)^((m-1)/2);
                r = (d*(1-(b/c))*a);
                A(i) = r;
                i=i+1;
            end
        R = sum(A);
        R2 = R;
        n = n + 2;
    end

u_a = (16/(pi^3))*R2;
U(ii) = u_a;
Z(ii) = z_o;
ii=ii+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GAUSS-SEIDEL LINE SUCCESSIVE OVER-RELAXATION (PSOR)

%DIRICHLET BOUNDARY CONDITIONS
u = zeros(N_z,N_y);

for k = 1:(N_z-1);
    u(k,1) = 0;
end

for j = 2:(N_y-1);
    u(1,j) = 0;
end



%INITIAL CONDITIONS
uo=zeros(N_z,N_y);

nn = (N_y-1);

%BUILDING THE VECTORS FOR TDMAsolver
tic
for it=1:1000
for k = 2:(N_z-1)
        j = 2;
        a(1,1) = 0;
        b(1,1) = -2*(1 + bb^2);
        c(1,1) = w;
        d(1,1) = -(2*(1-w)*(1+bb^2)*uo(k,j)) - w*bb^2*(uo(k+1,j)+u(k-1,j)) - dy^2;
        for i = 2:(nn-1);
            a(i,1) = w;
            b(i,1) = -2*(1 + bb^2);
            c(i,1) = w;
            j=j+1;
            d(i,1) = -(2*(1-w)*(1+bb^2)*uo(k,j)) - w*bb^2*(uo(k+1,j)+u(k-1,j)) - dy^2;
        end
        a(nn,1) = -1;
        b(nn,1) = 1;
        c(nn,1) = 0;
        d(nn,1) = 0;
        
        u(k,(2:N_y)) = TDMAsolver(a,b,c,d);
end

k = N_z;
a = zeros((N_y-2),1);
b = ones((N_y-2),1);
c = zeros((N_y-2),1);
d = zeros((N_y-2),1);
for tt = 1:(N_y-2);
    tt;
    d(tt,1) = u((k-1),(tt+1));
end

u(k,(2:(N_y-1))) = TDMAsolver(a,b,c,d);

uo = u;
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTTING THE RESULTS 

ind = (1+(y_o/dy));

Z;
U;
U_LSOR = transpose(uo(:,ind));

plot(Z,U,Z,U_LSOR,'o')











