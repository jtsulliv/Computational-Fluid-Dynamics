%PROJECT 4 
%GAUSS-SEIDEL POINT SUCCESSIVE OVER-RELAXATION (PSOR)



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
bb = dy/dz

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



%GAUSS-SEIDEL POINT SUCCESSIVE OVER-RELAXATION (PSOR)

%DIRICHLET BOUNDARY CONDITIONS
u = zeros(N_z,N_y);

for j = 2:N_z;
     u(j,1) = 0;
end
 
for k = 1:N_y-1;
     u(N_z,k) = 0;
end


%INITIAL CONDITIONS
uo=zeros(N_z,N_y);


R = 1;                          %INITIALIZING RMS ERROR FOR RESIDUALS VECTOR
error = 0.00001;                %ERROR TOLERANCE FOR RESIDUALS VECTOR

tic
for it = 1:10000
    
    for k = 2:(N_y-1);
        for j = 2:(N_z-1);
            u(j,k) = (1/(2*(1+bb^2)))*(uo(j+1,k)+u(j-1,k)+(bb^2)*(uo(j,k+1)+u(j,k-1))+dy^2);
        end
    end


%NEUMANN BOUNDARY CONDITION y BOUNDARY (FLUX)
%BACKWARD DIFFERENCE
    k = N_y;
    for j = 2:(N_z-1);
        u(j,k) = (1/3) *  ( (4 * u(j,k-1)) - u(j,k-2) )  ;
    end

%NEUMANN BOUNDARY CONDITION z Boundary
%FORWARD DIFFERENCE
    j = 1;
    for k = 2:(N_y-1);
        u(j,k) = (1/3) *  ( (4 * u(j+1,k)) - u(j+2,k) )  ;
    end
uo = u;

end
toc




%CONVERGENCE TEST COMPUTING RESIDUAL VECTOR 
% iii=1;
% 
% j = (y_o/dy)+1;
% E = zeros(1,length(2:N_z));
% for k = 2:N_z;
%     E(iii) = (uo(j,k) - u(j,k))^2;
%     iii = iii+1;
% end
% 
% R = sqrt(sum(E));
% uo = u;
yyy = (y_o/dy)+1

u_25 = uo(:,yyy);

Z2=[z:-dz:0];

U_t=transpose(u_25);

diff = abs(U_t(1) - U(length(U)))

plot(Z,U,Z2,u_25)





