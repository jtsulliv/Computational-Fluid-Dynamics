%PROJECT 4 
%GAUSS-SEIDEL ALGORITHM
%WITHOUT SUCCESSIVE OVER-RELAXATION



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


RES = 1;                          %INITIALIZING RMS ERROR FOR RESIDUALS VECTOR
error = 0.0000001;                %ERROR TOLERANCE FOR RESIDUALS VECTOR
it_num = 0;
it_num_x = zeros(1, 100);
RES_Y = zeros(1, 100);


tic
while RES > error
    it_num = it_num + 1;
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
    
%RESIDUAL CHECK
R = zeros (length(2:(N_z-1)),length(2:(N_y-1)));
 
for k = 2:(N_y-1);
    for j = 2:(N_z-1);
        R(j,k) = (u(j,k) - uo(j,k))^2;
    end
end

RR = sqrt( sum( sum(R) ) );


R = zeros (length(2:(N_z-1)),length(2:(N_y-1)));

for k = 2:(N_y-1);
    for j = 2:(N_z-1);
        R(j,k) = (uo(j,k))^2;
    end
end

RRR = sqrt( sum( sum(R) ) );


RES = RR/RRR;
it_num_x (it_num) = it_num;
RES_Y (it_num) = RES;
uo = u;

end
toc

RES
it_num


yyy = (y_o/dy)+1;

u_25 = uo(:,yyy);

Z2=[z:-dz:0];



%COMPUTING THE ERROR OF THE NUMERICAL SCHEME WITH RESPECT TO THE ANALYTICAL
%SOLUTION BY COMPUTING THE NORM OF THE CONVERGED SOLUTION TO
%THE ANALYTICAL SOLUTION
U_t=transpose(u_25);

U_25 = zeros(1,(length(U_t)));

for ixy = 1:length(U_t);
    U_25(ixy) = U_t(length(U_t)-(ixy-1));
end


R_num = zeros(1,length(U_25));
R_den = zeros(1, length(U_25));
for i = 1:length(U_25);
    R_num(i) = (U_25(i) - U(i))^2;
    R_den(i) = (U(i))^2;
end


R_num_2 = sqrt ( sum (R_num) );
R_den_2 = sqrt ( sum (R_den) );

error_wrt_analytical = (R_num_2/R_den_2) * 100;

fprintf('The grid used is a %d X %d (y X z) equally spaced grid \n',N_y,N_z)
fprintf('The grid step sizes are %d (dy) and %d (dz) \n',dy,dz)
fprintf('The algorithm required %d iterations to achieve convergence \n',it_num)
fprintf('There is a %d percent error with respect to the analytical solution\n',error_wrt_analytical)

plot(Z,U,Z2,u_25,'O')
xlabel('z*')
ylabel('u*(y*=.25,z*)')
legend('Analytical Solution','Gauss-Seidel Solution')

figureHandle = gcf;
% %# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')


