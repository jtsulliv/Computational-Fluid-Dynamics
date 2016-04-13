%FINAL PROJECT
%INCOMPRESSIBLE NAVIER-STOKES SOLVER

clear;
clc;


%PARAMETERS 
x = 1;
y = 1;
Re = 400;
Nx = 81;
Ny = Nx;
dx = x/(Nx-1);
dy = y/(Ny-1);
dt = 0.001
w = 1.73;
b = dx/dy;
n=0;
%STEP 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IN STEP ONE, INITIAL CONDITIONS AND BOUNDARY CONDITIONS ARE SET FOR u, v
%VORTICITY, AND STREAMFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SETTING INITIAL CONDITIONS FOR VELOCITY

uo = zeros(Ny,Nx);
i = 1;
for j = 1:Nx;
    uo(i,j) = 1;
end

u = zeros(Ny,Nx);
i = 1 ;
for j = 1:Nx;
    u(i,j) = 1;
end

vo = zeros(Ny,Nx);
v = zeros(Ny,Nx);


%VORTICITY AND STREAMFUNCTION B.C.'s and I.C.'s
zo=zeros(Ny,Nx);
z=zeros(Ny,Nx);
phi_o = zeros(Ny,Nx);
phi = zeros(Ny,Nx);
phi_oo = zeros(Ny,Nx);


RES_2 = 1;
error2 = 0.00001;
t=0;

while t < 30;

%STEP 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IN STEP 2, THE VORTICITY TRANSPORT IS SOLVED USING EXPLICIT FTCS
%!!!!!!!!!!!!!!NEED TO UPDATE B.C.'s INSIDE THE LOOP BECAUSE THEY ARE OF 
%NEUMANN TYPE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:(Ny-1);
    for j = 2:(Nx-1);
        z(i,j) = zo(i,j) + dt * ( -( ( uo(i+1,j)*zo(i+1,j) - uo(i-1,j)*zo(i-1,j) ) / (2*dx) ) - ( ( vo(i,j+1)*zo(i,j+1) - vo(i,j-1)*zo(i,j-1) ) / (2*dy) ) + (1/Re) * ( ( (zo(i+1,j) - 2*zo(i,j) + zo(i-1,j))/(dx^2) ) + ( (zo(i,j+1) - 2*zo(i,j) + zo(i,j-1))/(dy^2) ) ) ); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%STEP 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IN STEP 3, THE PSOR METHOD IS USED TO COMPUTE STREAMFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = 0.00001;
RES = 1;
it_num = 0;
while RES > error;
    
for i = 2:(Ny-1);
    for j = 2:(Nx-1);
        phi_oo(i,j) = (1-w) * phi_o(i,j) + (w/(2*(1+b^2))) * (dx^2 * z(i,j) + phi_o(i+1,j) + phi_oo(i-1,j) + b^2*(phi_o(i,j+1) + phi_oo(i,j-1)));
    end
end


%CONVERGENCE CHECK
%COMPUTING THE RESIDUAL VECTOR FOR COMPARISON TO THE ERROR TOLERANCE
R = zeros (length(2:(Ny-1)),length(2:(Nx-1)));
 
for i = 2:(Ny-1);
    for j = 2:(Nx-1);
        R(i,j) = (phi_oo(i,j) - phi_o(i,j))^2;
    end
end

RR = sqrt( sum( sum(R) ) );

for i = 2:(Ny-1);
    for j = 2:(Nx-1);
        R(i,j) = (phi_o(i,j))^2;
    end
end

RRR = sqrt( sum( sum(R) ) );

RES = RR/RRR;
phi_o = phi_oo;
it_num = it_num + 1;
end
it_num;
phi = phi_oo;

%STEP 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTE u and v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:(Ny-1);
    for j = 2:(Nx-1);
        u(i,j) = ( phi(i+1,j) - phi(i-1,j) ) / 2*dy;
    end
end



for i = 2:(Ny-1);
    for j = 2:(Nx-1);
        v(i,j) = -1*(( phi(i,j+1) - phi(i,j-1) ) / 2*dx);
    end
end

%COMPUTE B.C.'s for v



%VORTICITY BOUNDARY CONDITIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%VORTICITY BOUNDARY CONDITION FOR LEFT WALL
%SECOND ORDER CENTRAL DIFFERENCE
j = 1;
for i = 2:(Ny-1);
    z(i,j) = (-2/(dx^2))*phi(i,j+1);
end

%VORTICITY BOUNDARY CONDITION FOR RIGHT WALL
%SECOND ORDER CENTRAL DIFFERENCE
j = Nx;
for i = 2:(Ny-1);
    z(i,j) = (-2/(dx^2))*phi(i,j-1);
end

%VORTICITY BOUNDARY CONDITION FOR TOP WALL
%SECOND ORDER CENTRAL DIFFERENCE
i = 1;
for j = 2:(Nx-1);
    z(i,j) = (-phi(i+1,j)+dy)*(2/(dy^2));
end

%VORTICITY BOUNDARY CONDITION FOR BOTTOM WALL
%FIRST ORDER BACKWARD DIFFERENCE
i = Ny;
for j = 2:(Nx-1);
    z(i,j) = (-2/(dy^2))*phi(i-1,j);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



zo = z;
uo = u;
vo = v;
n=n+1;
t=n*dt
it_num = 0;
end


%PLOTTING GHIA's RESULTS
%Re = 400
Y_ghia = [0 0.0547 0.0625 0.0703 0.1016 0.1719 0.2813 0.4531 0.5000 0.6172 0.7344 0.8516 0.9531 0.9609 0.9688 0.9766 1.0000]
U_ghia = [0 -0.08186 -0.09266 -0.10338 -0.14612 -0.24299 -0.32726 -0.17119 -0.11477 0.02135 0.16256 0.29093 0.55892 0.61756 0.68439 0.75837 1.0000]
plot(U_ghia,Y_ghia,'o')

u(:,3)

