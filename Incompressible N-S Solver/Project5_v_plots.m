%FINAL PROJECT
%INCOMPRESSIBLE NAVIER-STOKES SOLVER

clear;
clc;


%PARAMETERS
%PARAMETERS 
x = 1;
y = 1;
Re = 400;
Nx = 21;
Ny = Nx;
dx = x/(Nx-1);
dy = y/(Ny-1);
dt = 0.002;                  %FROM THE STABILITY CONDITIONS FOR FTCS--ZIKANOV
w = 1.6;
b = dx/dy;
n=0;

%INITIAL CONDITIONS AND BOUNDARY CONDITIONS
u = zeros(Nx,Ny);               
j = Ny;
for i=1:Nx;
    u(i,j) = 1;
end

uo = zeros(Nx,Ny);
j = Ny;
for i=1:Nx;
    uo(i,j) = 1;
end

v = zeros(Nx,Ny);
vo = zeros(Ny,Nx);
z = zeros(Nx,Ny);               %VORTICITY
zo = zeros(Nx,Ny);              %VORTICITY
psi = zeros(Nx,Ny);             %STREAMFUNCTION
psi_o = zeros(Nx,Ny);           %STREAMFUNCTION



t=0;
error2 = 0.00001;
RES2 = 1;
tic

while RES2 > error2
%EXPLICIT SCHEME TO COMPUTE VORTICITY AT t^n+1

for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        z(i,j) = zo(i,j) + dt*...
            ( -( (uo(i+1,j)*zo(i+1,j) - uo(i-1,j)*zo(i-1,j))/(2*dx) )...
            -( (vo(i,j+1)*zo(i,j+1) - vo(i,j-1)*zo(i,j-1))/(2*dy) )...
            + (1/Re)* ( ((zo(i+1,j) - 2*zo(i,j) + zo(i-1,j))/(dx^2))...
            + ((zo(i,j+1) - 2*zo(i,j) + zo(i,j-1))/(dy^2)) ) ) ;
    end
end





%COMPUTE STREAMFUNCTION AT t^n+1
%PSOR ITERATIVE METHOD

error = 0.00001;
RES = 1;
it_num = 0;

while RES > error

for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        psi(i,j) = (1-w)*psi_o(i,j)+(w/(2*(1+b^2))) *...
            (dx^2*zo(i,j) + psi_o(i+1,j) + psi(i-1,j)...
            +b^2 * (psi_o(i,j+1) + psi(i,j-1)) );
    end
end


%PSOR CONVERGENCE CHECK
%COMPUTING THE RESIDUAL VECTOR FOR COMPARISON TO THE ERROR TOLERANCE
R = zeros (length(2:(Nx-1)),length(2:(Ny-1)));
 
for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        R(i,j) = (psi(i,j) - psi_o(i,j))^2;
    end
end

RR = sqrt( sum( sum(R) ) );

for i = 2:(Ny-1);
    for j = 2:(Nx-1);
        R(i,j) = (psi_o(i,j))^2;
    end
end

RRR = sqrt( sum( sum(R) ) );

RES = RR/RRR;
psi_o = psi;
it_num = it_num + 1;

end

it_num;



%COMPUTE VELOCITY
for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        u(i,j) = (psi(i,j+1) - psi(i,j-1))/(2*dy);
    end
end

for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        v(i,j) = -((psi(i+1,j) - psi(i-1,j))/(2*dx));
    end
end


%COMPUTE VORTICITY B.C.'s
%LEFT WALL (1,j)
for j = 2:(Ny-1);
    z(1,j) = -2*psi(2,j)/dx^2;
end

%RIGHT WALL (Nx,j)
for j = 2:(Ny-1);
    z(Nx,j) = -2*psi(Nx-1,j)/dx^2;
end

%TOP WALL (i,Ny)
for i = 2:(Nx-1);
    z(i,Ny) = -2*psi(i,Ny-1)/(dy^2) - 2/dy;
end

%BOTTOM WALL (i,1)
for i = 2:(Nx-1);
    z(i,1) = -2*psi(i,2)/dy^2;
end





%VORTICITY CONVERGENCE CHECK
%COMPUTING THE RESIDUAL VECTOR FOR COMPARISON TO THE ERROR TOLERANCE

if t>=1
    R = zeros (length(2:(Nx-1)),length(2:(Ny-1)));
 
for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        R(i,j) = (z(i,j) - zo(i,j))^2;
    end
end

RR = sqrt( sum( sum(R) ) );

for i = 2:(Ny-1);
    for j = 2:(Nx-1);
        R(i,j) = (zo(i,j))^2;
    end
end

RRR = sqrt( sum( sum(R) ) );

RES2 = RR/RRR;
else
    RES2=1;
end



vo = v;
uo = u;
zo = z;
n=n+1;
t=n*dt;


end
toc



V_21 = v(:,(((Ny - 1)/2)+1))
X_21 = [0:dx:x]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nx=Ny=21
%PARAMETERS
%PARAMETERS 
x = 1;
y = 1;
Re = 400;
Nx = 41;
Ny = Nx;
dx = x/(Nx-1);
dy = y/(Ny-1);
dt = 0.002;                  %FROM THE STABILITY CONDITIONS FOR FTCS--ZIKANOV
w = 1.6;
b = dx/dy;
n=0;

%INITIAL CONDITIONS AND BOUNDARY CONDITIONS
u = zeros(Nx,Ny);               
j = Ny;
for i=1:Nx;
    u(i,j) = 1;
end

uo = zeros(Nx,Ny);
j = Ny;
for i=1:Nx;
    uo(i,j) = 1;
end

v = zeros(Nx,Ny);
vo = zeros(Ny,Nx);
z = zeros(Nx,Ny);               %VORTICITY
zo = zeros(Nx,Ny);              %VORTICITY
psi = zeros(Nx,Ny);             %STREAMFUNCTION
psi_o = zeros(Nx,Ny);           %STREAMFUNCTION



t=0;
error2 = 0.00001;
RES2 = 1;
tic

while RES2 > error2
%EXPLICIT SCHEME TO COMPUTE VORTICITY AT t^n+1

for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        z(i,j) = zo(i,j) + dt*...
            ( -( (uo(i+1,j)*zo(i+1,j) - uo(i-1,j)*zo(i-1,j))/(2*dx) )...
            -( (vo(i,j+1)*zo(i,j+1) - vo(i,j-1)*zo(i,j-1))/(2*dy) )...
            + (1/Re)* ( ((zo(i+1,j) - 2*zo(i,j) + zo(i-1,j))/(dx^2))...
            + ((zo(i,j+1) - 2*zo(i,j) + zo(i,j-1))/(dy^2)) ) ) ;
    end
end





%COMPUTE STREAMFUNCTION AT t^n+1
%PSOR ITERATIVE METHOD

error = 0.00001;
RES = 1;
it_num = 0;

while RES > error

for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        psi(i,j) = (1-w)*psi_o(i,j)+(w/(2*(1+b^2))) *...
            (dx^2*zo(i,j) + psi_o(i+1,j) + psi(i-1,j)...
            +b^2 * (psi_o(i,j+1) + psi(i,j-1)) );
    end
end


%PSOR CONVERGENCE CHECK
%COMPUTING THE RESIDUAL VECTOR FOR COMPARISON TO THE ERROR TOLERANCE
R = zeros (length(2:(Nx-1)),length(2:(Ny-1)));
 
for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        R(i,j) = (psi(i,j) - psi_o(i,j))^2;
    end
end

RR = sqrt( sum( sum(R) ) );

for i = 2:(Ny-1);
    for j = 2:(Nx-1);
        R(i,j) = (psi_o(i,j))^2;
    end
end

RRR = sqrt( sum( sum(R) ) );

RES = RR/RRR;
psi_o = psi;
it_num = it_num + 1;

end

it_num;



%COMPUTE VELOCITY
for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        u(i,j) = (psi(i,j+1) - psi(i,j-1))/(2*dy);
    end
end

for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        v(i,j) = -((psi(i+1,j) - psi(i-1,j))/(2*dx));
    end
end


%COMPUTE VORTICITY B.C.'s
%LEFT WALL (1,j)
for j = 2:(Ny-1);
    z(1,j) = -2*psi(2,j)/dx^2;
end

%RIGHT WALL (Nx,j)
for j = 2:(Ny-1);
    z(Nx,j) = -2*psi(Nx-1,j)/dx^2;
end

%TOP WALL (i,Ny)
for i = 2:(Nx-1);
    z(i,Ny) = -2*psi(i,Ny-1)/(dy^2) - 2/dy;
end

%BOTTOM WALL (i,1)
for i = 2:(Nx-1);
    z(i,1) = -2*psi(i,2)/dy^2;
end





%VORTICITY CONVERGENCE CHECK
%COMPUTING THE RESIDUAL VECTOR FOR COMPARISON TO THE ERROR TOLERANCE

if t>=1
    R = zeros (length(2:(Nx-1)),length(2:(Ny-1)));
 
for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        R(i,j) = (z(i,j) - zo(i,j))^2;
    end
end

RR = sqrt( sum( sum(R) ) );

for i = 2:(Ny-1);
    for j = 2:(Nx-1);
        R(i,j) = (zo(i,j))^2;
    end
end

RRR = sqrt( sum( sum(R) ) );

RES2 = RR/RRR;
else
    RES2=1;
end



vo = v;
uo = u;
zo = z;
n=n+1;
t=n*dt;


end
toc


V_41 = v(:,(((Ny - 1)/2)+1))
X_41 = [0:dx:x]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nx=Ny=41
%PARAMETERS
%PARAMETERS 


x = 1;
y = 1;
Re = 400;
Nx = 81;
Ny = Nx;
dx = x/(Nx-1);
dy = y/(Ny-1);
dt = 0.002;                  %FROM THE STABILITY CONDITIONS FOR FTCS--ZIKANOV
w = 1.6;
b = dx/dy;
n=0;

%INITIAL CONDITIONS AND BOUNDARY CONDITIONS
u = zeros(Nx,Ny);               
j = Ny;
for i=1:Nx;
    u(i,j) = 1;
end

uo = zeros(Nx,Ny);
j = Ny;
for i=1:Nx;
    uo(i,j) = 1;
end

v = zeros(Nx,Ny);
vo = zeros(Ny,Nx);
z = zeros(Nx,Ny);               %VORTICITY
zo = zeros(Nx,Ny);              %VORTICITY
psi = zeros(Nx,Ny);             %STREAMFUNCTION
psi_o = zeros(Nx,Ny);           %STREAMFUNCTION



t=0;
error2 = 0.00001;
RES2 = 1;
tic

while RES2 > error2
%EXPLICIT SCHEME TO COMPUTE VORTICITY AT t^n+1

for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        z(i,j) = zo(i,j) + dt*...
            ( -( (uo(i+1,j)*zo(i+1,j) - uo(i-1,j)*zo(i-1,j))/(2*dx) )...
            -( (vo(i,j+1)*zo(i,j+1) - vo(i,j-1)*zo(i,j-1))/(2*dy) )...
            + (1/Re)* ( ((zo(i+1,j) - 2*zo(i,j) + zo(i-1,j))/(dx^2))...
            + ((zo(i,j+1) - 2*zo(i,j) + zo(i,j-1))/(dy^2)) ) ) ;
    end
end





%COMPUTE STREAMFUNCTION AT t^n+1
%PSOR ITERATIVE METHOD

error = 0.00001;
RES = 1;
it_num = 0;

while RES > error

for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        psi(i,j) = (1-w)*psi_o(i,j)+(w/(2*(1+b^2))) *...
            (dx^2*zo(i,j) + psi_o(i+1,j) + psi(i-1,j)...
            +b^2 * (psi_o(i,j+1) + psi(i,j-1)) );
    end
end


%PSOR CONVERGENCE CHECK
%COMPUTING THE RESIDUAL VECTOR FOR COMPARISON TO THE ERROR TOLERANCE
R = zeros (length(2:(Nx-1)),length(2:(Ny-1)));
 
for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        R(i,j) = (psi(i,j) - psi_o(i,j))^2;
    end
end

RR = sqrt( sum( sum(R) ) );

for i = 2:(Ny-1);
    for j = 2:(Nx-1);
        R(i,j) = (psi_o(i,j))^2;
    end
end

RRR = sqrt( sum( sum(R) ) );

RES = RR/RRR;
psi_o = psi;
it_num = it_num + 1;

end

it_num;



%COMPUTE VELOCITY
for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        u(i,j) = (psi(i,j+1) - psi(i,j-1))/(2*dy);
    end
end

for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        v(i,j) = -((psi(i+1,j) - psi(i-1,j))/(2*dx));
    end
end


%COMPUTE VORTICITY B.C.'s
%LEFT WALL (1,j)
for j = 2:(Ny-1);
    z(1,j) = -2*psi(2,j)/dx^2;
end

%RIGHT WALL (Nx,j)
for j = 2:(Ny-1);
    z(Nx,j) = -2*psi(Nx-1,j)/dx^2;
end

%TOP WALL (i,Ny)
for i = 2:(Nx-1);
    z(i,Ny) = -2*psi(i,Ny-1)/(dy^2) - 2/dy;
end

%BOTTOM WALL (i,1)
for i = 2:(Nx-1);
    z(i,1) = -2*psi(i,2)/dy^2;
end





%VORTICITY CONVERGENCE CHECK
%COMPUTING THE RESIDUAL VECTOR FOR COMPARISON TO THE ERROR TOLERANCE

if t>=1
    R = zeros (length(2:(Nx-1)),length(2:(Ny-1)));
 
for i = 2:(Nx-1);
    for j = 2:(Ny-1);
        R(i,j) = (z(i,j) - zo(i,j))^2;
    end
end

RR = sqrt( sum( sum(R) ) );

for i = 2:(Ny-1);
    for j = 2:(Nx-1);
        R(i,j) = (zo(i,j))^2;
    end
end

RRR = sqrt( sum( sum(R) ) );

RES2 = RR/RRR;
else
    RES2=1;
end



vo = v;
uo = u;
zo = z;
n=n+1;
t=n*dt;


end
toc


V_81 = v(:,(((Ny - 1)/2)+1))
X_81 = [0:dx:x]


X_ghia = [0.0000 0.0625 0.0703 0.0781 0.0938 0.1563 0.2266 0.2344 0.5000 0.8047 0.8594 0.9063 0.9453 0.9531 0.9609 0.9688 1.0000]
V_ghia = [0.00000 0.18360 .19713 0.20920 0.22965 0.28124 0.30203 0.30174 0.05186 -0.38598 -0.44993 -.23827 -0.22847 -0.19254 -0.15663 -0.12146 0.00000]

% plot(X_21,V_21,X_ghia,V_ghia,'o')


plot(X_21,V_21,X_41,V_41,X_81,V_81,X_ghia,V_ghia,'o')
ylabel('v velocity at geometric center')
xlabel('x location')
title('v velocity profile at geometric center')
legend('21x21','41x41','81x81','Ghia')

figureHandle = gcf;
% %# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')


