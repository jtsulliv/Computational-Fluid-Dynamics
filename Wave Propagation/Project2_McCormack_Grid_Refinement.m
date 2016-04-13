%PROJECT 2 - FTCS SCHEME
%CONVECTION DIFFUSION EQUATION
%McCORMACK GRID REFINEMENT STUDY


clear;
clc;


%PARAMETERS
t_f = .05;                       %FINAL TIME 
a = 2.5;                        %CONVECTION TERM CONSTANT [m/s]
alpha = 0.005;                  %DIFFUSION TERM CONSTANT [m^2/s]
x_f = 1;                        %FINAL LENGTH [m]
xo = 0.2;

%TRIAL 1
dx = 0.01;                       %DISCRETIZATION IN x [m]
dt = 0.001;                   %DISCRETIZATION IN t [s]
c = (a*dt)/dx;                          %COURANT NUMBER
d = (alpha*dt)/(dx^2);           %FOURIER NUMBER     
N = (x_f/dx);                    %NUMBER OF NODES IN x
T = floor(t_f/dt);             %NUMBER OF NODES IN t              

fprintf('THE COMPUTATIONAL GRID IS %d x %d (N x t)',N,T)
dx
dt

%MACCORMACK METHOD FOR CONVECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SECOND ORDER CENTRAL SPACE DIFFUSION SCHEME%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INITIAL CONDITIONS
for i=1:N
    x = (i-1)*dx;
    X(i)=x;
    if x < 0.2;
        u_oo(i) = 1;
    else
        if x == 0.2;
            u_oo(i) = 0.5;
        else
            if x > 0.2;
                u_oo(i) = 0.0;
            end
        end
    end
end;
j = 1;                          %TIME INDEX
t = 0;
X1=X;


while t < t_f;
    t = j*dt;                  
    u_m(1) = 1;                   %BOUNDARY CONDITION [m/s]
    u_m(N) = 0;                   %BOUNDARY CONDITION [m/s]
    u_c(1)=1;
    u_c(N)=0;
    %PREDICTOR
    for i = 2:(N-1);
        u_c(i) = u_oo(i)-c*(u_oo(i+1)-u_oo(i))+d*(u_oo(i+1)-2*u_oo(i)+u_oo(i-1));
    end
    %CORRECTOR
    for i = 2:(N-1);
        u_m(i) = .5* (u_oo(i)+u_c(i) - c*(u_c(i)-u_c(i-1)) + d*(u_c(i+1)-2*u_c(i)+u_c(i-1)));
    end
    u_oo = u_m;
    j = j+1;
end;
u_oo1=u_oo;



%TRIAL 2
dx = 0.005;                       %DISCRETIZATION IN x [m]
dt = 0.0005;                   %DISCRETIZATION IN t [s]
c = (a*dt)/dx                          %COURANT NUMBER
d = (alpha*dt)/(dx^2);           %FOURIER NUMBER     
N = (x_f/dx);                    %NUMBER OF NODES IN x
T = floor(t_f/dt);             %NUMBER OF NODES IN t              

fprintf('THE COMPUTATIONAL GRID IS %d x %d (N x t)',N,T)
dx
dt


%MACCORMACK METHOD FOR CONVECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SECOND ORDER CENTRAL SPACE DIFFUSION SCHEME%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INITIAL CONDITIONS
for i=1:N
    x = (i-1)*dx;
    X(i)=x;
    if x < 0.2;
        u_oo(i) = 1;
    else
        if x == 0.2;
            u_oo(i) = 0.5;
        else
            if x > 0.2;
                u_oo(i) = 0.0;
            end
        end
    end
end;
j = 1;                          %TIME INDEX
t = 0;
X2=X;


while t < t_f;
    t = j*dt;                  
    u_m(1) = 1;                   %BOUNDARY CONDITION [m/s]
    u_m(N) = 0;                   %BOUNDARY CONDITION [m/s]
    u_c(1)=1;
    u_c(N)=0;
    %PREDICTOR
    for i = 2:(N-1);
        u_c(i) = u_oo(i)-c*(u_oo(i+1)-u_oo(i))+d*(u_oo(i+1)-2*u_oo(i)+u_oo(i-1));
    end
    %CORRECTOR
    for i = 2:(N-1);
        u_m(i) = .5* (u_oo(i)+u_c(i) - c*(u_c(i)-u_c(i-1)) + d*(u_c(i+1)-2*u_c(i)+u_c(i-1)));
    end
    u_oo = u_m;
    j = j+1;
end;
u_oo2=u_oo;


%TRIAL 3
dx = 0.0025;                       %DISCRETIZATION IN x [m]
dt = 0.00025;                   %DISCRETIZATION IN t [s]
c = (a*dt)/dx                          %COURANT NUMBER
d = (alpha*dt)/(dx^2);           %FOURIER NUMBER     
N = (x_f/dx);                    %NUMBER OF NODES IN x
T = floor(t_f/dt);             %NUMBER OF NODES IN t              

fprintf('THE COMPUTATIONAL GRID IS %d x %d (N x t)',N,T)
dx
dt


%MACCORMACK METHOD FOR CONVECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SECOND ORDER CENTRAL SPACE DIFFUSION SCHEME%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INITIAL CONDITIONS
for i=1:N
    x = (i-1)*dx;
    X(i)=x;
    if x < 0.2;
        u_oo(i) = 1;
    else
        if x == 0.2;
            u_oo(i) = 0.5;
        else
            if x > 0.2;
                u_oo(i) = 0.0;
            end
        end
    end
end;
j = 1;                          %TIME INDEX
t = 0;
X3=X;


while t < t_f;
    t = j*dt;                  
    u_m(1) = 1;                   %BOUNDARY CONDITION [m/s]
    u_m(N) = 0;                   %BOUNDARY CONDITION [m/s]
    u_c(1)=1;
    u_c(N)=0;
    %PREDICTOR
    for i = 2:(N-1);
        u_c(i) = u_oo(i)-c*(u_oo(i+1)-u_oo(i))+d*(u_oo(i+1)-2*u_oo(i)+u_oo(i-1));
    end
    %CORRECTOR
    for i = 2:(N-1);
        u_m(i) = .5* (u_oo(i)+u_c(i) - c*(u_c(i)-u_c(i-1)) + d*(u_c(i+1)-2*u_c(i)+u_c(i-1)));
    end
    u_oo = u_m;
    j = j+1;
end;
u_oo3=u_oo;







%%%%%%%%%%%%%%%%%ANALYTICAL SOLUTION CALCULATION%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt=t_f;

for i=1:N
    x = (i-1)*dx;
    u_a = 1-(1/2)*(1+erf((x-xo-a*tt)/(2*sqrt(alpha*tt))));
    U_a(i)=u_a;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





plot(X,U_a,X1,u_oo1,X2,u_oo2,':',X3,u_oo3,'.')
title({'MacCormack Grid Refinement';'Velocity Profile'})
xlabel({'Distance from boundary';'[m]'})
ylabel({'Velocity';'[m/s]'})
legend('Analytical Solution','Grid 1','Grid 2','Grid 3')



figureHandle = gcf;
% %# make all text in the figure to size 14 and bold
 set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')




    


