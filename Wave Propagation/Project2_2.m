%PROJECT 2 - FTCS SCHEME
%CONVECTION DIFFUSION EQUATION
%PROGRAM LOOPS THROUGH IN TIME FROM 0.05s to 0.20s


clear;
clc;


%PARAMETERS
for t_f = .05:0.05:0.2;                       %FINAL TIME 
a = 2.5;                        %CONVECTION TERM CONSTANT [m/s]
alpha = 0.005;                  %DIFFUSION TERM CONSTANT [m^2/s]
x_f = 1;                        %FINAL LENGTH [m]
xo = 0.2;

%STABILITY CONDITION FROM CELL Re
c = 0.2;                         %COURANT NUMBER
Re = .5;                          %CELL Re
dx = (Re*alpha)/a;               %DISCRETIZATION IN x [m]
dt = (c*dx)/a;                   %DISCRETIZATION IN t [s]
d = (alpha*dt)/(dx^2);           %FOURIER NUMBER     
N = (x_f/dx);                    %NUMBER OF NODES IN x
T = floor(t_f/dt);             %NUMBER OF NODES IN t              

fprintf('THE COMPUTATIONAL GRID IS %d x %d (N x t)',N,T)

%%%%%%%%%%%%%%%%%%%%FTCS NUMERICAL SCHEME%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIAL CONDITIONS
for i=1:N
    x = (i-1)*dx;
    X(i)=x;
    if x < 0.2;
        u_o(i) = 1;
    else
        if x == 0.2;
            u_o(i) = 0.5;
        else
            if x > 0.2;
                u_o(i) = 0.0;
            end
        end
    end
end;
j = 1;                          %TIME INDEX
t = 0;


%FTCS SCHEME
while t < t_f;
    t = j*dt;                  
    u(1) = 1;                   %BOUNDARY CONDITION [m/s]
    u(N) = 0;                   %BOUNDARY CONDITION [m/s]
    for i = 2:(N-1);
        u(i) = u_o(i)-(c/2)*(u_o(i+1)-u_o(i-1))+d*(u_o(i+1)-2*(u_o(i))+u_o(i-1));
    end
    u_o = u;
    j = j+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%ANALYTICAL SOLUTION CALCULATION%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt=t_f;

for i=1:N
    x = (i-1)*dx;
    u_a = 1-(1/2)*(1+erf((x-xo-a*tt)/(2*sqrt(alpha*tt))));
    U_a(i)=u_a;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%FIRST ORDER UPWIND FOR CONVECTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FTCS FOR DIFFUSION SCHEME%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INITIAL CONDITIONS
for i=1:N
    x = (i-1)*dx;
    X(i)=x;
    if x < 0.2;
        u_2(i) = 1;
    else
        if x == 0.2;
            u_2(i) = 0.5;
        else
            if x > 0.2;
                u_2(i) = 0.0;
            end
        end
    end
end;
j = 1;                          %TIME INDEX
t = 0;

while t < t_f;
    t = j*dt;                  
    u(1) = 1;                   %BOUNDARY CONDITION [m/s]
    u(N) = 0;                   %BOUNDARY CONDITION [m/s]
    for i = 2:(N-1);
        u(i) = u_2(i)-c*(u_2(i)-u_2(i-1))+d*(u_2(i+1)-2*(u_2(i))+u_2(i-1));
    end
    u_2 = u;
    j = j+1;
end;



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
plot(X,U_a,X,u_o,X,u_2,X,u_oo)
title({'Velocity Profile';'c=0.2, d=0.4'})
xlabel({'Distance from boundary';'[m]'})
ylabel({'Velocity';'[m/s]'})
legend('Analytical Solution','FTCS Convection Diffusion','Upwind Convection FTCS Diffusion','MacCormack')
hold on

t_f
dx
dt
c
d

end




figureHandle = gcf;
% %# make all text in the figure to size 14 and bold
 set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')




    


