%PROJECT 1 - ANALYTICAL SOLUTION
%PART 6 - COMPARISON

clear;
clc;

%ANALYTICAL SOLUTION
%PARAMETERS
h = 0.04;                   %distance between two parallel plates extended to infinity [m]
nu = 0.000217;              %kinematic viscosity [m^2/s]
U_0 = 40;                   %velocity of lower plate [m/s]                     
k=1;
d=0.1356
t_final = 0.18
%NUMERICAL SOLUTION        %DISTANCE BETWEEN THE PLATES [m]
N=41;    

del_y = h/(N-1)             %DISCRETIZATION IN y [m]
del_t=(d*(del_y^2))/nu    %STABILITY CONDITION FOR THE FTCS SCHEME

%INITIAL CONDITIONS
u_o(1) = 40;                 
u_o(2:N) = 0;
j = 1;                       %TIME INDEX
t = 0;
while t<t_final;
    t=j*del_t;
    u(1)=40;                 %BOUNDARY CONDITION [m/s]
    u(N)=0;                  %BOUNDARY CONDITION [m/s]
    for i=2:(N-1);
        u(i)=u_o(i)+d*(u_o(i+1)-2*u_o(i)+u_o(i-1));
    end
    u_o=u;
    j=j+1;
end
Y_axis=0:del_y:h;

%TIME = .18 [s]
del_y = h/(N-1);                  %grid point increment
eta_1 = h/(2*sqrt(nu*t_final));
i=1;
    for y=0:del_y:h;
        Y(i) = y;
        eta  = y/(2*sqrt(nu*t_final)); 
            for n = 1:10;
                A(n) = erfc((2*(n-1)*eta_1)+eta);
                B(n) = erfc((2*(n)*eta_1)-eta);
            end
        X=sum(A);
        Z=sum(B);
        U(i)=U_0*(X-Z);
        i=i+1;
    end
    
U
Y
u_o
Y_axis

h = 0.04;                   %distance between two parallel plates extended to infinity [m]
nu = 0.000217;              %kinematic viscosity [m^2/s]
U_0_2 = 40;                   %velocity of lower plate [m/s]                     
k=1;
d=0.1356
t_final = 1.08

%NUMERICAL SOLUTION        %DISTANCE BETWEEN THE PLATES [m]
N=41;    

del_y = h/(N-1)             %DISCRETIZATION IN y [m]
del_t=(d*(del_y^2))/nu    %STABILITY CONDITION FOR THE FTCS SCHEME

%INITIAL CONDITIONS
u_o_2(1) = 40;                 
u_o_2(2:N) = 0;
j = 1;                       %TIME INDEX
t = 0;
while t<t_final;
    t=j*del_t;
    u_2(1)=40;                 %BOUNDARY CONDITION [m/s]
    u_2(N)=0;                  %BOUNDARY CONDITION [m/s]
    for i=2:(N-1);
        u_2(i)=u_o_2(i)+d*(u_o_2(i+1)-2*u_o_2(i)+u_o_2(i-1));
    end
    u_o_2=u_2;
    j=j+1;
end
Y_axis=0:del_y:h;

%TIME = .18 [s]
del_y = h/(N-1);                  %grid point increment
eta_1 = h/(2*sqrt(nu*t_final));
i=1;
    for y=0:del_y:h;
        Y(i) = y;
        eta  = y/(2*sqrt(nu*t_final)); 
            for n = 1:10;
                A(n) = erfc((2*(n-1)*eta_1)+eta);
                B(n) = erfc((2*(n)*eta_1)-eta);
            end
        X=sum(A);
        Z=sum(B);
        U_2(i)=U_0_2*(X-Z);
        i=i+1;
    end
    



plot(U,Y,u_o,Y_axis,'+',U_2,Y,u_o_2,Y_axis,'+')
title('Velocity profile of the FTCS Scheme and Analytical Solution')
ylabel('Distance from bottom plate [m]')
xlabel('Velocity [m/s]')
legend('Analytical Solution t = 0.18 s','Numerical Solution t = 0.18 s','Analytical Solution t = 1.08 s','Numerical Solution t = 1.08 s')



figureHandle = gcf;
%# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')

