%PROJECT 1 - FINITE DEIFFERENCE - DIFFUSION EQUATION
clear;
clc;

%FTCS NUMERICAL SCHEME 1 del_t = 0.01 s
%INPUT PARAMETERS
nu = 0.000217;               %KINEMATIC VISCOSITY [m^2/s]
h = 0.04;                    %DISTANCE BETWEEN THE PLATES [m]
N = 11;                      %NUMBER OF NODES IN THE y DIRECTION
del_t = 0.01                %DISCRETIZATION IN TIME [s]
t_final = 4                %END TIME [s]
U_0 = 40;

del_y = h/(N-1)             %DISCRETIZATION IN y [m]
d = (nu*del_t)/(del_y^2)    %STABILITY CONDITION FOR THE FTCS SCHEME

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

%FTCS NUMERICAL SCHEME 2 - del_t = 0.037 s 
%INPUT PARAMETERS
nu = 0.000217;               %KINEMATIC VISCOSITY [m^2/s]
h = 0.04;                    %DISTANCE BETWEEN THE PLATES [m]
N = 11;                      %NUMBER OF NODES IN THE y DIRECTION
del_t_2 = 0.037                %DISCRETIZATION IN TIME [s]
t_final = 4                %END TIME [s]
U_0 = 40;

del_y = h/(N-1)             %DISCRETIZATION IN y [m]
d_2 = (nu*del_t_2)/(del_y^2)    %STABILITY CONDITION FOR THE FTCS SCHEME

%INITIAL CONDITIONS
u_o_2(1) = 40;                 
u_o_2(2:N) = 0;
j = 1;                       %TIME INDEX
t = 0;
while t<t_final;
    t=j*del_t_2;
    u_2(1)=40;                 %BOUNDARY CONDITION [m/s]
    u_2(N)=0;                  %BOUNDARY CONDITION [m/s]
    for i=2:(N-1);
        u_2(i)=u_o_2(i)+d_2*(u_o_2(i+1)-2*u_o_2(i)+u_o_2(i-1));
    end
    u_o_2=u_2;
    j=j+1;
end
Y_axis=0:del_y:h;

%FTCS NUMERICAL SCHEME 3 - del_t = 0.04 s
%INPUT PARAMETERS
nu = 0.000217;               %KINEMATIC VISCOSITY [m^2/s]
h = 0.04;                    %DISTANCE BETWEEN THE PLATES [m]
N = 11;                      %NUMBER OF NODES IN THE y DIRECTION
del_t_3 = 0.04                %DISCRETIZATION IN TIME [s]
t_final = 4                %END TIME [s]
U_0 = 40;

del_y = h/(N-1)             %DISCRETIZATION IN y [m]
d_3 = (nu*del_t_3)/(del_y^2)    %STABILITY CONDITION FOR THE FTCS SCHEME

%INITIAL CONDITIONS
u_o_3(1) = 40;                 
u_o_3(2:N) = 0;
j = 1;                       %TIME INDEX
t = 0;
while t<t_final;
    t=j*del_t_3;
    u_3(1)=40;                 %BOUNDARY CONDITION [m/s]
    u_3(N)=0;                  %BOUNDARY CONDITION [m/s]
    for i=2:(N-1);
        u_3(i)=u_o_3(i)+d_3*(u_o_3(i+1)-2*u_o_3(i)+u_o_3(i-1));
    end
    u_o_3=u_3;
    j=j+1;
end
Y_axis=0:del_y:h;





plot(u_o,Y_axis,u_o_2,Y_axis,u_o_3,Y_axis)
xlabel('Velocity (m/s)');
ylabel('Distance from the bottom plate (m)');
legend('Trial 1','Trial 2','Trial 3')
title('Velocity profile at t=4 s')

figureHandle = gcf;
%# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')
        

