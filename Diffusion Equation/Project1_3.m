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



%FTCS NUMERICAL SCHEME 4 N=6
%INPUT PARAMETERS
nu = 0.000217;               %KINEMATIC VISCOSITY [m^2/s]
h = 0.04;                    %DISTANCE BETWEEN THE PLATES [m]
N_4 = 6;                      %NUMBER OF NODES IN THE y DIRECTION
del_t = 0.01                %DISCRETIZATION IN TIME [s]
t_final = 4                %END TIME [s]
U_0 = 40;

del_y_4 = h/(N_4-1)             %DISCRETIZATION IN y [m]
d_4 = (nu*del_t)/(del_y_4^2)    %STABILITY CONDITION FOR THE FTCS SCHEME

%INITIAL CONDITIONS
u_o_4(1) = 40;                 
u_o_4(2:N_4) = 0;
j = 1;                       %TIME INDEX
t = 0;
while t<t_final;
    t=j*del_t;
    u_4(1)=40;                 %BOUNDARY CONDITION [m/s]
    u_4(N_4)=0;                  %BOUNDARY CONDITION [m/s]
    for i=2:(N_4-1);
        u_4(i)=u_o_4(i)+d_4*(u_o_4(i+1)-2*u_o_4(i)+u_o_4(i-1));
    end
    u_o_4=u_4;
    j=j+1;
end
Y_axis_4=0:del_y_4:h;

%FTCS NUMERICAL SCHEME 5 N=21
%INPUT PARAMETERS
nu = 0.000217;               %KINEMATIC VISCOSITY [m^2/s]
h = 0.04;                    %DISTANCE BETWEEN THE PLATES [m]
N_5 = 21;                      %NUMBER OF NODES IN THE y DIRECTION
del_t = 0.01                %DISCRETIZATION IN TIME [s]
t_final = 4                %END TIME [s]
U_0 = 40;

del_y_5 = h/(N_5-1)             %DISCRETIZATION IN y [m]
d_5 = (nu*del_t)/(del_y_5^2)    %STABILITY CONDITION FOR THE FTCS SCHEME

%INITIAL CONDITIONS
u_o_5(1) = 40;                 
u_o_5(2:N_5) = 0;
j = 1;                       %TIME INDEX
t = 0;
while t<t_final;
    t=j*del_t;
    u_5(1)=40;                 %BOUNDARY CONDITION [m/s]
    u_5(N_5)=0;                  %BOUNDARY CONDITION [m/s]
    for i=2:(N_5-1);
        u_5(i)=u_o_5(i)+d_5*(u_o_5(i+1)-2*u_o_5(i)+u_o_5(i-1));
    end
    u_o_5=u_5;
    j=j+1;
end
Y_axis_5=0:del_y_5:h;


plot(u_o,Y_axis,u_o_4,Y_axis_4,'-.or',u_o_5,Y_axis_5,'+')
xlabel('Velocity (m/s)');
ylabel('Distance from the bottom plate (m)');
legend('Trial 1', 'Trial 4','Trial 5')

title('Velocity profile at t=4 s')

figureHandle = gcf;
%# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')

