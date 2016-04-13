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
t_final = 4
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

for i=1:length(Y_axis)
    u_a(i)=-1000*Y_axis(i)+40
end




plot(u_o,Y_axis,u_a,Y_axis,'+')
title('Velocity profile of the FTCS Scheme and Analytical Solution at steady state')
ylabel('Distance from bottom plate [m]')
xlabel('Velocity [m/s]')
legend('Numerical Solution','Analytical Solution')



figureHandle = gcf;
%# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')

