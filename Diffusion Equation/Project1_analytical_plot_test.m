%PROJECT 1 - ANALYTICAL SOLUTION


clear;
clc;

%Part 1)
%PARAMETERS
h = 0.04;                   %distance between two parallel plates extended to infinity [m]
nu = 0.000217;              %kinematic viscosity [m^2/s]
U_0 = 40;                   %velocity of lower plate [m/s]
N=11;                     %number of grid points


%TIME = .1 [s]
for t=0.1:.1:4;
del_y = h/(N-1);                  %grid point increment
eta_1 = h/(2*sqrt(nu*t));
i=1
    for y=0:del_y:h
        Y(i) = y
        eta  = y/(2*sqrt(nu*t)) 
            for n = 1:10
                A(n) = erfc((2*(n-1)*eta_1)+eta)
                B(n) = erfc((2*(n)*eta_1)-eta)
            end
        X=sum(A)
        Z=sum(B)
        U(i)=U_0*(X-Z)
        i=i+1
    end     
    
plot(U,Y)
hold on
end

 
title('Velocity profile at different times')
xlabel('Velocity (m/s)')
ylabel('Distance from the bottom plate (m)')
        

