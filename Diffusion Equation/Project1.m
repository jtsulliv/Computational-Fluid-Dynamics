% Project 1

clear;
clc;

%Part 1)

h = 0.40;                   %distance between two parallel plates extended to infinity [m]
nu = 0.000217;              %kinematic viscosity [m^2/s]
U_0 = 40;                   %velocity of lower plate [m/s]
tol = 0.00001;              %convergence criteria
num=10;                     %number of grid points
inc=h/num;                  %grid point increment

%y = 0.05 m
%this section computes the analytical solution of u(0.05,t) until steady-state
y=0.04;


t=0.1;
eta = y/(2*sqrt(nu*t));
eta_1 = h/(2*sqrt(nu*t));
n = 0;
a = erfc((2*n*eta_1)+eta);
b = erfc((2*(n+1)*eta_1)-eta);
u_old = U_0*(a-b);

for n = 1:10
    a = a+erfc((2*n*eta_1)+eta);
    b = b+erfc((2*(n+1)*eta_1)-eta);
    u_old = U_0*(a-b);
end
t
u_old


t=0.2;
eta = y/(2*sqrt(nu*t));
eta_1 = h/(2*sqrt(nu*t));
n = 0;
a = erfc((2*n*eta_1)+eta);
b = erfc((2*(n+1)*eta_1)-eta);
u_new = U_0*(a-b);

for n = 1:10
    a = a+erfc((2*n*eta_1)+eta);
    b = b+erfc((2*(n+1)*eta_1)-eta);
    u_new = U_0*(a-b);
end
t
u_new


while abs(u_old-u_new)>tol
    t=t+0.1;
    eta = y/(2*sqrt(nu*t));
    eta_1 = h/(2*sqrt(nu*t));
    n=0;
    a = erfc((2*n*eta_1)+eta);
    b = erfc((2*(n+1)*eta_1)-eta);
    u_old=u_new;
    for n = 1:10
        a = a+erfc((2*n*eta_1)+eta);
        b = b+erfc((2*(n+1)*eta_1)-eta);
        u_new = U_0*(a-b);
    end
end
u_new
t



       
        
        

