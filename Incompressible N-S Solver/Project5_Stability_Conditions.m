%PROJECT 5
%STABILITY CONDITIONS


%PARAMETERS 
x = 1;
y = 1;
Re = 400;
Nx = 11;
Ny = Nx;
dx = x/(Nx-1);
dy = y/(Ny-1);



dt1 = (1/(2*(1/Re))) * ( (1/(dx^2)) + (1/(dy^2)) )^-1



dt2 = (2*(1/Re))/(1^2 + 1^2)


x=[1 2 3]
y=[1 2 3]
plot(x,y)
title({'You can do it','with a cell array'})
