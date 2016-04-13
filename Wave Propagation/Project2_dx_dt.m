%PROJECT 2 - FTCS SCHEME
%CONVECTION DIFFUSION EQUATION

clear;
clc;


%PARAMETERS
a = 2.5;                        %CONVECTION TERM CONSTANT [m/s]
alpha = 0.005;                  %DIFFUSION TERM CONSTANT [m^2/s]
t_f = 0.2;                      %FINAL TIME 
x_f = 1                         %FINAL LENGTH [m]
N = 21                          %NUMBER OF NODES IN x

c = 0.3
dx = (1*alpha)/a
dt = (c*dx)/a
d = (alpha*dt)/(dx^2)
Re = (a*dx)/alpha


 



    


