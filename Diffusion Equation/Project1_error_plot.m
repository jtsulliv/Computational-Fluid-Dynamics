%PROJECT 1 - ANALYTICAL SOLUTION
%PART 5 ERROR PLOT

clear;
clc;

%ANALYTICAL SOLUTION
%PARAMETERS
h = 0.04;                   %distance between two parallel plates extended to infinity [m]
nu = 0.000217;              %kinematic viscosity [m^2/s]
U_0 = 40;                   %velocity of lower plate [m/s]                     
k=1;
d=0.001
t_final = 0.18
%NUMERICAL SOLUTION        %DISTANCE BETWEEN THE PLATES [m]
for N=11:51;    

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

for i=1:N
    num(i)=(u_o(i)-U(i))^2
    den(i)=(U(i))^2
end
A=sqrt(sum(num))
B=sqrt(sum(den))
error=(A/B)*100


DEL_Y(k)=del_y
ERR(k)=error
k=k+1
end

plot(DEL_Y,ERR)
title('Error vs. \Delta y')
xlabel('\Delta y [m]')
ylabel('Error')



figureHandle = gcf;
%# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')

