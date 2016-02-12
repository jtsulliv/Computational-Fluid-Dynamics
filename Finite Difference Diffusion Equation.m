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
t=.1;
del_y = h/(N-1);                  %grid point increment
eta_1 = h/(2*sqrt(nu*t));
i=1;
    for y=0:del_y:h
        Y(i) = y;
        eta  = y/(2*sqrt(nu*t)) ;
            for n = 1:10
                A(n) = erfc((2*(n-1)*eta_1)+eta);
                B(n) = erfc((2*(n)*eta_1)-eta);
            end
        X=sum(A);
        Z=sum(B);
        U(i)=U_0*(X-Z);
        i=i+1;
    end
    
%TIME = .2 [s]
t=.2;
del_y = h/(N-1);                  %grid point increment
eta_1 = h/(2*sqrt(nu*t));
i=1;
    for y=0:del_y:h
        Y(i) = y;
        eta  = y/(2*sqrt(nu*t)); 
            for n = 1:10
                A(n) = erfc((2*(n-1)*eta_1)+eta);
                B(n) = erfc((2*(n)*eta_1)-eta);
            end
        X=sum(A);
        Z=sum(B);
        U_2(i)=U_0*(X-Z);
        i=i+1;
    end    

    
%TIME = .3 [s]
t=.3;
del_y = h/(N-1);                  %grid point increment
eta_1 = h/(2*sqrt(nu*t));
i=1;
    for y=0:del_y:h
        Y(i) = y;
        eta  = y/(2*sqrt(nu*t)) ;
            for n = 1:10
                A(n) = erfc((2*(n-1)*eta_1)+eta);
                B(n) = erfc((2*(n)*eta_1)-eta);
            end
        X=sum(A);
        Z=sum(B);
        U_3(i)=U_0*(X-Z);
        i=i+1;
    end        

    
%TIME = .4 [s]
t=.4;
del_y = h/(N-1);                  %grid point increment
eta_1 = h/(2*sqrt(nu*t));
i=1;
    for y=0:del_y:h
        Y(i) = y;
        eta  = y/(2*sqrt(nu*t)) ;
            for n = 1:10
                A(n) = erfc((2*(n-1)*eta_1)+eta);
                B(n) = erfc((2*(n)*eta_1)-eta);
            end
        X=sum(A);
        Z=sum(B);
        U_4(i)=U_0*(X-Z);
        i=i+1;
    end        

    
%TIME = .5 [s]
t=.5;
del_y = h/(N-1);                  %grid point increment
eta_1 = h/(2*sqrt(nu*t));
i=1;
    for y=0:del_y:h
        Y(i) = y;
        eta  = y/(2*sqrt(nu*t)) ;
            for n = 1:10
                A(n) = erfc((2*(n-1)*eta_1)+eta);
                B(n) = erfc((2*(n)*eta_1)-eta);
            end
        X=sum(A);
        Z=sum(B);
        U_5(i)=U_0*(X-Z);
        i=i+1;
    end    
    

%TIME = 1 [s]
t=1;
del_y = h/(N-1);                  %grid point increment
eta_1 = h/(2*sqrt(nu*t));
i=1;
    for y=0:del_y:h
        Y(i) = y;
        eta  = y/(2*sqrt(nu*t)) ;
            for n = 1:10
                A(n) = erfc((2*(n-1)*eta_1)+eta);
                B(n) = erfc((2*(n)*eta_1)-eta);
            end
        X=sum(A);
        Z=sum(B);
        U_6(i)=U_0*(X-Z);
        i=i+1;
    end      
    
    
%TIME = 2 [s]
t=2;
del_y = h/(N-1);                  %grid point increment
eta_1 = h/(2*sqrt(nu*t));
i=1;
    for y=0:del_y:h
        Y(i) = y;
        eta  = y/(2*sqrt(nu*t)) ;
            for n = 1:10
                A(n) = erfc((2*(n-1)*eta_1)+eta);
                B(n) = erfc((2*(n)*eta_1)-eta);
            end
        X=sum(A);
        Z=sum(B);
        U_7(i)=U_0*(X-Z);
        i=i+1;
    end
    
    
%TIME = 3 [s]
t=3;
del_y = h/(N-1);                  %grid point increment
eta_1 = h/(2*sqrt(nu*t));
i=1;
    for y=0:del_y:h
        Y(i) = y;
        eta  = y/(2*sqrt(nu*t)) ;
            for n = 1:10
                A(n) = erfc((2*(n-1)*eta_1)+eta);
                B(n) = erfc((2*(n)*eta_1)-eta);
            end
        X=sum(A);
        Z=sum(B);
        U_8(i)=U_0*(X-Z);
        i=i+1;
    end       

    
%TIME = 4 [s]
t=4;
del_y = h/(N-1);                  %grid point increment
eta_1 = h/(2*sqrt(nu*t));
i=1;
    for y=0:del_y:h
        Y(i) = y;
        eta  = y/(2*sqrt(nu*t)) ;
            for n = 1:10
                A(n) = erfc((2*(n-1)*eta_1)+eta);
                B(n) = erfc((2*(n)*eta_1)-eta);
            end
        X=sum(A);
        Z=sum(B);
        U_9(i)=U_0*(X-Z);
        i=i+1;
    end
    
    
plot(U,Y,U_2,Y,U_3,Y,U_4,Y,U_5,Y,U_6,Y,U_7,Y,U_8,Y,U_9,Y)
title('Velocity profile at different times')
xlabel('Velocity (m/s)')
ylabel('Distance from the bottom plate (m)')
legend('t=0.1 s','t=0.2 s','t=0.3 s','t=0.4 s','t=0.5 s','t=1 s','t=2 s','t=3 s','t=4 s')

figureHandle = gcf;
%# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')
