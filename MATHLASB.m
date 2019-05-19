%% Adams-Bashworth-Moulton Predictor, corrector and other methods
%% initialization
%%
syms x y;
% Step size
h=10^-4;
% Range of x
xspan=[0,0.5];
steps = 0.5/h;
% Initial values
x = zeros(steps,1);     %%Set up x and y matrices.
y = zeros(steps,1);
% Define Initial Condition.
y(1)=1;x(1)=x0;
%% Define f(x,y) 
%%
f= @(x,y)x.^2+y.^3;
%% Exact Solution
%%
[Xe,Yexact]=ode45(f,xspan,y0);%%Use the ode45 funtion to solve the ODE initial value problem.
Xe=Xe(1:4:end);%%Xe is the number of t. Take one after 4 steps.
Yexact=Yexact(1:4:end);%%Take one after 4 steps until the end.
Exact_Sol=vpa([Xe Yexact],5)%%Change all the numbers into decimals with 5 significant numbers.
%% Euler's Method
%%
for j=2:steps+1%%   When using the Adam method, we need four initial value. We can use Euler or Heun to calculate. Notice that "steps+1"
               %%   should be replaced by 4.
    
    x(j,1)=x(j-1,1)+h;
    
    y(j,1)=y(j-1,1)+h*f(x(j-1,1),y(j-1,1));
    
end
euler=vpa([x y],5)
%% Heun's Method
%%
for j=2:steps+1
    x(j,1)=x(j-1,1)+h;
    
    k1(j-1,1)=h*f( x(j-1), y(j-1) );
    
    k2(j-1,1)=h*f( x(j-1)+h, y(j-1)+k1(j-1));
    
    y(j,1)=y(j-1)+0.5*(k1(j-1)+k2(j-1));
end
heuns=vpa([x y],5)
%% RK 4th order method
%%
for j=2:steps+1
    x(j,1)=x(j-1)+h;
    
    k1(j-1,1)=h*f( x(j-1), y(j-1) );
    
    k2(j-1,1)=h*f( x(j-1)+h/2, y(j-1)+0.5*k1(j-1) );
    
    k3(j-1,1)=h*f( x(j-1)+h/2, y(j-1)+0.5*k2(j-1) );
    
    k4(j-1,1)=h*f( x(j-1)+h, y(j-1)+k3(j-1) );
    
    y(j,1)=y(j-1)+(1/6)*(k1(j-1)+2*k2(j-1)+2*k3(j-1)+k4(j-1));
end
rk4=vpa([x y],5)
%% Adam-Bashworth predictor
%%
for k=5:steps+1
x(k,1)=x(k-1)+h;
   
y(k,1)=y(k-1) +(h/24)*( -9*f(x(k-4),y(k-4)) +37*f(x(k-3),y(k-3))...
                        -59*f(x(k-2),y(k-2)) +55*f(x(k-1),y(k-1)));
end
p=y;
Adam_Bashworth=vpa([x p],5)
%% Adam-Moulton corrector
%%
for k=5:steps+1
x(k,1)=x(k-1)+h;
   
y(k,1)=y(k-1) +(h/24)*( f(x(k-3),y(k-3)) -5*f(x(k-2),y(k-2))...
                        +19*f(x(k-1),y(k-1)) +9*f(x(k),p(k)));
end
Adam_Moulton=vpa([x y],5)
%%
