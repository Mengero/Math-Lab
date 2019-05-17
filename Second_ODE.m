%% Find & Plot Exact Solution
%%
dydt = @(t,y)[y(2);y1*(1-y(1))*y(2)+y(1)];
f= @(t,y,Dy)y*(1-y)*Dy-y;
[Xe,Yexact]=ode45(@vdp1,[0 1],[1;1]);%%Use the ode45 funtion to solve the ODE initial value problem.
plot(Xe,Yexact(:,1),'-o',Xe,Yexact(:,2),'-o')
title('Solution of van der Pol Equation (\mu = 1) using ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')
Yexact = Yexact(:,1);
Exact_Sol=vpa([Xe Yexact],5)
%%
syms x y;
%% Step size
h=10^-4;
%% Range of x
xspan=[0,1];
steps = (xspan(2)-xspan(1))/h;
%% Initial values
x = zeros(steps,1);
y1 = zeros(steps,1);
y2 = zeros(steps,1);
y1(1)=1;x(1)=0;y2(1)=1;
%% Define Equation Set
%%
f= @(x,y1,y2)y1*(1-y1)*y2+y1;
g = @(x,y1,y2)y2;
%% Euler's Method
%%
for j=2:steps+1
    
    x(j,1)=x(j-1,1)+h;
    
    y2(j,1)=y2(j-1,1)+h*f(x(j-1,1),y1(j-1,1),y2(j-1,1));
    
    y1(j,1)=y1(j-1,1)+h*g(x(j-1,1,y1(j-1,1),y2(j-1,1)));
    
end
euler=vpa([x y],5)
%% Heun's Method
%%
for j=2:steps+1
    x(j,1)=x(j-1,1)+h;
    
    k1(j-1,1)=h*f(x(j-1,1),y1(j-1,1),y2(j-1,1));
    
    k2(j-1,1)=h*f(x(j-1,1),y1(j-1,1),y2(j-1,1));
    
    y2(j,1)=y2(j-1)+0.5*(k1(j-1)+k2(j-1));
    
    p1(j-1,1)=h*g(x(j-1,1),y1(j-1,1),y2(j-1,1));
    
    p2(j-1,1)=h*g(x(j-1,1),y1(j-1,1),y2(j-1,1));
    
    y1(j,1)=y1(j-1)+0.5*(p1(j-1)+p2(j-1));
end
heuns=vpa([x y1 y2],5)
