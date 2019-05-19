%% ReadMe: To solve 2nd-ODE, we need to transfer the equation into euqation set consisting of two 1st-ODE.
%%         That is, suppose y1 = y, y2 = y' (This is the first equation). Then substitute them into the original equation to get another
%%         equation. 

%%  For further information: https://ww2.mathworks.cn/help/matlab/math/solve-nonstiff-odes.html

%% Find & Plot Exact Solution
%%
dydt = @(t,y)[y(2);y1*(1-y(1))*y(2)-y(1)];
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
f= @(x,y1,y2)y1*(1-y1)*y2-y1;
g = @(x,y1,y2)y2;
%% Euler's Method
%%
for j=2:steps+1
    
    x(j,1)=x(j-1,1)+h;
    
    y1(j,1)=y1(j-1,1)+h*f(x(j-1,1),y1(j-1,1),y2(j-1,1));
    
    y2(j,1)=y2(j-1,1)+h*g(x(j-1,1),y1(j-1,1),y2(j-1,1));
    
end
euler=vpa([x y1 y2],5)
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
%% Adam-Bashworth predictor
%%
for k=5:steps+1
x(k,1)=x(k-1)+h;
   
y2(k,1)=y2(k-1) +(h/24)*( -9*f(x(k-4),y1(k-4),y2(k-4)) +37*f(x(k-3),y1(k-3),y2(k-3))...
                        -59*f(x(k-2),y1(k-2),y2(k-2)) +55*f(x(k-1),y1(k-1),y2(k-1)));
                    
y1(k,1)=y1(k-1) +(h/24)*( -9*g(x(k-4),y1(k-4),y2(k-4)) +37*g(x(k-3),y1(k-3),y2(k-3))...
                        -59*g(x(k-2),y1(k-2),y2(k-2)) +55*g(x(k-1),y1(k-1),y2(k-1)));
end
p1=y1;
p2=y2;
Adam_Bashworth=vpa([x p1 p2],5)
%% Adam-Moulton corrector
%%
for k=5:steps+1
x(k,1)=x(k-1)+h;
   
y2(k,1)=y2(k-1) +(h/24)*( f(x(k-3),y1(k-3),y2(k-3)) -5*f(x(k-2),y1(k-2),y2(k-2))...
                        +19*f(x(k-1),y1(k-1),y2(k-1)) +9*f(x(k),p1(k),p2(k)));
y1(k,1)=y1(k-1) +(h/24)*( g(x(k-3),y1(k-3),y2(k-3)) -5*g(x(k-2),y1(k-2),y2(k-2))...
                        +19*g(x(k-1),y1(k-1),y2(k-1)) +9*g(x(k),p1(k),p2(k)));

end
Adam_Moulton=vpa([x y1 y2],5)
plot(x,y1(:,1),'-o',x,y2(:,1),'-o')
title('Solution of van der Pol Equation (\mu = 1) using ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')
