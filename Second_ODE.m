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