clc
clear
close all

format long g

xtar = 20;
r1 = 200;
r2 = 50;
rair = 1.225;
cd = 0.47 + 2/80;
g = 9.81;

t0 = 0;
tf = 10;
dt = 0.001;

x0 = 0;
y0 = 0;

psif = 0;
Ff = 0;

u0 = 8; %b1
v0 = 20; %b2
E = 0.1; %b3
M = 1.3;
epsilon = 10^(-5);
kill = false;
count = 1;

for p = 1:20 %ALM
    
    for k = 1:500
R = (0.75/pi)^(1/3) * (0.5/r1 + E/r2)^(1/3);
if imag(R) ~= 0
    kill = true;
    break;
end
A = pi * R^2;
m = 0.5+E;

F1 = @(u,v) -(2*m)^(-1)*rair*cd*A*sqrt(u^2 + v^2)*u;
F2 = @(u,v) -g - (2*m)^(-1)*rair*cd*A*sqrt(u^2 + v^2)*v;

% [x1,y1] = RK4(F1,F2,x0,y0,u0,v0,0.1,tf);
% [x2,y2] = RK4(F1,F2,x0,y0,u0,v0,0.01,tf);
% [x3,y3] = RK4(F1,F2,x0,y0,u0,v0,0.001,tf);
% [x4,y4] = RK4(F1,F2,x0,y0,u0,v0,0.0001,tf);
% 
% figure;
% plot(x1,y1)
% hold on
% plot(x2,y2)
% plot(x3,y3)
% plot(x4,y4)
% hold off
% grid on
% xlabel('x [m]')
% ylabel('y [m]')
% legend('Δt=0.1','Δt=0.01','Δt=0.001','Δt=0.0001')

[x,y,u,v,T,k1u,k1v,k2u,k2v,k3u,k3v,k4u,k4v] = RK4(F1,F2,x0,y0,u0,v0,dt,tf);
Tp(count)=T;
count=count+1;


[Df1,Df2] = FD(F1,F2,x0,y0,u0,v0,dt,tf,epsilon,M,xtar);

E = E + epsilon;
R = (0.75/pi)^(1/3) * (0.5/r1 + E/r2)^(1/3);
A = pi * R^2;
m = 0.5+E;
F1 = @(u,v) -(2*m)^(-1)*rair*cd*A*sqrt(u^2 + v^2)*u;
F2 = @(u,v) -g - (2*m)^(-1)*rair*cd*A*sqrt(u^2 + v^2)*v;

[x5,y5,~,~,T5] = RK4(F1,F2,x0,y0,u0,v0,dt,tf);

E = E - 2*epsilon;
R = (0.75/pi)^(1/3) * (0.5/r1 + E/r2)^(1/3);
A = pi * R^2;
m = 0.5+E;
F1 = @(u,v) -(2*m)^(-1)*rair*cd*A*sqrt(u^2 + v^2)*u;
F2 = @(u,v) -g - (2*m)^(-1)*rair*cd*A*sqrt(u^2 + v^2)*v;

[x6,y6,~,~,T6] = RK4(F1,F2,x0,y0,u0,v0,dt,tf);

Df3 = ((T5 - T6) + M*((x5(end)-xtar)^2 - (x6(end) -xtar)^2) + (y5(end) - y6(end)))/(2 * epsilon);

E = E + epsilon;

if E - 0.01 * Df3 < 0
    kill = true;
    break;
end
E = E - 0.01  * Df3;
u0 = u0 - 0.01 * Df1;
v0 = v0 - 0.01 * Df2;
    end
    if kill 
        break;
    end
    if M <30
        M = M*1.03;
    end
end


function [Df1,Df2] = FD(F1,F2,x0,y0,u0,v0,step,tf,epsilon,M,xtar)
    Df1 = 0;
    Df2 = 0;
    [x1,y1,~,~,T1] = RK4(F1,F2,x0,y0,u0+epsilon,v0,step,tf);
    [x2,y2,~,~,T2] = RK4(F1,F2,x0,y0,u0-epsilon,v0,step,tf);
    Df1 = ((T1-T2) + M*((x1(end)-xtar)^2 - (x2(end) -xtar)^2) + (y1(end) - y2(end)))/(2*epsilon);
    [x3,y3,~,~,T3] = RK4(F1,F2,x0,y0,u0,v0+epsilon,step,tf);
    [x4,y4,~,~,T4] = RK4(F1,F2,x0,y0,u0,v0-epsilon,step,tf);
    Df2 = ((T3-T4) + + M*((x3(end)-xtar)^2 - (x4(end) -xtar)^2) + (y3(end) - y4(end)))/(2*epsilon);
end


function [x,y,u,v,T,k1u,k1v,k2u,k2v,k3u,k3v,k4u,k4v] = RK4(F1,F2,x0,y0,u0,v0,step,tf)
   
x = zeros(tf/step,1);
x(1) = x0;
y = zeros(tf/step,1);
y(1) = y0;
u = zeros(tf/step,1);
u(1) = u0; %b1
v = zeros(tf/step,1);
v(1) = v0; %b2
k1x = zeros(tf/step,1);
k2x = zeros(tf/step,1);
k3x = zeros(tf/step,1);
k4x = zeros(tf/step,1);

k1y = zeros(tf/step,1);
k2y = zeros(tf/step,1);
k3y = zeros(tf/step,1);
k4y = zeros(tf/step,1);

k1u = zeros(tf/step,1);
k2u = zeros(tf/step,1);
k3u = zeros(tf/step,1);
k4u = zeros(tf/step,1);

k1v = zeros(tf/step,1);
k2v = zeros(tf/step,1);
k3v = zeros(tf/step,1);
k4v = zeros(tf/step,1);

for i = 1:(tf/step - 1)
    
    k1x(i) = step * u(i);
    k1y(i) = step * v(i);
    k1u(i) = step * F1(u(i),v(i));
    k1v(i) = step * F2(u(i),v(i));

    k2x(i) = step * (u(i) + k1u(i)/2);
    k2y(i) = step * (v(i) + k1v(i)/2);
    k2u(i) = step * F1(u(i) + k1u(i)/2,v(i) + k1v(i)/2);
    k2v(i) = step * F2(u(i) + k1u(i)/2,v(i) + k1v(i)/2);

    k3x(i) = step * (u(i) + k2u(i)/2);
    k3y(i) = step * (v(i) + k2v(i)/2);
    k3u(i) = step * F1(u(i) + k2u(i)/2,v(i) + k2v(i)/2);
    k3v(i) = step * F2(u(i) + k2u(i)/2,v(i) + k2v(i)/2);

    k4x(i) = step * (u(i) + k3u(i));
    k4y(i) = step * (v(i) + k3v(i));
    k4u(i) = step * F1(u(i) + k3u(i)/2,v(i) + k3v(i)/2);
    k4v(i) = step * F2(u(i) + k3u(i)/2,v(i) + k3v(i)/2);

    x(i+1) = x(i) + (k1x(i) + 2*k2x(i) + 2*k3x(i) + k4x(i))/6;
    y(i+1) = y(i) + (k1y(i) + 2*k2y(i) + 2*k3y(i) + k4y(i))/6;
    u(i+1) = u(i) + (k1u(i) + 2*k2u(i) + 2*k3u(i) + k4u(i))/6;
    v(i+1) = v(i) + (k1v(i) + 2*k2v(i) + 2*k3v(i) + k4v(i))/6;

    if y(i+1) < 0
        break;
    end

end

x = x(1:i+1);
y = y(1:i+1);
u = u(1:i+1);
v = v(1:i+1);
k1u = k1u(1:i+1);
k2u = k2u(1:i+1);
k3u = k3u(1:i+1);
k4u = k4u(1:i+1);
k1v = k1v(1:i+1);
k2v = k2v(1:i+1);
k3v = k3v(1:i+1);
k4v = k4v(1:i+1);
T = (i+1)*step;

disp(T)

end