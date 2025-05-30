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
ymax = 11;
ymin = 8;

t0 = 0;
tf = 10;
dt = 0.001;


x0 = 0;
y0 = 0;

psif = 0;
Ff = 0;

u0 = 8; %b1
v0 = 20; %b2
E = 5; %b3
M = 1.3;
count=1;

kill = false;
for p = 1:20 %ALM
    
    for k = 1:100
R = (0.75/pi)^(1/3) * (0.5/r1 + E/r2)^(1/3);
if imag(R) ~= 0
    kill = true;
    break;
end
A = pi * R^2;
m = 0.5+E;

F1 = @(u,v) -(2*m)^(-1)*rair*cd*A*sqrt(u^2 + v^2)*u;
F2 = @(u,v) -g - (2*m)^(-1)*rair*cd*A*sqrt(u^2 + v^2)*v;

[x,y,u,v,T,k1u,k1v,k2u,k2v,k3u,k3v,k4u,k4v] = RK4(F1,F2,x0,y0,u0,v0,dt,tf);
Tp(count)=T;
count = count+1;
[~, index] = min(abs(x - 10));         
ywall = y(index);             

B = -(rair * cd * A)/(2*m);

uprime = @(u,v) B * sqrt(u^2 + v^2)*u;
vprime = @(u,v) -g + B * sqrt(u^2 + v^2) * v;
vel = @(u,v) sqrt(u^2 + v^2);

F3 = @(Psi,F,Psiprime,Fprime,u,v,uprime,vprime,vel) -B*(u^2/vel + vel)*Psiprime - (B/vel)*(-(u*uprime + v*vprime)*(u^2/vel^2) + 3*u*uprime + v*vprime)*Psi - (B*u*v/vel)*Fprime - (B/vel)*(-(u*uprime + v*vprime)*(u*v/vel^2) + v*uprime + u*vprime)*F;
F4 = @(Psi,F,Psiprime,Fprime,u,v,uprime,vprime,vel) -B*(v^2/vel + vel)*Fprime - (B/vel)*(-(u*uprime + v*vprime)*(v^2/vel^2) + 3*v*vprime + u*uprime)*F - (B*u*v/vel)*Psiprime - (B/vel)*(-(u*uprime + v*vprime)*(u*v/vel^2) + v*uprime + u*vprime)*Psi;
F5 = @(Psi,F,u,v,vel) -((rair * cd * A * vel)/(2*m^2)) * (Psi * u + F * v);
F6 = @(Psi,F,u,v,vel) ((rair * cd * vel)/(2 * m * r2)) * (pi/6)^(1/3) * (0.5/r1 + E/r2)^(-1/3) * (Psi * u + F * v);

Psiprimef = 2*M*(x(end)-xtar);
% Fprimef = -(1+2*M*(x(end)-xtar)*u(end))/(v(end));  %Χωρις τοιχο
if ywall - R < ymin
    Fprimef = -(1+2*M*(x(end)-xtar)*u(end))/(v(end)) +   2 * M * (ywall - R - ymin);
elseif ywall + R > ymax
    Fprimef = -(1+2*M*(x(end)-xtar)*u(end))/(v(end)) + 2 * M * (ywall + R - ymax);
else
    Fprimef = -(1+2*M*(x(end)-xtar)*u(end))/(v(end));
end
[Psi,F,Psiprime,Fprime,Df3] = RK4_fae(F3,F4,psif,Ff,Psiprimef,Fprimef,dt,u,v,uprime,vprime,vel,F5,F6,k1u,k1v,k2u,k2v,k3u,k3v);

Df1 = -Psi(1);
Df2 = -F(1);

if E - (0.25-0.012*p) * Df3 < 0
    kill = true;
    break;
end
E = E - (0.25-0.012*p)  * Df3;
u0 = u0 - (0.25-0.012*p) * Df1;
v0 = v0 - (0.25-0.012*p) * Df2;
Ep(count)=E;
up(count)=u0;
vp(count)=v0;
    end
    if kill 
        break;
    end
    if M <30
        M = M*1.03;
    end
end

b = [u0 v0 E]
figure(1)
plot(x,y,'LineWidth',3)
grid on
xlabel('x [m]')
ylabel('y [m]')
title('Τροχιά μάζας')


function [Psi,F,Psiprime,Fprime,Df3] = RK4_fae(F1,F2,Psif,Ff,Psiprimef,Fprimef,step,u,v,uprime,vprime,vel,F5,F6,k1u,k1v,k2u,k2v,k3u,k3v)
    
    size = length(u);
    Psi = zeros(size,1);
    Psi(size) = Psif;
    F = zeros(size,1);
    F(size) = Ff;
    Psiprime = zeros(size,1);
    Psiprime(size) = Psiprimef; 
    Fprime = zeros(size,1);
    Fprime(size) = Fprimef; 
    Df3 = 0;


     for i = size:-1:2

        k1psi = -step * Psiprime(i);
        k1f = -step * Fprime(i);
        k1uprime = -step*uprime(u(i),v(i));
        k1vprime = -step*vprime(u(i),v(i));
        k1psiprime = -step * F1(Psi(i),F(i),Psiprime(i),Fprime(i),u(i),v(i),uprime(u(i),v(i)),vprime(u(i),v(i)),vel(u(i),v(i)));
        k1fprime = -step * F2(Psi(i),F(i),Psiprime(i),Fprime(i),u(i),v(i),uprime(u(i),v(i)),vprime(u(i),v(i)),vel(u(i),v(i)));

        k2psi = -step * (Psi(i) + k1psiprime/2); 
        k2f = -step * (Fprime(i) + k1fprime/2);
        k2uprime = -step*uprime(u(i)+k1uprime/2,v(i)+k1vprime/2);
        k2vprime = -step*vprime(u(i)+k1uprime/2,v(i)+k1vprime/2);
        k2psiprime = -step * F1(Psi(i) + k1psi/2,F(i)+k1f/2,Psiprime(i)+k1psiprime/2,Fprime(i)+k1fprime/2,u(i)+k1u(i),v(i)+k1v(i),uprime(u(i)+k1uprime/2,v(i)+k1vprime/2),vprime(u(i)+k1uprime/2,v(i)+k1vprime/2),vel(u(i)+k1uprime/2,v(i)+k1vprime/2));
        k2fprime = -step * F2(Psi(i)+k1psi/2,F(i)+k1f/2,Psiprime(i)+k1psiprime/2,Fprime(i)+k1fprime/2,u(i)+k1u(i),v(i)+k1v(i),uprime(u(i)+k1uprime/2,v(i)+k1vprime/2),vprime(u(i)+k1uprime/2,v(i)+k1vprime/2),vel(u(i)+k1uprime/2,v(i)+k1vprime/2));

        k3psi = -step * (Psi(i) + k2psiprime/2); 
        k3f = -step * (Fprime(i) + k2fprime/2);
        k3uprime = -step*uprime(u(i)+k2uprime/2,v(i)+k2vprime/2);
        k3vprime = -step*vprime(u(i)+k2uprime/2,v(i)+k2vprime/2);
        k3psiprime = -step * F1(Psi(i) + k2psi/2,F(i)+k2f/2,Psiprime(i)+k2psiprime/2,Fprime(i)+k2fprime/2,u(i)+k2u(i),v(i)+k2v(i),uprime(u(i)+k2uprime/2,v(i)+k2vprime/2),vprime(u(i)+k2uprime/2,v(i)+k2vprime/2),vel(u(i)+k2uprime/2,v(i)+k2vprime/2));
        k3fprime = -step * F2(Psi(i)+k2psi/2,F(i)+k2f/2,Psiprime(i)+k2psiprime/2,Fprime(i)+k2fprime/2,u(i)+k2u(i),v(i)+k2v(i),uprime(u(i)+k2uprime/2,v(i)+k2vprime/2),vprime(u(i)+k2uprime/2,v(i)+k2vprime/2),vel(u(i)+k2uprime/2,v(i)+k2vprime/2));

        k4psi = -step * (Psi(i) + k3psiprime); 
        k4f = -step * (Fprime(i) + k3fprime);
        k4psiprime = -step * F1(Psi(i) + k3psi,F(i)+k3f,Psiprime(i)+k3psiprime,Fprime(i)+k3fprime,u(i)+k3u(i),v(i)+k3v(i),uprime(u(i)+k3uprime/2,v(i)+k3vprime/2),vprime(u(i)+k3uprime/2,v(i)+k3vprime/2),vel(u(i)+k3uprime/2,v(i)+k3vprime/2));
        k4fprime = -step * F2(Psi(i)+k3psi,F(i)+k3f,Psiprime(i)+k3psiprime,Fprime(i)+k3fprime,u(i)+k3u(i),v(i)+k3v(i),uprime(u(i)+k3uprime/2,v(i)+k3vprime/2),vprime(u(i)+k3uprime/2,v(i)+k3vprime/2),vel(u(i)+k3uprime/2,v(i)+k3vprime/2));

        Psi(i-1) = Psi(i) + (k1psi + 2 * k2psi + 2 * k3psi + k4psi)/6;
        F(i-1) = F(i) + (k1f + 2 * k2f + 2 * k3f + k4f)/6;
        Psiprime(i-1) = Psiprime(i) + (k1psiprime + 2 * k2psiprime + 2 * k3psiprime + k4psiprime)/6;
        Fprime(i-1) = Fprime(i) + (k1fprime + 2 * k2fprime + 2 * k3fprime + k4fprime)/6;

        Df3 = Df3 + step * (F5(Psi(i),F(i),u(i),v(i),vel(u(i),v(i))) + F5(Psi(i-1),F(i-1),u(i-1),v(i-1),vel(u(i-1),v(i-1))))/2 + step * (F6(Psi(i),F(i),u(i),v(i),vel(u(i),v(i))) + F6(Psi(i-1),F(i-1),u(i-1),v(i-1),vel(u(i-1),v(i-1))))/2;
     end

%disp(Df3)

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