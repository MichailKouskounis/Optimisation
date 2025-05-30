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
for p = 1:20 %ALM
    
    for k = 1:500
R = (0.75/pi)^(1/3) * (0.5/r1 + E/r2)^(1/3);
if imag(R) ~= 0
    kill = true;
    break;
end
A = pi * R^2;
m = 0.5+E;

F1 = @(vx,vy) -(2*m)^(-1)*rair*cd*A*sqrt(vx^2 + vy^2)*vx;
F2 = @(vx,vy) -g - (2*m)^(-1)*rair*cd*A*sqrt(vx^2 + vy^2)*vy;

[x,y,u,v,T] = RK4(F1,F2,x0,y0,u0,v0,dt,tf);

[Df1,Df2] =  cvm(F1,F2,x0,y0,u0,v0,dt,tf,epsilon,M,R);

E = complex(E,epsilon);
R = (0.75/pi)^(1/3) * (0.5/r1 + E/r2)^(1/3);
A = pi * R^2;
m = 0.5+E;
F1 = @(vx,vy) -(2*m)^(-1)*rair*cd*A*sqrt(vx^2 + vy^2)*vx;
F2 = @(vx,vy) -g - (2*m)^(-1)*rair*cd*A*sqrt(vx^2 + vy^2)*vy;

[x5,y5,~,~,T5] = RK4(F1,F2,x0,y0,u0,v0,dt,tf);

[~, index] = min(abs(x5 - 10));         
ywall = y5(index);             

if ywall - R < 8
    Df3 = imag((T5 + M*((x5(end)-xtar)^2) + M*(ywall-R-8)^2 + y5(end)))/(epsilon);
elseif ywall + R > 11
    Df3 = imag((T5 + M*((x5(end)-xtar)^2) + M*(ywall+R-11)^2 + y5(end)))/(epsilon);
else
    Df3 = imag((T5 + M*((x5(end)-xtar)^2) + y5(end)))/(epsilon);
end
fprintf('x(T) = %f + %fi p = %d k = %d \n',real(x(end)),imag(x(end)),p,k)
E = real(E);

if E - 0.01 * Df3 < 0
    kill = true;
    break;
end
u0 = u0 - 0.01 * Df1;
v0 = v0 - 0.01 * Df2;
E = E - 0.01 * Df3;
    end
    if kill 
        break;
    end
    if M <30
        M = M*1.02;
    end
end



function [df1,df2] = cvm(F1,F2,x0,y0,u0,v0,step,tf,epsilon,M,R)
    df1 = 0;
    df2 = 0;

    [x1,y1,~,~,T1] = RK4(F1,F2,x0,y0,complex(u0,epsilon),v0,step,tf);
    [~, index1] = min(abs(x1 - 10)); 
    ywall = y1(index1);             
 
    if ywall - R < 8
       df1 = imag((T1 + M*((x1(end)-20)^2) + M*(ywall-R-8)^2 + y1(end)))/(epsilon);
    elseif ywall + R > 11
        df1 = imag((T1 + M*((x1(end)-20)^2) + M*(ywall+R-11)^2 + y1(end)))/(epsilon);
    else
        df1 = imag((T1 + M*((x1(end)-20)^2) + y1(end)))/(epsilon);
    end

    [x2,y2,~,~,T2] = RK4(F1,F2,x0,y0,u0,complex(v0,epsilon),step,tf);
    [~, index2] = min(abs(x2 - 10)); 
    ywall = y2(index2);             

    if ywall - R < 8
       df2 = imag((T2 + M*((x2(end)-20)^2) + M*(ywall-R-8)^2 + y2(end)))/(epsilon);
    elseif ywall + R > 11
        df2 = imag((T2 + M*((x2(end)-20)^2) + M*(ywall+R-11)^2 + y2(end)))/(epsilon);
    else
        df2 = imag((T2 + M*((x2(end)-20)^2) + y2(end)))/(epsilon);
    end
   
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

% disp(T)

end