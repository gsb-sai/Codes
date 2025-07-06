function [t,z] = cascade_UV()
tspan = 0:0.5:200;
z0 = [0,0,0,0,0,0,0];

%T = 0.1;
 U = 0.4;
 V = 0.5;
kd1 = 0.3;
kd2 = 0.9;
kd3 = 0.86;
kd4 = 0.78;
kd5 = 0.6; %0.66
kd6 = 0.8; %0.8
kd7 = 0.7; %0.7

k1 = 0.6; km1 = 0.15; 
k2 = 0.85; km2 = 0.04; 
k3 = 0.75; km3 = 0.049; 
k4 = 0.66; km4 = 0.3; 
k5 = 0.83; km5 = 0.43; 
k6 = 0.89; km6 = 0.39; %0.8
k7 = 0.69; km7 = 0.3; %0.7
k8 = 0.87; km8 = 0.19; %0.6
k9 = 0.92; km9 = 0.2; % 0.7
n =2;

function fa = s1(t)
%if (t >= 2 && t < 4.5) || (t >= 10 && t < 12) || (t >= 18 && t < 20) || (t >= 26 && t < 28) || (t >= 60 && t < 100)|| (t >= 130 && t < 170)
if (t>=2 && t<=60) %|| (t >= 8 && t <= 10) || (t >= 13 && t <= 15) || (t >= 18 && t <= 20)  || (t >= 30 && t<=70)
    fa = 1;
else 
    fa = 0;
end
end

function fb = s2(t)
%if (t >= 6.5 && t < 8.5) || (t >= 14.5 && t < 16.5) || (t >= 22 && t < 28) || (t >= 32. && t < 34) || (t >= 40 && t < 80)|| (t >= 110 && t < 170)
if (t>=2 && t<=60) %|| (t >= 8 && t <= 10) || (t >= 13 && t <= 15) || (t >= 18 && t <= 20)  || (t >= 30 && t<=70)
    fb = 1;
else 
    fb = 0;
end
end

figure(1);
subplot(9,1,1)
fplot(@(t) s1(t),'k','LineWidth',1.2);hold on
xlim([0 100])
ylim([0 1.5])
ylabel('U','FontSize',12,'FontName','Times New Roman')
figure(1);
subplot(9,1,2)
fplot(@(t) s2(t),'k','LineWidth',1.2);hold on
xlim([0 100])
ylim([0 1.5])
ylabel('V','FontSize',12,'FontName','Times New Roman')
[t,z] = ode23s(@cascade_UV,tspan,z0);
subplot(9,1,3)
plot(t,z(:,1),'k','LineWidth',1.2);
xlim([0 100])
%xlabel('Time,(A.U)','FontSize',12,'FontName','Times New Roman')
ylabel('[X1],(A.U)','FontSize',12,'FontName','Times New Roman')
subplot(9,1,4)
plot(t,z(:,2),'k', 'LineWidth',1.2);
xlim([0 100])
%xlabel('Time,(A.U)','FontSize',12,'FontName','Times New Roman')
ylabel('[Y1],(A.U)','FontSize',12,'FontName','Times New Roman')
subplot(9,1,5)
plot(t,z(:,3),'k', 'LineWidth',1.2);
xlim([0 100])
%xlabel('Time,(A.U)','FontSize',12,'FontName','Times New Roman')
ylabel('[Z1],(A.U)','FontSize',12,'FontName','Times New Roman')
subplot(9,1,6)
plot(t,z(:,4),'k', 'LineWidth',1.2);
xlim([0 100])
%xlabel('Time,(A.U)','FontSize',12,'FontName','Times New Roman')
ylabel('[Y2],(A.U)','FontSize',12,'FontName','Times New Roman')
subplot(9,1,7)
plot(t,z(:,5),'k', 'LineWidth',1.2);
xlim([0 100])
%xlabel('Time,(A.U)','FontSize',12,'FontName','Times New Roman')
ylabel('[Z2],(A.U)','FontSize',12,'FontName','Times New Roman')
subplot(9,1,8)
plot(t,z(:,6),'k', 'LineWidth',1.2);
xlim([0 100])
%xlabel('Time,(A.U)','FontSize',12,'FontName','Times New Roman')
ylabel('[Y3],(A.U)','FontSize',12,'FontName','Times New Roman')
subplot(9,1,9)
plot(t,z(:,7),'k', 'LineWidth',1.2);
xlim([0 100])
xlabel('Time,(A.U)','FontSize',12,'FontName','Times New Roman')
ylabel('[Z3],(A.U)','FontSize',12,'FontName','Times New Roman')


function f = cascade_UV(t,z)
% z1 = x1; z2 = y1; z3 = z1; z4 = y2; z5 = z2; z6 = y3; z7 = z3;
f(1) = U*s1(t) - kd1*z(1); 
f(2) = V*s2(t) + (k1*z(1)^n)/(km1^n + z(1)^n) - kd2*z(2); 
f(3) = ((k2*z(2)^n)/(km2^n + z(2)^n)) * ((k3*(z(1)^n)/(km3^n + z(1)^n))) - kd3*z(3); % AND Logic
f(4) = (k4*z(3)^n)/(km4^n + z(3)^n) - kd4*z(4); 
f(5) =  ((k5*z(4)^n)/(km5^n + z(4)^n)) * ((k6*(z(3)^n)/(km6^n + z(3)^n))) - kd5*z(5); % AND Logic
f(6) = (k7*z(5)^n)/(km7^n + z(5)^n) - kd6*z(6); 
f(7) =  ((k8*z(6)^n)/(km8^n + z(6)^n)) * ((k9*(z(5)^n)/(km9^n + z(5)^n))) - kd7*z(7); % AND Logic
f= f';

end 
end