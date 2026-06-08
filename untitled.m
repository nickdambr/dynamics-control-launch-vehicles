clc
clear all
close all

mu_c =6.0;
mu_alfa = 4.0;

%% LTI system

s = tf('s');

G = mu_c / ( s^2 -mu_alfa);

G2 = [mu_c/( s^2 - mu_alfa)  s*mu_c/(s^2-mu_alfa)];

num = [mu_c];
den = [1 0 -mu_alfa];
G1 = tf(num, den);

%zpk

% Space-State representation 
% EXAMPLE: x'' = u, y=x
A= [0 1
    0 0];
B= [0; 1];
C = [1 0];
D = 0;

sys = ss(A, B, C, D);


% LV1dof th'' = mu_alfa * th + mu_c * delta
% x = th
% u = delta
% x' =        0 * x + 1 * x' +    0 * u
% x'' = mu_alfa * x + 0 * x' + mu_c * u
% y = x


A= [0 1
    mu_alfa 0];
B= [0; mu_c];
C = [1 0];
D = 0;

sys_LV = ss(A,B,C,D);
sys_LV.StateName = {'th', 'thdot'};
sys_LV.InputName = {'delta'};
sys_LV.OutputName = {'th_{meas}'};

G_LV = tf(sys_LV);

% Controller

Kp = 1;
Kd = 0.1;


G_ctrl = Kp + s * Kd;


% Open_Loop

L = G_LV * G_ctrl;



%% Frequency-Domain tools


% Display the bode plot
figure()
bode(L)


%figure()
%nyquistplot(L)

figure()
hold on
nichols(L)


num = [mu_c * 0.7];
den = [1 0 -mu_alfa*1.3];
G_LV2 = tf(num,den);

L2 = G_LV2 * G_ctrl;

nichols(L2)
legend('G','G2')
hold on

G_ctrl3 = (Kp + s*Kd) * 2;


L4 = G_LV *G_ctrl_rob;
L3 = G_LV2 * G_ctrl3;

nichols(L3)
legend('G_nom, K_nom','G2, K_nom','G_2_robust')