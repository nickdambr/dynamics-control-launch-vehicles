clc
clear all

mu_c = 6;
mu_alfa = 4;

Kp = 2;
Kd = 1;

s = tf('s');

% th'' = mu:c é delta + mu_alfa * th


A =[0 1
    mu_alfa 0
    ]
B = [0;
    mu_c]

C = [1,0] % y = th

D = [0]
    

sysLV = ss(A,B,C,D);

sysLV.InputName = {'delta'};
sysLV.OutputName = {'th_meas'};
sysLV.StateName = {'th', 'thdot'};

Krigid = Kp + Kd*s; % Krigid(s) = delta(s)/e(s)
Krigid.InputName  = {'th_err'};
Krigid.OutputName = {'delta'};

% Closed Loop system

err_name = {'th_err'};
ref_name = {'th_ref'};
meas_name = {'th_meas'};

sumJunction = sumblk('%s = %s - %s', err_name, ref_name, meas_name);

T0 = connect(sumJunction, Krigid, sysLV, ref_name, {'th_meas','delta'}, {'delta'});

L = getLoopTransfer(T0, 'delta', -1)

figure(1)
hold on
grid on
nicholsplot(L)