%% MAE 200 Controls Double Inverted Pendulum
%dynamics
% x = [x th1 th2 xDot th1Dot th2Dot]
l1 = 1;
l2 = .5;
mc = 10;
m1 = 1;
m2 = 0.5;
g = 9.81;
I1 = 1/3 *m1*l1^2;
I2 = 1/3 * m2*l2^2;
T2 = 15; %total time to run simulation

%% Control of cart during phase 2 (stablization)
% x = [x th1 th2 xDot th1Dot th2Dot]

E = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 mc+m1+m2 -m1*l1 -m2*l2;
     0 0 0 -m1*l1 I1+m1*l1^2 0;
     0 0 0 -m2*l2 0 I2+m2*l2^2;];
N = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0;
     0 m1*g*l1 0 0 0 0;
     0 0 m2*g*l2 0 0 0;];

A = E^-1*N;
B = (E^-1*[0 0 0 1 0 0]');

C = diag([1 1 1 0 0 0]);
D = [0 0 0 0 0 0]';

Qc = diag([100 1 1 1 100 100]);  % State weight matrix for energy function

Qe = diag([1 1 1 1 1 1]);  % Spectral density for disturbance
   
Rc = 3;     % Control authority weight for Energy function
Re = 9;     % Spectral density for noise






X= icare(A,B,Qc,Rc);     % X = Ricatti Equation Energy value
[P,Le] = icare(A',C,Qe,Re);  % P = covariance from Ricatti Eqn

K = -Rc^-1*B'*X;    %Calculate optimal K
L = -P*C'*Qe^-1;    %Calculate optimal L
x=zeros(6,length(t));
x0=x_k_6(:,end);

x_hat=zeros(6,length(t));
x = RK4_x_LQR(A,B,K,T2,x0); % March LQR equation using RK4
xhat=RK4_xhat_LQG(A,B,C,K,Le,T2,x);     %March LQG equation using RK4

xtotalLQG = [x_hat_sim, xhat];  % Append states from Phase A to Phase B LQG
xtotalLQR = [xsim, x];  %Append states from phase A to phase B LQR

%% Plotting
sysFS = ss(A,B,C,D);
sysKF = ss((A-Le'*C),[B Le'],eye(6),0*[D Le']);

clf

t = 0:size(xtotalLQG,2)-1;
figure(1)
hold on
for i = 1:6
plot(t/1000,xtotalLQG(i,:))
end
legend('x (m)','theta1 (rad)','theta2 (rad)','xdot (m/s)','theta1dot (rad/s)','theta2dot (rad/s)');
xlabel('time (s)')
title('LQG Steady State Errors vs Time for Full System');

t = 0:size(xtotalLQR,2)-1;
figure(2)
hold on
for i = 1:6
plot(t/1000,xtotalLQR(i,:))
end
legend('x (m)','theta1 (rad)','theta2 (rad)','xdot (m/s)','theta1dot (rad/s)','theta2dot (rad/s)');
xlabel('time (s)')
title('LQR Steady State Errors vs Time');

figure(3)
hold on
t=0:size(x_hat_sim,2)-1;
for i = 1:6
plot(t/1000,x_hat_sim(i,:))
plot(t/1000,x_k_6)
end
legend('x (m)','theta1 (rad)','theta2 (rad)','xdot (m/s)','theta1dot (rad/s)','theta2dot (rad/s)');
xlabel('time (s)')
title('LQG Steady State Errors  vs Time Swing Up phase ');













