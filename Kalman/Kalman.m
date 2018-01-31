

clc
clear
close all

%% user defined

Fs=1;
    
N=1e2;

sigma_x=0.5;

sigma_z=0.1;

n=4; % state dimension
l=1; % input dimension
m=2; % measurement dimension

v_x=10;
a_x=0.2; % acceleration

v_y=2; 
a_y=0.1;


%% Intermediate calculations

Ts=1/Fs;

k=0:Ts:(N-1)*Ts;

Q=(sigma_x^2)*eye(n);
R=(sigma_z^2)*eye(m);

%% Simulation

%%% Data Generation

x_GT=zeros(1,N);
x_dot_GT=v_x*ones(1,N);
y_GT=zeros(1,N);
y_dot_GT=v_y*ones(1,N);

for kk=2:N
    
    x_GT(kk)=x_GT(kk-1)+x_dot_GT(kk-1)*Ts+a_x*(Ts^2)/2;
    x_dot_GT(kk)=x_dot_GT(kk-1)+a_x*Ts;
    
%     y_GT(kk)=y_GT(kk-1)+y_dot_GT(kk-1)*Ts;
%     y_dot_GT(kk)=y_dot_GT(kk-1);

    y_GT(kk)=y_GT(kk-1)+y_dot_GT(kk-1)*Ts+a_y*(Ts^2)/2;
    y_dot_GT(kk)=y_dot_GT(kk-1)+a_y*Ts;
    
end

X_GT=[x_GT;y_GT;x_dot_GT;y_dot_GT];

X_Noised=X_GT+sigma_x*randn(n,N);

Z_GT=[x_GT;y_GT];

Z_Noised=Z_GT+sigma_z*randn(m,N);

%%% Kalman filter 

u=[a_x;a_y]; % input vector

F=[1 0 Ts 0;0 1 0 Ts;0 0 1 0; 0 0 0 1];
G=[Ts^2/2 0; 0 Ts^2/2; Ts 0; 0 Ts];

%G=[Ts^2/2; 0; Ts; 0];
H=[1 0 0 0;0 1 0 0];

%x_est_apr=zeros(n,N);
X_est_apr=zeros(n,N);
x_est_postr=zeros(n,1);

P_est_postr=eye(n);

for kk=2:N

     %%% 1) Time update equations
    
    x_est_apr=F*x_est_postr+G*u;
    
    X_est_apr(:,kk)= x_est_apr;
    
    P_est_apr=F*P_est_postr*F'+Q;
    
    %%%% 2) Measurement correction equations
    
    K=P_est_apr*H'*inv(H*P_est_apr*H'+R);
    
    z=Z_Noised(:,kk);
    
    x_est_postr=x_est_apr+K*(z-H*x_est_apr);
    
    P_est_postr=(eye(n)-K*H)*P_est_apr;
  
end

%% Display

figure
set(gcf,'windowstyle','docked')


subplot(2,2,1)
plot(X_GT(1,:))
hold on
plot(X_Noised(1,:))
hold on
plot(X_est_apr(1,:))
title('x axis location')
grid minor
xlabel('time [Ts]')
ylabel('Location [m]')
legend('Ground Truth','Noised','Predicted','Location','best')

subplot(2,2,2)
plot(X_GT(2,:))
hold on
plot(X_Noised(2,:))
hold on
plot(X_est_apr(2,:))
title('y axis location')
grid minor
xlabel('time [Ts]')
ylabel('Location [m]')
legend('Ground Truth','Noised','Predicted','Location','best')

subplot(2,2,3)
plot(X_GT(3,:))
hold on
plot(X_Noised(3,:))
hold on
plot(X_est_apr(3,:))
title('x axis velocity')
grid minor
xlabel('time [Ts]')
ylabel('Velocity [m/sec]')
legend('Ground Truth','Noised','Predicted','Location','best')

subplot(2,2,4)
plot(X_GT(4,:))
hold on
plot(X_Noised(4,:))
hold on
plot(X_est_apr(4,:))
title('y axis velocity')
grid minor
xlabel('time [Ts]')
ylabel('Velocity [m/sec]')
legend('Ground Truth','Noised','Predicted','Location','best')











