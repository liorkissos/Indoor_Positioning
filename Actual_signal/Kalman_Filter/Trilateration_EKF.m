
%%% based on the article "Ultra-Wideband-Based Localization for
%%% Quadcopter Navigation" by Kexin Guo & Zhirong Qiu
%%% Relevant material:
%%% 1) DecaRange ARM source files
%%% 2) DecaRangeRTLS_ARM_Source_Code_Guide.pdf
%%% 3) based on the article "Ultra-Wideband-Based Localization for Quadcopter Navigation" by Kexin Guo & Zhirong Qiu

clc
clear
close all

%% User defined parameters


%%% 1) System model- Noises

% the bigger the state eqaution noise, the more we rely on the measurements
% and thus get closer to pure Trilateration.

% Measurement equation noise: standard deviation (v does not mean velocity)
sigma_v=1.5; % [m^2].
% State equation noise: standard deviation of the acceleration.  Guo&Qiu (page 5 left): maximum possible value
sigma_acc_x=0.1; % [m/sec^2].  acceleration
sigma_acc_y=0.1; % [m/sec^2].  acceleration

% % Measurement equation noise: standard deviation (v does not mean velocity)
% sigma_v=0.2; % [m^2].
% % State equation noise: standard deviation of the acceleration.  Guo&Qiu (page 5 left): maximum possible value
% sigma_acc_x=5.0; % [m/sec^2].  acceleration
% sigma_acc_y=5.0; % [m/sec^2].  acceleration

% % Measurement equation noise: standard deviation (v does not mean velocity)
% sigma_v=1.5; % [m^2].
% % State equation noise: standard deviation of the acceleration.  Guo&Qiu (page 5 left): maximum possible value
% sigma_acc_x=20.0; % [m/sec^2].  acceleration
% sigma_acc_y=20.0; % [m/sec^2].  acceleration

% problem dimensions
n=4; % state vector dimension: position x axis, position y axis, velocity x axis, velocity y axis
m=3; % measurement vector dimension: square distance from 3 anchors

%%% 2) Anchors positioning
% X-Y [meters]
p_anch_1_x=0;
p_anch_1_y=1.8;

p_anch_2_x=3;
p_anch_2_y=1.2;

p_anch_3_x=0;
p_anch_3_y=0;


% Z (height) [meters]

% Anchors
h1=2.28;
h2=2.28;
h3=2.28;

h_tag=1.43; % tag


%%% 3) Decawave configuration
Fs=1/280e-3; % the time interval between samples is 280msec. % parameter retrieved from ARM source code syscalls.c (portGetTickCount) and compiler.h (CLOCKS_PER_SEC)


%%  Recoreded file parsing

Record=readtable('col.txt');

%%% Ranges
R1_3D_mm=hex2dec(table2cell(Record(:,3)));
R2_3D_mm=hex2dec(table2cell(Record(:,4)));
R3_3D_mm=hex2dec(table2cell(Record(:,5)));

% conversion to  2D via Pythagoras
R1_2D_mm=sqrt(R1_3D_mm.^2-((h1-h_tag)*1e3).^2);
R2_2D_mm=sqrt(R2_3D_mm.^2-((h2-h_tag)*1e3).^2);
R3_2D_mm=sqrt(R3_3D_mm.^2-((h3-h_tag)*1e3).^2);


% mm's to meters
R1_2D_m=R1_2D_mm/1e3;
R2_2D_m=R2_2D_mm/1e3;
R3_2D_m=R3_2D_mm/1e3;

%% Intermediate calculations

%%% Measurements vector
Z_Noised=[R1_2D_m'.^2;R2_2D_m'.^2;R3_2D_m'.^2]; % measurement vector is the square horizontal distance from the 3 anchors which is regarded as the noisy measurement

% Time
t_ticks=hex2dec(table2cell(Record(:,9))); % the DEBUG row see (DecaRangeRTLS_ARM_source.pdf)
t_ticks=t_ticks-t_ticks(2); % the record starts to be regular only from the 2nd sample and on
Ts=1/Fs; % the reading time interval in seconds

t=t_ticks*Ts; % the time axis in seconds

% length of recording
N=size(Record,1);


%% Kalman filter
%%%%%%%%%%%%

%% Kalman problem formulation and Initialization

%%% Initializations

% 0) Calculating the initial guess: Guo& Qiu article via linear LS problem (UWB based
% localization article.) equation (15)

A_LS=2*[p_anch_1_x-p_anch_2_x , p_anch_1_y-p_anch_2_y ;...
    p_anch_1_x-p_anch_3_x , p_anch_1_y-p_anch_3_y];

d1=R1_2D_m(1);
d2=R2_2D_m(1);
d3=R3_2D_m(1);

m1=d2^2-d1^2-(p_anch_2_x^2-p_anch_1_x^2+p_anch_2_y^2-p_anch_2_y^2);
m2=d3^2-d1^2-(p_anch_3_x^2-p_anch_1_x^2+p_anch_3_y^2-p_anch_1_y^2);

p_tag_init=A_LS\[m1;m2]; % LS solution

% 1) State vector
X_est_apr=zeros(n,N);

x_est_postr=[p_tag_init;0;0]; %we assume the initial velocities are zero
X_est_postr=[x_est_postr,zeros(n,N-1)];

% 2) Covariance matrix of
P_est_postr=eye(n); % estimation error Covariance matrix

% 3)  Anchors position

p_anch_1=[p_anch_1_x;p_anch_1_y];
p_anch_2=[p_anch_2_x;p_anch_2_y];
p_anch_3=[p_anch_3_x;p_anch_3_y];

%%% Kalman equations

% 1) State equation

F=[1 0 Ts 0; 0 1 0 Ts;0 0 1 0;0 0 0 1]; % state equation matrix

f=@(x,w) F*x+w; % we model acceleration as noise, hence we do not need the input vector;  f=@(x,u,w) F*x+G*u+w;

%%% 2) Measurement equation
h=@(p,p_anch,v) (p-p_anch)'*(p-p_anch)+v; % p and p_anch are 2x1 column vectors

%%% Jacobians

% 1) current state wrt to former state
A=F;
% 2) current state wrt current state noise
W=eye(n);

% 3) current measurement wrt current state
p_anch_x=[p_anch_1_x;p_anch_2_x;p_anch_3_x];
p_anch_y=[p_anch_1_y;p_anch_2_y;p_anch_3_y];

H=2*[x_est_postr(1)*ones(m,1)-p_anch_x,...
    x_est_postr(2)*ones(m,1)-p_anch_y,...
    zeros(m,1),zeros(m,1)];

% 4) current measurment wrt current measurement noise
V=eye(m);

%%% Covariance matrices

%%% Model noises covariance matrices

% 1)  The state equation noise covariance matrix

% the Q(1,1) (and Q(2,1)) need to be sigma_acc_x^2*Ts^4/4
% but sometimes it does not converge unless we change to sigma_acc_x^2*Ts^4/2

Q=[sigma_acc_x^2*Ts^4/4, 0 , sigma_acc_x^2*Ts^3/2, 0;...
    0, sigma_acc_y^2*Ts^4/4, 0 , sigma_acc_y^2*Ts^3/2;...
    sigma_acc_x^2*Ts^3/2    , 0 , sigma_acc_x^2*Ts^2, 0;...
    0, sigma_acc_y^2*Ts^3/2, 0, sigma_acc_y^2*Ts^2         ]; % equation (7)

%%% Do not erase the following 4 lines!
% Q=[(sigma_acc_x^2)*(Ts^4)/2, 0 , sigma_acc_x^2*Ts^3/2, 0;...
%     0, (sigma_acc_y^2)*(Ts^4)/2, 0 , sigma_acc_y^2*Ts^3/2;...
%     sigma_acc_x^2*Ts^3/2    , 0 , sigma_acc_x^2*Ts^2, 0;...
%     0, sigma_acc_y^2*Ts^3/2, 0, sigma_acc_y^2*Ts^2         ]; % equation (7)

% 2) Measurement equation noise covariance ,atrix

R=(sigma_v^2)*eye(m);

%% Kalman Iteration

for kk=2:N
    
    %%% 1) Time update equations
    
    x_est_apr=f(x_est_postr,0); % apriori state vector estimator
    
    P_est_apr=A*P_est_postr*A'+W*Q*W'; % apriori estimation error covariance matrix
    
    %%% 2) Measurement correction equations
    
    H=2*[x_est_postr(1)*ones(m,1)-p_anch_x,... % Jacobian of measurement wrt state calculation
        x_est_postr(2)*ones(m,1)-p_anch_y,...
        zeros(m,1),zeros(m,1)];
    
    K=P_est_apr*H'/(H*P_est_apr*H'+V*R*V); % the K factor (K=P_est_apr*H'*inv(H*P_est_apr*H'+V*R*V);)
    
    z=Z_Noised(:,kk);
    h_vec=[h(x_est_apr(1:2),p_anch_1,0);h(x_est_apr(1:2),p_anch_2,0);h(x_est_apr(1:2),p_anch_3,0)];
    
    x_est_postr=x_est_apr+K*(z-h_vec); % posteriori state vector estimator
    
    P_est_postr=(eye(n)-K*H)*P_est_apr; % posteriori estimation error estimator
    
    %%% 3) Estimators saving
    
    X_est_apr(:,kk)= x_est_apr;
    X_est_postr(:,kk)= x_est_postr;
    
end

%% Display


figure
set(gcf,'windowstyle','docked')
%plot(X_est_apr(1,:),X_est_apr(2,:))
plot(X_est_postr(1,:),X_est_postr(2,:))
title({['Tag Positioning estimation: EKF'],[' \sigma^v=',num2str(sigma_v),'[m^2]'],[' \sigma^{acc}=',num2str(sigma_acc_x),'[m/sec^2]  ']})
grid minor
xlabel('x[meter]')
ylabel('y[meter]')


% subplot(2,2,1)
% % plot(X_GT(1,:))
% % hold on
% % plot(X_Noised(1,:))
% % hold on
% plot(X_est_apr(1,:))
% title('x axis location')
% grid minor
% xlabel('time [Ts]')
% ylabel('Location [m]')
% legend('Ground Truth','Noised','Predicted','Location','best')
%
% subplot(2,2,2)
% % plot(X_GT(2,:))
% % hold on
% % plot(X_Noised(2,:))
% % hold on
% plot(X_est_apr(2,:))
% title('y axis location')
% grid minor
% xlabel('time [Ts]')
% ylabel('Location [m]')
% legend('Ground Truth','Noised','Predicted','Location','best')
%
% subplot(2,2,3)
% % plot(X_GT(3,:))
% % hold on
% % plot(X_Noised(3,:))
% % hold on
% plot(X_est_apr(3,:))
% title('x axis velocity')
% grid minor
% xlabel('time [Ts]')
% ylabel('Velocity [m/sec]')
% legend('Ground Truth','Noised','Predicted','Location','best')
%
% subplot(2,2,4)
% % plot(X_GT(4,:))
% % hold on
% % plot(X_Noised(4,:))
% % hold on
% plot(X_est_apr(4,:))
% title('y axis velocity')
% grid minor
% xlabel('time [Ts]')
% ylabel('Velocity [m/sec]')
% legend('Ground Truth','Noised','Predicted','Location','best')




shg








