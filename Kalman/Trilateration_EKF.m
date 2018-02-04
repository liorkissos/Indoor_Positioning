
%%%% based on the article "Ultra-Wideband-Based Localization for
%%%% Quadcopter Navigation" by Kexin Guo & Zhirong Qiu
clc
clear
close all

%% user defined


%Fs=1;
Fs=1/280e-3; % the time interval between samples is 280msec


% System model

sigma_v=0.3; % measurement equiation noise standard deviation (v does not mean velocity)

acc_x=0.05; % acceleration
acc_y=0.05;


% Anchors positioning
p_anch_1_x=0;
p_anch_1_y=1.8;

p_anch_2_x=3;
p_anch_2_y=1.2;

p_anch_3_x=0;
p_anch_3_y=0;


%%% height
h_tag=1.43; % tag

% Anchors
h1=2.28;
h2=2.28;
h3=2.28;

%%% dimensions
n=4; % state dimension
l=2; % input dimension
m=3; % measurement dimension

%% Parse the recoreded file

Record=readtable('col.txt');

%%% Ranges
R1_3D_mm=hex2dec(table2cell(Record(:,3)));
R2_3D_mm=hex2dec(table2cell(Record(:,4)));
R3_3D_mm=hex2dec(table2cell(Record(:,5)));

% conversion to  2D via Pythagoras
R1_2D_mm=sqrt(R1_3D_mm.^2-((h1-h_tag)*1e3).^2);
R2_2D_mm=sqrt(R2_3D_mm.^2-((h2-h_tag)*1e3).^2);
R3_2D_mm=sqrt(R3_3D_mm.^2-((h3-h_tag)*1e3).^2);

% R1_2D_m=R1_2D_mm/1e3-0.24;
% R2_2D_m=R2_2D_mm/1e3-0.24;
% R3_2D_m=R3_2D_mm/1e3-0.12;

R1_2D_m=R1_2D_mm/1e3;
R2_2D_m=R2_2D_mm/1e3;
R3_2D_m=R3_2D_mm/1e3;

Z_Noised=[R1_2D_m'.^2;R2_2D_m'.^2;R3_2D_m'.^2];

% Time
t_ticks=hex2dec(table2cell(Record(:,9))); % the DEBUG row see (DecaRangeRTLS_ARM_source.pdf)
t_ticks=t_ticks-t_ticks(2); % the record starts to be regular only from the 2nd sample and on
T_tick=1/Fs; % the reading time interval in seconds

t=t_ticks*T_tick; % the time axis in seconds

%
N=size(Record,1);


%% Intermediate calculations

Ts=1/Fs;

k=0:Ts:(N-1)*Ts;

% switch Model

sigma_acc_x=2*acc_x; % page 5 left. maximum possible value
sigma_acc_y=2*acc_y; % page 5 left. maximum possible value

%%%  The state equation noise covariance. the Q(1,1) needs to be sigma_acc_x^2*Ts^4/4
%%% but it does not converge unless we change to sigma_acc_x^2*Ts^4/2

%                 Q=[sigma_acc_x^2*Ts^4/4, 0 , sigma_acc_x^2*Ts^3/2, 0;...
%                       0, sigma_acc_y^2*Ts^4/4, 0 , sigma_acc_y^2*Ts^3/2;...
%                       sigma_acc_x^2*Ts^3/2    , 0 , sigma_acc_x^2*Ts^2, 0;...
%                       0, sigma_acc_y^2*Ts^3/2, 0, sigma_acc_y^2*Ts^2         ]; % equation (7)
%

Q=[(sigma_acc_x^2)*(Ts^4)/2, 0 , sigma_acc_x^2*Ts^3/2, 0;...
    0, (sigma_acc_y^2)*(Ts^4)/2, 0 , sigma_acc_y^2*Ts^3/2;...
    sigma_acc_x^2*Ts^3/2    , 0 , sigma_acc_x^2*Ts^2, 0;...
    0, sigma_acc_y^2*Ts^3/2, 0, sigma_acc_y^2*Ts^2         ]; % equation (7)

R=(sigma_v^2)*eye(m);


%% Simulation
%%%%%%%%%%%%%%%%%%%%%%

%% Data Generation

% r1=sqrt((3*0.6)^2+(2*0.6)^2);
% r2=sqrt((3*0.6)^2+(2*0.6)^2);
% r3=1.2;


% Z_Noised=[R1_2D_m'.^2;R2_2D_m'.^2;R3_2D_m'.^2];


%% Kalman filter

%%%% calculating the initial guess: Guo& Qiu article (UWB based
%%%% localization....)

%A=2*[x1-x2 , y1-y2 ; x1-x3 , y1-y3];
A=2*[p_anch_1_x-p_anch_2_x , p_anch_1_y-p_anch_2_y ;...
    p_anch_1_x-p_anch_3_x , p_anch_1_y-p_anch_3_y];
d1=R1_2D_m(1);
d2=R2_2D_m(1);
d3=R3_2D_m(1);

% m1=d2^2-d1^2-(x2^2-x1^2+y2^2-y1^2);
% m2=d3^2-d1^2-(x3^2-x1^2+y3^2-y1^2);

m1=d2^2-d1^2-(p_anch_2_x^2-p_anch_1_x^2+p_anch_2_y^2-p_anch_2_y^2);
m2=d3^2-d1^2-(p_anch_3_x^2-p_anch_1_x^2+p_anch_3_y^2-p_anch_1_y^2);

p_tag_init=A\[m1;m2];

%%% Initialization
X_est_apr=zeros(n,N);

%x_est_postr=[1.2;0;0;0];
x_est_postr=[p_tag_init;0;0]; %we assume the initial velocities are zero
X_est_postr=[x_est_postr,zeros(n,N-1)];


P_est_postr=eye(n); % estimation error Covariance matrix

u=[acc_x;acc_y]; % input vector= acceleration

F=[1 0 Ts 0; 0 1 0 Ts;0 0 1 0;0 0 0 1]; % state equation matrix

%%% State functions

p_anch_x=[p_anch_1_x;p_anch_2_x;p_anch_3_x];
p_anch_y=[p_anch_1_y;p_anch_2_y;p_anch_3_y];

p_anch=[p_anch_x,p_anch_y];

p_anch_1=[p_anch_1_x;p_anch_1_y];
p_anch_2=[p_anch_2_x;p_anch_2_y];
p_anch_3=[p_anch_3_x;p_anch_3_y];


%%% State equation
u=[acc_x;acc_y]; % input vector= acceleration

F=[1 0 Ts 0; 0 1 0 Ts;0 0 1 0;0 0 0 1]; % state equation matrix

f=@(x,w) F*x+w; % we model acceleration as noise, hence we do not need the input vector;  f=@(x,u,w) F*x+G*u+w;

%%% Measurement equation
h=@(p,p_anch,v) (p-p_anch)'*(p-p_anch)+v; % p and p_anch are 2x1 column vectors


%%% Jacobians

A=F; % current state wrt to former state
W=eye(n); % current state wrt current state noise

H=2*[x_est_postr(1)*ones(m,1)-p_anch_x,... % current measurement wrt current state
    x_est_postr(2)*ones(m,1)-p_anch_y,...
    zeros(m,1),zeros(m,1)];

V=eye(m); % current measurment wrt current measurement noise

for kk=2:N
    
    %%% 1) Time update equations
    
    x_est_apr=f(x_est_postr,0);
    
    P_est_apr=A*P_est_postr*A'+W*Q*W';
    
    %%% 2) Measurement correction equations
    
    H=2*[x_est_postr(1)*ones(m,1)-p_anch_x,...
        x_est_postr(2)*ones(m,1)-p_anch_y,...
        zeros(m,1),zeros(m,1)];
    
    % K=P_est_apr*H'*inv(H*P_est_apr*H'+V*R*V);
    K=P_est_apr*H'/(H*P_est_apr*H'+V*R*V);
    
    z=Z_Noised(:,kk);
    h_vec=[h(x_est_apr(1:2),p_anch_1,0);h(x_est_apr(1:2),p_anch_2,0);h(x_est_apr(1:2),p_anch_3,0)];
    
    x_est_postr=x_est_apr+K*(z-h_vec);
    
    P_est_postr=(eye(n)-K*H)*P_est_apr;
    
    X_est_apr(:,kk)= x_est_apr;
    X_est_postr(:,kk)= x_est_postr;
    
end

%% Display


figure
set(gcf,'windowstyle','docked')
%plot(X_est_apr(1,:),X_est_apr(2,:))
plot(X_est_postr(1,:),X_est_postr(2,:))
title('Tag Positioning estimation (Trilateration)')
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








