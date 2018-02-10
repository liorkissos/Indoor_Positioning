
%%%% based on the article "Ultra-Wideband-Based Localization for
%%%% Quadcopter Navigation" by Kexin Guo & Zhirong Qiu
clc
clear
close all

%% user defined

%Model='Known Acceleration';
Model='UnKnown Acceleration';

Fs=1;

N=5e2; % length of sampled vector

% Noises
sigma_w=1; % state equation noise standard deviation. irrelevant in the case of UnKnown acceleration

sigma_v=10; % measurement equiation noise standard deviation (v does not mean velocity)

n=4; % state dimension
l=2; % input dimension
m=3; % measurement dimension

acc_x=1; % acceleration
acc_y=1;
% acc_x=0.5; % acceleration
% acc_y=0.25;

vel_y=1;  % velocity
vel_x=1;

% Anchors positioning
p_anch_1_x=0;
p_anch_1_y=0;

p_anch_2_x=10;
p_anch_2_y=0;

p_anch_3_x=5;
p_anch_3_y=10;


%% Intermediate calculations

Ts=1/Fs;

k=0:Ts:(N-1)*Ts;

switch Model
    case 'Known Acceleration'
        Q=(sigma_w^2)*eye(n);
    case 'UnKnown Acceleration'
        sigma_acc_x=10*acc_x; % page 5 left. maximum possible value
        sigma_acc_y=10*acc_y; % page 5 left. maximum possible value
 
%%%  The state equation noise covariance. the Q(1,1) needs to be sigma_acc_x^2*Ts^4/4
%%% but it does not converge unless we change to sigma_acc_x^2*Ts^4/2
        
%                 Q=[sigma_acc_x^2*Ts^4/4, 0 , sigma_acc_x^2*Ts^3/2, 0;...
%                       0, sigma_acc_y^2*Ts^4/4, 0 , sigma_acc_y^2*Ts^3/2;...
%                       sigma_acc_x^2*Ts^3/2    , 0 , sigma_acc_x^2*Ts^2, 0;...
%                       0, sigma_acc_y^2*Ts^3/2, 0, sigma_acc_y^2*Ts^2         ]; % equation (7)
        

        Q=[(sigma_acc_x^2)*(Ts^4)/2, 0 , sigma_acc_x^2*Ts^3/2, 0;...
               0, (sigma_acc_y^2)*(Ts^4)/2, 0 , sigma_acc_y^2*Ts^3/2;...
              sigma_acc_x^2*Ts^3/2    , 0 , sigma_acc_x^2*Ts^2, 0;...
              0, sigma_acc_y^2*Ts^3/2, 0, sigma_acc_y^2*Ts^2         ]; % equation (7)
end

%Q=(sigma_x^2)*eye(n);
R=(sigma_v^2)*eye(m);

p_anch_1=[p_anch_1_x;p_anch_1_y];
p_anch_2=[p_anch_2_x;p_anch_2_y];
p_anch_3=[p_anch_3_x;p_anch_3_y];

%error(' data generation: lines 64-78 need to be reviewed. algo: 116-139 needs rto be reviewed as well ')

%% Simulation
%%%%%%%%%%%%%%%%%%%%%%

%% Data Generation

x_GT=zeros(1,N);
x_dot_GT=vel_x*ones(1,N);
y_GT=zeros(1,N);
y_dot_GT=vel_y*ones(1,N);

d_GT_1=zeros(1,N);
d_GT_2=zeros(1,N);
d_GT_3=zeros(1,N);

for kk=2:N
    
    %%% State parameters
    
    x_GT(kk)=x_GT(kk-1)+x_dot_GT(kk-1)*Ts+acc_x*(Ts^2)/2;
    x_dot_GT(kk)=x_dot_GT(kk-1)+acc_x*Ts;
    
    y_GT(kk)=y_GT(kk-1)+y_dot_GT(kk-1)*Ts+acc_y*(Ts^2)/2;
    y_dot_GT(kk)=y_dot_GT(kk-1)+acc_y*Ts;
    
    
    %%% Measurement parameters
    p_GT=[x_GT(kk);y_GT(kk)];
    
    d_GT_1(kk)=(p_GT-p_anch_1)'*(p_GT-p_anch_1);
    d_GT_2(kk)=(p_GT-p_anch_2)'*(p_GT-p_anch_2);
    d_GT_3(kk)=(p_GT-p_anch_3)'*(p_GT-p_anch_3);
    
end

X_GT=[x_GT;y_GT;x_dot_GT;y_dot_GT];

switch Model
    case 'Known Acceleration'
        X_Noised=X_GT+sigma_w*randn(n,N);
    case 'UnKnown Acceleration'
        X_Noised=X_GT;
end

%X_Noised=X_GT+sigma_w*randn(n,N);

Z_GT=[d_GT_1;d_GT_2;d_GT_3];

Z_Noised=Z_GT+sigma_v*randn(m,N);


%% Kalman filter

%%%
X_est_apr=zeros(n,N);
x_est_postr=zeros(n,1);


P_est_postr=eye(n); % Covariance matrix

u=[acc_x;acc_y]; % input vector= acceleration

F=[1 0 Ts 0; 0 1 0 Ts;0 0 1 0;0 0 0 1]; % state equation matrix
G=[Ts^2/2 0;0 Ts^2/2;Ts 0;0 Ts]; % measurement eqauation matrix

%%% State functions

p_anch_x=[p_anch_1_x;p_anch_2_x;p_anch_3_x];
p_anch_y=[p_anch_1_y;p_anch_2_y;p_anch_3_y];

p_anch=[p_anch_x,p_anch_y];

p_anch_1=[p_anch_1_x;p_anch_1_y];
p_anch_2=[p_anch_2_x;p_anch_2_y];
p_anch_3=[p_anch_3_x;p_anch_3_y];


%%% State equation

switch Model
    case 'Known Acceleration'
        f=@(x,u,w) F*x+G*u+w;
    case 'UnKnown Acceleration'
        %  f=@(x,u,w) F*x+G*u+w;
        f=@(x,u,w) F*x+w;
end
%f=@(x,u,w) F*x+G*u+w;
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
    
    % x_est_apr1=F*x_est_postr
    x_est_apr=f(x_est_postr,u,0);
    %    x_est_apr=f(x_est_postr,zeros(2,1),0)
    
    X_est_apr(:,kk)= x_est_apr;
    
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


shg








