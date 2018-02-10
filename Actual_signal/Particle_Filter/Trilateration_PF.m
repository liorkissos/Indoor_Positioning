
%%% based on the article "Ultra-Wideband-Based Localization for
%%% Quadcopter Navigation" by Kexin Guo & Zhirong Qiu
%%% Relevant material:
%%% 1) DecaRange ARM source files
%%% 2) DecaRangeRTLS_ARM_Source_Code_Guide.pdf
%%% 3) based on the article "Ultra-Wideband-Based Localization for Quadcopter Navigation" by Kexin Guo & Zhirong Qiu

clc
clear
close all

debug_flag=0

%% User defined parameters


Nth=300; % Resmapling threshold

N_s=500; % grid length

%%% 1) System model- Noises

% the bigger the state eqaution noise, the more we rely on the measurements
% and thus get closer to pure Trilateration.

% Measurement equation noise: standard deviation (v does not mean velocity)
sigma_v=1.5; % [m^2].
% State equation noise: standard deviation of the acceleration.  Guo&Qiu (page 5 left): maximum possible value
sigma_acc_x=1.5; % [m/sec^2].  acceleration
sigma_acc_y=1.5; % [m/sec^2].  acceleration

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


%%% TEMP
R1_2D_m=R1_2D_m(1:20);
R2_2D_m=R2_2D_m(1:20);
R3_2D_m=R3_2D_m(1:20);

N=length(R1_2D_m);



%% Particel Filter formulation and Initialization

%%% Initializations

% % 0) Calculating the initial guess: Guo& Qiu article via linear LS problem (UWB based
% % localization article.) equation (15)
% 
% A_LS=2*[p_anch_1_x-p_anch_2_x , p_anch_1_y-p_anch_2_y ;...
%     p_anch_1_x-p_anch_3_x , p_anch_1_y-p_anch_3_y];
% 
% d1=R1_2D_m(1);
% d2=R2_2D_m(1);
% d3=R3_2D_m(1);
% 
% m1=d2^2-d1^2-(p_anch_2_x^2-p_anch_1_x^2+p_anch_2_y^2-p_anch_2_y^2);
% m2=d3^2-d1^2-(p_anch_3_x^2-p_anch_1_x^2+p_anch_3_y^2-p_anch_1_y^2);
% 
% p_tag_init=A_LS\[m1;m2]; % LS solution



% 3)  Anchors position

p_anch_1=[p_anch_1_x;p_anch_1_y];
p_anch_2=[p_anch_2_x;p_anch_2_y];
p_anch_3=[p_anch_3_x;p_anch_3_y];

%%% HMM model probabilities

% 1) State equation

F=[1 0 Ts 0; 0 1 0 Ts;0 0 1 0;0 0 0 1]; % state equation matrix

%f=@(x,w) F*x+w; % we model acceleration as noise, hence we do not need the input vector;  f=@(x,u,w) F*x+G*u+w;
mu_xk=@(F,xk_1) F*xk_1;
p_xk_given_xk_1=@(xk,F,xk_1,Q) mvnpdf(xk,mu_xk(F,xk_1),Q); % need the pdf of a random VECTOR


%%% 2) Measurement equation
%h=@(p,p_anch,v) (p-p_anch)'*(p-p_anch)+v; % p and p_anch are 2x1 column vectors

mu_zk=@ (xk,p_anch_1,p_anch_2,p_anch_3)...
    [ (xk-p_anch_1)'*(xk-p_anch_1);  (xk-p_anch_2)'*(xk-p_anch_2)    ;   (xk-p_anch_3)'*(xk-p_anch_3)     ];

p_zk_given_xk=@(zk,mu_zk,Cov_zk) mvnpdf(zk,mu_zk,Cov_zk);

%%% General Gausssian random

Gen_xk=@(mu_xk,Cov_xk) mvnrnd (mu_xk,Cov_xk);


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

%% Particle filtering

X_particles=zeros(n,N_s,N);
W_particles=zeros(n,N_s,N);

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

% MAP estimator initialization
x_est_MAP=[[p_tag_init;0;0],zeros(n,N-1)]; 

% Particles of positioning are initialized to be randomly distributed around the supposed
% GT. velocity is supposed 0
xk=[normrnd(p_tag_init(1),0.2,[1,N_s]);normrnd(p_tag_init(2),0.15,[1,N_s]);zeros(2,N_s)]; %  a row contains a single particle of the state vector, and we have N_s particles

wk=(1/N_s)*ones(1,N_s); % uniform distribution. it is p(xk | xk_1) so even if xk and xk_1 are vector the conditional probability is a scalar



for k=2:N
   
    xk_1=xk;
    wk_1=wk;
    zk=Z_Noised(:,k);

    for i=1:N_s
        
        mu_xk_vec= mu_xk(F,xk_1(:,i)); % mu_xk=@(F,xk_1) F*xk_1;
        xk(:,i)=Gen_xk(mu_xk_vec,Q); %Gen_xk=@(mu_xk,Cov_xk) mvnrnd (mu_xk,Cov_xk);
        
        mu_zk_vec=mu_zk(xk(1:2,i),p_anch_1,p_anch_2,p_anch_3); %mu_zk=@ (xk,p_anch_1,p_anch_2,p_anch_3)
        wk(i)=wk_1(i)*p_zk_given_xk(zk,mu_zk_vec,R); %=@(zk,mu_zk,R) mvnpdf(zk,mu_zk,R);(zk,xk(1:2,i),p_anch_1,p_anch_2,p_anch_3,R);
        
    end
   
    if debug_flag && 0
        close all
        figure
        set(gcf,'windowstyle','docked')
        plot(wk)
        shg
    end
    
    wk=wk/sum(wk);
    
    aaa=find(wk==max(wk));
    x_est_MAP(:,k)=xk(:,aaa);
    
    Neff=(sum(wk.^2))^(-1);
    
    %%% Resampling
    if Neff<Nth
        
        if debug_flag
            [xk_sorted,I]=sort(xk(1,:)); %only x axis position demonstration
            wk_sorted=wk(I);
        end
        
        [xk]=resample_mat(xk,wk);
        wk=(1/N_s)*ones(1,N_s);
        
        if debug_flag
            close all
            figure
            set(gcf,'windowstyle','docked')
            subplot(2,1,1)
            plot(xk_sorted,wk_sorted)
            grid minor
            xlabel('xk sorted')
            ylabel('wk')
            title('PDF of xk before reampling')
            
            subplot(2,1,2)
            hist(xk(1,:))
            title('histogram of xk after reampling')
            xlabel('xk resmpled')
            shg
        end
        
        
    end
    
    X_particles(:,:,k)=xk;
    W_particles(:,:,k)=repmat(wk,[n 1]); % the repmat is not really necessary. done just to be consistent with X_particles
    
    
end

%% Kalman Iteration

% for kk=2:N
%
%     %%% 1) Time update equations
%
%     x_est_apr=f(x_est_postr,0); % apriori state vector estimator
%
%     P_est_apr=A*P_est_postr*A'+W*Q*W'; % apriori estimation error covariance matrix
%
%     %%% 2) Measurement correction equations
%
%     H=2*[x_est_postr(1)*ones(m,1)-p_anch_x,... % Jacobian of measurement wrt state calculation
%         x_est_postr(2)*ones(m,1)-p_anch_y,...
%         zeros(m,1),zeros(m,1)];
%
%     K=P_est_apr*H'/(H*P_est_apr*H'+V*R*V); % the K factor (K=P_est_apr*H'*inv(H*P_est_apr*H'+V*R*V);)
%
%     z=Z_Noised(:,kk);
%     h_vec=[h(x_est_apr(1:2),p_anch_1,0);h(x_est_apr(1:2),p_anch_2,0);h(x_est_apr(1:2),p_anch_3,0)];
%
%     x_est_postr=x_est_apr+K*(z-h_vec); % posteriori state vector estimator
%
%     P_est_postr=(eye(n)-K*H)*P_est_apr; % posteriori estimation error estimator
%
%     %%% 3) Estimators saving
%
%     X_est_apr(:,kk)= x_est_apr;
%     X_est_postr(:,kk)= x_est_postr;
%
% end

%%  Analysis

%%% MAP estimator
% [M ,I]=max(W_particles);
% Ind=sub2ind(size(W_particles),I,1:size(W_particles,2));
% x_est_MAP=(X_particles(Ind))';

%%% MMSE estimator
A=X_particles.*W_particles;

x_est_MMSE=squeeze(sum(A,2)); % summation along the row; along the particle index



%%% Errors
%Err_MAP=sqrt((x_GT-x_est_MAP)'*(x_GT-x_est_MAP));

%Err_MMSE=sqrt((x_GT-x_est_MMSE)'*(x_GT-x_est_MMSE));

%% Display

figure
set(gcf,'windowstyle','docked')
% plot(x_GT)
% hold on
%plot(x_est_MMSE(:,1),x_est_MMSE(:,2))
% hold on
plot(x_est_MAP(:,1),x_est_MAP(:,2))
grid minor
%legend(['Ground Thruth'],['MAP estimator. Error=',num2str(Err_MAP),''],['MMSE estimator.Error=',num2str(Err_MMSE),''])

shg







