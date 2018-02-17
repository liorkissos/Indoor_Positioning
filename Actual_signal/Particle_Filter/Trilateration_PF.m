
%%% based on the article "Ultra-Wideband-Based Localization for
%%% Quadcopter Navigation" by Kexin Guo & Zhirong Qiu
%%% Relevant material:
%%% 1) DecaRange ARM source files
%%% 2) DecaRangeRTLS_ARM_Source_Code_Guide.pdf
%%% 3) "Ultra-Wideband-Based Localization for Quadcopter Navigation" by Kexin Guo & Zhirong Qiu
%%% 4) "A Tutorial on Particle Filters for Online Nonlinear/Non-Gaussian Bayesian Tracking"
%%% 5) "Indoor Pos PF" Lior Kissos summary

clc
clear
close all

debug_flag=0

%% User defined parameters

%%% Choosing between working with the ranges read from the UWB sensors or
%%% their square
%Observation='Range'
Observation='Range Squared'

N_s=400; % grid length

Nth=20; % Resmapling threshold. should be regarded as an effective number of particles having a non zero probaility

%%% 1) System model- Noises

% the bigger the state eqaution noise, the more we rely on the measurements
% and thus get closer to pure Trilateration.

switch Observation
    case 'Range'
        % Measurement equation noise: standard deviation (v does not mean velocity)
        sigma_v=0.5; % [m^2].
        % State equation noise: standard deviation of the acceleration.  Guo&Qiu (page 5 left): maximum possible value
        sigma_acc_x=0.25; % [m/sec^2].  acceleration
        sigma_acc_y=0.25; % [m/sec^2].  acceleration
    case 'Range Squared'
        % Measurement equation noise: standard deviation (v does not mean velocity)
        sigma_v=0.8; % [m^2].
        % State equation noise: standard deviation of the acceleration.  Guo&Qiu (page 5 left): maximum possible value
        sigma_acc_x=0.1; % [m/sec^2].  acceleration
        sigma_acc_y=0.1; % [m/sec^2].  acceleration
        
end


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

%% Sanity checks

if strcmp(Observation,'Range') && N_s<500
    warning('Number of particles in Range mode needs to be higher')
end

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

switch Observation
    case 'Range'
        Z_Noised=[R1_2D_m';R2_2D_m';R3_2D_m']; % measurement vector is the square horizontal distance from the 3 anchors which is regarded as the noisy measurement
    case 'Range Squared'
        Z_Noised=[R1_2D_m'.^2;R2_2D_m'.^2;R3_2D_m'.^2]; % measurement vector is the square horizontal distance from the 3 anchors which is regarded as the noisy measurement
        
end


% Time
t_ticks=hex2dec(table2cell(Record(:,9))); % the DEBUG row see (DecaRangeRTLS_ARM_source.pdf)
t_ticks=t_ticks-t_ticks(2); % the record starts to be regular only from the 2nd sample and on
Ts=1/Fs; % the reading time interval in seconds

t=t_ticks*Ts; % the time axis in seconds

% length of recording
N=size(Record,1);



%% Particle Filter formulation and Initialization


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
switch Observation
    case 'Range'
        mu_zk=@ (xk,p_anch_1,p_anch_2,p_anch_3)...
            [sqrt( (xk-p_anch_1)'*(xk-p_anch_1)); sqrt( (xk-p_anch_2)'*(xk-p_anch_2))    ;   sqrt((xk-p_anch_3)'*(xk-p_anch_3) ) ];
    case 'Range Squared'
        mu_zk=@ (xk,p_anch_1,p_anch_2,p_anch_3)...
            [ (xk-p_anch_1)'*(xk-p_anch_1);  (xk-p_anch_2)'*(xk-p_anch_2)    ;   (xk-p_anch_3)'*(xk-p_anch_3)     ];
end

p_zk_given_xk=@(zk,mu_zk,Cov_zk) mvnpdf(zk,mu_zk,Cov_zk);

%%% 3) General Gausssian random

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

% 2) Measurement equation noise covariance matrix

R=(sigma_v^2)*eye(m);

%% Particle filtering


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

% Initializations

% Particles of positioning are initialized to be randomly distributed around the supposed
% GT. velocity is supposed 0 and normally distributed
xk=[normrnd(p_tag_init(1),0.2,[1,N_s]);normrnd(p_tag_init(2),0.15,[1,N_s]);...
    normrnd(0,0.1,[1,N_s]);normrnd(0,0.1,[1,N_s])]; %  a row contains a single particle of the state vector, and we have N_s particles

wk=(1/N_s)*ones(1,N_s); % uniform distribution. it is p(xk | xk_1) so even if xk and xk_1 are vector the conditional probability is a scalar

X_particles=zeros(n,N_s,N);
W_particles=zeros(n,N_s,N);

X_particles(:,:,1)=xk;
W_particles(:,:,1)=repmat(wk,[n,1]);


for k=2:N
    
    xk_1=xk;
    wk_1=wk;
    zk=Z_Noised(:,k);
    
    for i=1:N_s
        
        %%% 1) Particles grid drawing out of p(xk|xk_1)
        mu_xk_vec= mu_xk(F,xk_1(:,i)); % mu_xk=@(F,xk_1) F*xk_1;
        xk(:,i)=Gen_xk(mu_xk_vec,Q); %Gen_xk=@(mu_xk,Cov_xk) mvnrnd (mu_xk,Cov_xk);
        
         %%% 2) Weights calculation
        mu_zk_vec=mu_zk(xk(1:2,i),p_anch_1,p_anch_2,p_anch_3); %mu_zk=@ (xk,p_anch_1,p_anch_2,p_anch_3)
        wk(i)=wk_1(i)*p_zk_given_xk(zk,mu_zk_vec,R); %=@(zk,mu_zk,R) mvnpdf(zk,mu_zk,R);(zk,xk(1:2,i),p_anch_1,p_anch_2,p_anch_3,R);
        
    end

    %%% 3) Weights Normalization
    wk=wk/sum(wk);

    %%% 4) Resampling
    Neff=(sum(wk.^2))^(-1);
    
    if Neff<Nth
        
        if debug_flag
            [xk_sorted,I]=sort(xk(1,:)); %only x axis position demonstration
            wk_sorted=wk(I);
        end
        
        [xk]=resample_mat(xk,wk);
        wk=(1/N_s)*ones(1,N_s); %assign uniform distribution to the resampled particles 
        
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



%%  Analysis

%%% MAP estimator
w_particles=W_particles(1,:,:);
[M ,I]=max(squeeze(w_particles),[],1);
for k=1:N
    x_est_MAP(:,k)=X_particles(:,I(k),k);
end

%%% MMSE estimator
A=X_particles.*W_particles; % value by probability; x*p(x|z)

x_est_MMSE=squeeze(sum(A,2)); % summation along the row; along the particle index; E(x|z)=sum(x*p(x|z))


%% Display

% figure
% set(gcf,'windowstyle','docked')
% plot(x_est_MMSE(1,:),x_est_MMSE(2,:))
% hold on
% plot(x_est_MAP(1,:),x_est_MAP(2,:))
% grid minor
% legend(['MMSE estimator'],['MAP estimator'])

figure
set(gcf,'windowstyle','docked')
plot(x_est_MMSE(1,:),x_est_MMSE(2,:))
grid minor
title({['Tag Positioning MMSE estimation: Particle Filter '],[' # of particles=',num2str(N_s),'. Resampling Threshold=',num2str(Nth),''],...
    ['Noises:  \sigma^v=',num2str(sigma_v),'[m^2]. \sigma^{acc}=',num2str(sigma_acc_x),'[m/sec^2]  ']})


shg







