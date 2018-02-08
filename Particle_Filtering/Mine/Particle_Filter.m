clc
clear
close all

debug=0

%% User defined parameters

Nth=50; % Resmapling threshold

N_s=1000; % grid length


sigma_v=1; % state vector noise
sigma_n=1; % measurement vector


N=1e3;
%% State and Measurment functions and probabilities


%%% State equation
f=@ (xk_1,vk_1,k) xk_1/2+25*xk_1/(1+xk_1^2)+8*cos(1.2*k)+vk_1; % State function

mu_xk=@(xk_1,k) xk_1/2+25*xk_1/(1+xk_1^2)+8*cos(1.2*k);
p_xk_given_xk_1=@(xk,xk_1,k,sigma_v) normpdf(xk,mu_xk(xk_1,k),sigma_v); % state probability


%%% Measurements equation
h=@(xk,n_k) xk^2/20+n_k; % measurement function

mu_zk=@(xk) xk^2/20;
p_zk_given_xk=@(zk,xk,sigma_n) normpdf(zk,mu_zk(xk),sigma_n); % state probability

%%% General Gausssian random
Gen_prob=@(mu,sigma) normrnd(mu,sigma);


%% Signal and observation generation



x_GT=zeros(N,1);
x_Noised=zeros(N,1);
z=zeros(N,1);

%x(1)=0;
%n1=sigma_n*randn;
z(1)=h(x_GT(1),sigma_n*randn);

for k=2:N
    
    vk_1=sigma_v*randn;
    x_Noised(k)=f(x_GT(k-1),vk_1,k);
    
    x_GT(k)=f(x_GT(k-1),0,k);
    
    nk=sigma_n*randn;
    z(k)=h(x_GT(k),nk);
    
end


%% Particle filtering

X_particles=zeros(N_s,N);
W_particles=zeros(N_s,N);

xk=zeros(N_s,1);
wk=(1/N_s)*ones(N_s,1); % uniform distribution

for k=2:N
    
    
    xk_1=xk;
    wk_1=wk;
    zk=z(k);
    
    
    for i=1:N_s
        
        mu=mu_xk(xk_1(i),k); % mu_k=@(xk_1,k) xk_1/2+25*xk_1/(1+xk_1^2)+8*cos(1.2*k);
        xk(i)=Gen_prob(mu,sigma_v); % Gen_prob=@(mu,sigma) normrnd(mu,sigma);
        
        wk(i)=wk_1(i)*p_zk_given_xk(zk,xk(i),sigma_n);%p_zk_given_xk=@(zk,xk,sigma_n) normpdf(zk,mu_zk(xk),sigma_n); % state probability
        
    end
    
    wk=wk/sum(wk);
    
    Neff=(sum(wk.^2))^(-1);
    
    if Neff<Nth
        
        if debug
            [xk_sorted,I]=sort(xk);
            wk_sorted=wk(I); 
        end
        
        [xk]=resample(xk,wk);
        wk=(1/N_s)*ones(N_s,1);
        
        if debug
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
            hist(xk)
            title('histogram of xk after reampling')
            xlabel('xk resmpled')
            shg
        end
        
        
    end
    
    X_particles(:,k)=xk;
    W_particles(:,k)=wk;
    
    
end

%%  Analysis

%%% MAP estimator
[M ,I]=max(W_particles);
Ind=sub2ind(size(W_particles),I,1:size(W_particles,2));
x_est_MAP=(X_particles(Ind))';

%%% MMSE estimator
A=X_particles.*W_particles;

x_est_MMSE=(sum(A,1))';


%%% Errors
Err_MAP=sqrt((x_GT-x_est_MAP)'*(x_GT-x_est_MAP));

Err_MMSE=sqrt((x_GT-x_est_MMSE)'*(x_GT-x_est_MMSE));

%% Display

figure
set(gcf,'windowstyle','docked')
plot(x_GT)
hold on
plot(x_est_MAP)
hold on
plot(x_est_MMSE)
grid minor
legend(['Ground Thruth'],['MAP estimator. Error=',num2str(Err_MAP),''],['MMSE estimator.Error=',num2str(Err_MMSE),''])

shg