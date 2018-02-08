clc
clear
close all

%% User defined parameters

N=1e3;

N_s=400;

sigma_v=1;
sigma_n=1;



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
    
    X_particles(:,k)=xk;
    W_particles(:,k)=wk;
    
    %     [x_particles(:,k),w_particles(:,k)]=Particle_Filter(PF);
    
    %   PF.q_params=[x_particles(:,k),x_particles(:,k-1),k+1,sigma_v];
    
end

%%  Display

[M ,I]=max(W_particles);

Ind=sub2ind(size(W_particles),I,1:size(W_particles,2));

x_est_MAP=X_particles(Ind);

figure
set(gcf,'windowstyle','docked')
plot(x_GT)
hold on
plot(x_est_MAP)
grid minor
legend('Ground Thruth','Estimation')

shg