clc
clear
close all

%% User defined parameters

N=1e3;

N_s=200;

sigma_v=1;
sigma_n=1;



%% State and Measurment functions and probabilities

f=@ (xk_1,v_k_1,k) xk_1/2+25*xk_1/(1+xk_1^2)+8*cos(1.2*k)+v_k_1; % State function

h=@(xk,n_k) xk^2/20+n_k; % measurement function

%     PF.p_xk_given_xk_1=@(xk,xk_1,k,sigma_v) normpdf(xk,xk_1/2+25*xk_1/(1+xk_1^2)+8*cos(1.2*k),sigma_v); % state probability
%     PF.Gen_prob=@(mu,sigma) normrnd(mu,sigma);
%
%     PF.N_s=N_s;
mu_k=@(xk_1,k) xk_1/2+25*xk_1/(1+xk_1^2)+8*cos(1.2*k);
%    p_xk_given_xk_1=@(xk,xk_1,k,sigma_v) normpdf(xk,xk_1/2+25*xk_1/(1+xk_1^2)+8*cos(1.2*k),sigma_v); % state probability
p_xk_given_xk_1=@(xk,xk_1,k,sigma_v) normpdf(xk,mu_k(xk_1,k),sigma_v); % state probability

Gen_prob=@(mu,sigma) normrnd(mu,sigma);


%% Signal and observation generation

x(1)=0;

n_1=sigma_n*randn;
z(1)=h(x(1),n_1);

for k=2:N
    
    v_k_1=sigma_v*randn;
    x(k)=f(x(k-1),v_k_1,k);
    
    n_k=sigma_n*randn;
    z(k)=h(x(k),n_k);
    
end


%% Particle filtering

X_particles=zeros(N_s,N);
W_particles=zeros(N_s,N);

xk=zeros(N_s,1);
wk=zeros(N_s,1);


%PF.q_params=[x_particles(:,1),x_particles(:,1),1*ones(N_s,1),sigma_v*ones(N_s,1)];

for k=2:N
    
    %         PF.x_particles_k=x_particles(:,k-1);
    %         PF.w_particles_k=w_particles(:,k-1);
    %         PF.z_k=z(k);
    %         x_particles_k=x_particles(:,k-1);
    %         w_particles_k=w_particles(:,k-1);
    
    xk_1=xk;
    z_k=z(k);
    
    
    
    for i=1:N_s
        
        mu=mu_k(xk_1(i),k); % mu_k=@(xk_1,k) xk_1/2+25*xk_1/(1+xk_1^2)+8*cos(1.2*k);
        xk(i)=Gen_prob(mu,sigma_v); % Gen_prob=@(mu,sigma) normrnd(mu,sigma);
        
        
    end
    
    X_particles(:,k)=xk;
    
    %     [x_particles(:,k),w_particles(:,k)]=Particle_Filter(PF);
    
    %   PF.q_params=[x_particles(:,k),x_particles(:,k-1),k+1,sigma_v];
    
end

