function [ xk_p_1,wk_p_1 ] = Particle_Filter( PF )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%%

%xk=PF.x_particle_k;
%w_particle_k=PF.w_particle_k;
%z=PF.z_k;

p_xk_given_xk_1=PF.p_xk_given_xk_1;
Gen_prob=PF.Gen_prob;

xk=  PF.x_particles_k;
wk=   PF.w_particles_k;
zk=     PF.z_k;
q_params= PF.q_params;
%%%%%%%%%%%%%%%



xk_p_1=1;
wk_p_1=1;

%%

for i=1:PF.N_s
    
    xk=p_xk_given_xk_1(xk,xk);
    
%    PF.p_xk_given_xk_1=@(xk,xk_1,k,sigma_v)
    
end



end

