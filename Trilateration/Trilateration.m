clc
clear
close all


%%

h_tag=1.43;

h1=2.28;
h2=2.28;
h3=2.28;

w1=[0;1.8];
w2=[3;1.2];
w3=[0;0] ;
%
r1=sqrt((3*0.6)^2+(2*0.6)^2);
r2=sqrt((3*0.6)^2+(2*0.6)^2);
r3=1.2;
%
Anch1.w=w1;
Anch2.w=w2;
Anch3.w=w3;
%
% Anch1.r=r1;
% Anch2.r=r2;
% Anch3.r=r3;

%% Patse the file

Record=readtable('col.txt');

R1_3D_mm=hex2dec(table2cell(Record(:,3)));
R2_3D_mm=hex2dec(table2cell(Record(:,4)));
R3_3D_mm=hex2dec(table2cell(Record(:,5)));


R1_2D_mm=sqrt(R1_3D_mm.^2-((h1-h_tag)*1e3).^2);
R2_2D_mm=sqrt(R2_3D_mm.^2-((h2-h_tag)*1e3).^2);
R3_2D_mm=sqrt(R3_3D_mm.^2-((h3-h_tag)*1e3).^2);

R1_2D_m=R1_2D_mm/1e3-0.24;
R2_2D_m=R2_2D_mm/1e3-0.24;
R3_2D_m=R3_2D_mm/1e3-0.12;

N=size(Record,1);





%% Solve the set of equations

W_trace=zeros(N,2);
W_trace(1,:)=[0.6,0.6];

%options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
options=optimoptions('fsolve','Algorithm','levenberg-marquardt');

for kk=2:N
    Anch1.r=R1_2D_m(kk);
    Anch2.r=R2_2D_m(kk);
    Anch3.r=R3_2D_m(kk);
    
    
    fun = @(w_tag) Trilatertion_equations(w_tag,Anch1,Anch2,Anch3);
    
         w_tag_0 = [0.6;0.6];
    %w_tag_0=(W_trace(kk-1,:))';
    %   w_tag = fsolve(fun,w_tag_0);
    
    w_tag = fsolve(fun,w_tag_0,options);
    
    
    W_trace(kk,:)=w_tag;
    
end

figure
set(gcf,'windowstyle','docked')
plot(W_trace(:,1),W_trace(:,2))
grid minor
xlabel('x[m]')
ylabel('y[m]')


