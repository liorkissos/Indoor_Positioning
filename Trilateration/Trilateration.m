clc
clear
close all


%% User defined

GT=[1.25,0 ; 1.25,0.6 ; 1.25,1.2 ; 1.85,1.2 ; 2.4,1.2 ; 1.85,1.2 ; 1.85,0.6 ; 1.85, 0 ; 1.25, 0];

h_tag=1.43;

h1=2.28;
h2=2.28;
h3=2.28;

x1=0;
y1=1.8;

x2=3;
y2=1.2;

x3=0;
y3=0;

w1=[x1;y1];
w2=[x2;y2];
w3=[x3;y3] ;


% % %
% r1=sqrt((3*0.6)^2+(2*0.6)^2);
% r2=sqrt((3*0.6)^2+(2*0.6)^2);
% r3=1.2;
% % %

Anch1.w=w1;
Anch2.w=w2;
Anch3.w=w3;

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


%%%% calculating the initial guess

A=2*[x1-x2 , y1-y2 ; x1-x3 , y1-y3];
d1=R1_2D_m(1);
d2=R2_2D_m(1);
d3=R3_2D_m(1);

m1=d2^2-d1^2-(x2^2-x1^2+y2^2-y1^2);
m2=d3^2-d1^2-(x3^2-x1^2+y3^2-y1^2);

w_tag_0=A\[m1;m2];

%%% Sequential Trilateration 

W_init_guess=[w_tag_0' ; zeros(N-1,2)];
W_trace= zeros(N,2);
options=optimoptions('fsolve','Algorithm','levenberg-marquardt');
%options=optimoptions('fsolve','Algorithm','levenberg-marquardt','OptimalityTolerance',20e-6,'StepTolerance',20e-6);

for kk=1:N

    Anch1.r=R1_2D_m(kk);
    Anch2.r=R2_2D_m(kk);
    Anch3.r=R3_2D_m(kk);
    
    
    fun = @(w_tag) Trilateration_equations(w_tag,Anch1,Anch2,Anch3);
    
    %w_tag_0 = [0.6;0.6];
   % w_tag_0=(W_trace(kk-1,:))';
    %   w_tag = fsolve(fun,w_tag_0);
    
    w_tag = fsolve(fun,(W_init_guess(kk,:))',options);
    
    
    W_trace(kk,:)=real(w_tag);
    W_init_guess(kk+1,:)=real(w_tag);
   
end

%% Analysis

Diff_trace=(diff(W_trace,[],1));

Diff_trace_norm=Diff_trace(:,1).^2+Diff_trace(:,2).^2;



figure
set(gcf,'windowstyle','docked')
plot(W_trace(:,1),W_trace(:,2));
hold on
graph1=plot(GT(:,1),GT(:,2));
set(graph1,'LineWidth',2)
title('Tag trace')
grid minor
xlabel('x[meter]')
ylabel('y[meter]')
legend('Estimated','Ground Truth','Location','best')

figure
set(gcf,'windowstyle','docked')
plot(Diff_trace_norm)
title('Tag trace norm of differece between consecutive samples')
xlabel('Sample')
ylabel('norm of difference')
grid minor
