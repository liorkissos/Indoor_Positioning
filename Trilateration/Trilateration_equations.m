function  F  = Trilateration_equations(w_tag,Anch1,Anch2,Anch3 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

r1=Anch1.r;
r2=Anch2.r;
r3=Anch3.r;

w1=Anch1.w;
w2=Anch2.w;
w3=Anch3.w;

% w1=[0;1.8];
% w2=[3;1.2];
% w3=[0;0] ;
% 
% r1=sqrt((3*0.6)^2+(2*0.6)^2);
% r2=sqrt((3*0.6)^2+(2*0.6)^2);
% r3=1.2;

F(1)=r1^2-(w1-w_tag)'*(w1-w_tag);
F(2)=r2^2-(w2-w_tag)'*(w2-w_tag);
F(3)=r3^2-(w3-w_tag)'*(w3-w_tag);

end

