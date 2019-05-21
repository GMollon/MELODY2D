function [Cn]=Spectre(Ndesc,D0,D1,D2,D3,Decay1,D8,Decay2)
Cn=zeros(Ndesc,1);
Cn(1,1)=D0*Ndesc;
Cn(2,1)=D1*Cn(1,1);
Cn(3,1)=D2*Cn(1,1);
Cn(4,1)=D3*Cn(1,1);
A=Decay1;
B=log2(D3)-A*log2(3);
for i=5:8
    Cn(i,1)=2^(A*log2(i-1)+B)*Cn(1,1);
end
Cn(9,1)=D8*Cn(1,1);
A=Decay2;
B=log2(D8)-A*log2(8);
for i=10:Ndesc/2+1
    Cn(i,1)=2^(A*log2(i-1)+B)*Cn(1,1);
end
Cn(Ndesc/2+2:Ndesc,1)=flipud(Cn(2:Ndesc/2,1));