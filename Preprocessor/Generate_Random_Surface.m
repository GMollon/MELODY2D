
Period=1000;%10;
Nnodes=100;%200;
D2=0;%0;
D3=4;%0.5;
Decay1=-1.3;%-1.3;
D8=1.12;%0.1;
Decay2=-1.3;%-1.3;

Dx=Period/Nnodes;
X=[0:Dx:Period]';
[Cn]=Spectre(Nnodes,1,0,D2,D3,Decay1,D8,Decay2)/Nnodes;
Cn(1)=0;
%Cn(2:64)=0;
figure;loglog(Cn.^2);
Y=zeros(Nnodes+1,1);
for i=2:Nnodes/2
    Y=Y+Cn(i)*cos(2*X/Period*pi*(i-1)+2*pi*rand);
end

Y=smooth(Y,3)/4;

figure;plot(X,Y,'.-b');axis equal
RMS=sqrt(mean(Y(1:Nnodes).^2))
Amplitude=[min(Y),max(Y),max(Y)-min(Y)]

