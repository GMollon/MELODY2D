function [x,y,R]=Fourier_Recons(Teta,An,Bn,Xc,Yc)
Npoints=size(Teta,1);
Ndesc=size(An,1);
R=zeros(Npoints,1);

for nt=1:Npoints
    t=Teta(nt,1);
    for i=0:Ndesc/2-1
        %R(nt,1)=R(nt,1)+An(i+1,1)*cos(i*t)+Bn(i+1,1)*sin(i*t);
        if An(i+1,1)>0
            deltai=atan(Bn(i+1,1)/An(i+1,1));
        else
            deltai=atan(Bn(i+1,1)/An(i+1,1))+3.14159;
        end
        amplitudei=sqrt(An(i+1,1)^2+Bn(i+1,1)^2);
        R(nt,1)=R(nt,1)+amplitudei*cos(i*t-deltai);
    end
end
R=abs(R/Ndesc);

x=zeros(Npoints,1);
y=zeros(Npoints,1);
for i=1:Npoints
    x(i,1)=Xc+R(i,1)*cos(Teta(i,1));
    y(i,1)=Yc+R(i,1)*sin(Teta(i,1));
end
%plot(x,y,'.r');hold on;plot(x,y,'-r');axis equal
