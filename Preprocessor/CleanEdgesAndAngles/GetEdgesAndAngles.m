function [e,a]=GetEdgesAndAngles(x,y)
n=size(x,1);
e=[];
a=[];
for j=1:n
    e(j,1)=sqrt((x(j)-x(mod1(j+1,n)))^2+(y(j)-y(mod1(j+1,n)))^2);
    a(j,1)=angle((x(j)-x(mod1(j-1,n)))+sqrt(-1)*(y(j)-y(mod1(j-1,n))))-angle((x(mod1(j+1,n))-x(j))+sqrt(-1)*(y(mod1(j+1,n))-y(j)));
end
a(find(a<0))=a(find(a<0))+2*pi;
a=pi-a;