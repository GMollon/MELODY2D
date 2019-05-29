function [x,y]=SortVertices(x,y)
c=mean([x,y]);
if (x(2)-x(1))*(y(1)-c(2))-(y(2)-y(1))*(x(1)-c(1))<0
    x=flipud(x);
    y=flipud(y);
end