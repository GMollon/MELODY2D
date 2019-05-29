function [x,y]=CleanEdgesAndAngles(x,y,MINedge,MINangle)
[x,y]=SortVertices(x,y);
[e,a]=GetEdgesAndAngles(x,y);
n=size(x,1);
    
% Merging successive small segments
flagOK=0;
while flagOK==0
    flagOK=1;
    for j=1:n
        em=e(mod1(j-1,n));
        ep=e(mod1(j,n));
        if em<MINedge & ep<MINedge
            flagOK=0;
            if j==1
                x=x(2:n);
                y=y(2:n);
            elseif j==n
                x=x(1:n-1);
                y=y(1:n-1);
            else
                x=[x(1:j-1);x(j+1:n)];
                y=[y(1:j-1);y(j+1:n)];
            end
            [e,a]=GetEdgesAndAngles(x,y);
            n=size(x,1);
            break
        end
    end
end

% Cleaning Angles
while min(a/pi*180)<MINangle
    [m,j]=min(a);
    em=e(mod1(j-1,n));
    ep=e(mod1(j,n));
    if em>MINedge & ep>MINedge
        xm=x(j)-(x(j)-x(mod1(j-1,n)))/em*MINedge;
        ym=y(j)-(y(j)-y(mod1(j-1,n)))/em*MINedge;
        xp=x(j)+(x(mod1(j+1,n))-x(j))/ep*MINedge;
        yp=y(j)+(y(mod1(j+1,n))-y(j))/ep*MINedge;
    elseif em<MINedge & ep>MINedge
        xm=x(mod1(j-1,n));
        ym=y(mod1(j-1,n));
        xp=x(j)+(x(mod1(j+1,n))-x(j))/ep*MINedge;
        yp=y(j)+(y(mod1(j+1,n))-y(j))/ep*MINedge;
    elseif em>MINedge & ep<MINedge
        xm=x(j)-(x(j)-x(mod1(j-1,n)))/em*MINedge;
        ym=y(j)-(y(j)-y(mod1(j-1,n)))/em*MINedge;
        xp=x(mod1(j+1,n));
        yp=y(mod1(j+1,n));
    end
    if j==1
        x=[xm;xp;x(2:n)];
        y=[ym;yp;y(2:n)];
    elseif j==n
        x=[x(1:n-1);xm;xp];
        y=[y(1:n-1);ym;yp];
    else
        x=[x(1:j-1);xm;xp;x(j+1:n)];
        y=[y(1:j-1);ym;yp;y(j+1:n)];
    end
    [e,a]=GetEdgesAndAngles(x,y);
    n=size(x,1);
end
   
% Cleaning Edges
while min(e)<MINedge
    [m,j]=min(e);
    em=e(mod1(j-1,n));
    ep=e(mod1(j+1,n));
    xm=x(j)-(x(j)-x(mod1(j-1,n)))/em*MINedge;
    ym=y(j)-(y(j)-y(mod1(j-1,n)))/em*MINedge;
    xp=x(mod1(j+1,n))+(x(mod1(j+2,n))-x(mod1(j+1,n)))/ep*MINedge;
    yp=y(mod1(j+1,n))+(y(mod1(j+2,n))-y(mod1(j+1,n)))/ep*MINedge;
    if j==1
        x=[xm;xp;x(3:n)];
        y=[ym;yp;y(3:n)];
    elseif j==n-1
        x=[x(1:n-2);xm;xp];
        y=[y(1:n-2);ym;yp];
    elseif j==n
        x=[xp;x(2:n-1);xm];
        y=[yp;y(2:n-1);ym];
    else
        x=[x(1:j-1);xm;xp;x(j+2:n)];
        y=[y(1:j-1);ym;yp;y(j+2:n)];
    end
    [e,a]=GetEdgesAndAngles(x,y);
    n=size(x,1);
end