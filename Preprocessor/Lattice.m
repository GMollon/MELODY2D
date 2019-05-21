function Points=Lattice(Npoints,xmin,xmax,ymin,ymax,noise)
Points=[];
dist=sqrt((xmax-xmin)*(ymax-ymin)/Npoints);
nlines=ceil((ymax-ymin)/(dist*atan(pi/3)));
ncol=ceil((xmax-xmin)/dist);
ydist=(ymax-ymin)/nlines;
xdist=(xmax-xmin)/ncol;
for i=1:nlines
    y=ymin+ydist/2+(i-1)*ydist;
    if y>ymax
        break
    end
    if floor(i/2)==i/2
        for j=1:ncol
            x=xmin+xdist/4+(j-1)*xdist;
            if x>xmax
                break
            end
            Points=[Points;[x,y]];
        end
    else
        for j=1:ncol
            x=xmin+3*xdist/4+(j-1)*xdist;
            if x>xmax
                break
            end
            Points=[Points;[x,y]];
        end
    end
end
Points=Points+rand(size(Points,1),2)*noise*xdist;
for i=1:size(Points,1)
    if Points(i,1)<xmin
        Points(i,1)=Points(i,1)+(xmax-xmin);
    elseif Points(i,1)>xmax
        Points(i,1)=Points(i,1)-(xmax-xmin);
    end
    if Points(i,2)<ymin
        Points(i,2)=Points(i,2)+(ymax-ymin);
    elseif Points(i,2)>ymax
        Points(i,2)=Points(i,2)-(ymax-ymin);
    end
end
%plot(Points(:,1),Points(:,2),'.b');axis equal