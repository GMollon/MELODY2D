function Int=SegmentsIntersection(X1,Y1,X2,Y2,X3,Y3,X4,Y4)
Int=0;
%
%if max([X1,X2])<min([X3,X4])
%    return
%elseif max([Y1,Y2])<min([Y3,Y4])
%    return
%elseif max([X3,X4])<min([X1,X2])
%    return
%elseif max([Y3,Y4])<min([Y1,Y2])
%    return
%end
%
if norm([X1,Y1]+[X2,Y2]-[X3,Y3]-[X4,Y4])>norm([X2,Y2]-[X1,Y1])+norm([X4,Y4]-[X3,Y3])
    return
end
%
if ((Y2-Y1)*(X3-X1)+(X1-X2)*(Y3-Y1))*((Y2-Y1)*(X4-X1)+(X1-X2)*(Y4-Y1))>0
    return
end
if ((Y4-Y3)*(X1-X3)+(X3-X4)*(Y1-Y3))*((Y4-Y3)*(X2-X3)+(X3-X4)*(Y2-Y3))>0
    return
end
Int=1;