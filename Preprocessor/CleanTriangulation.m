function [Cleaned_Tri,Cleaned_Pos]=CleanTriangulation(Tri,Pos)

epsilon=1e-8;

Dx=max(Pos(:,1))-min(Pos(:,1));
Dy=max(Pos(:,2))-min(Pos(:,2));
D=sqrt(2*Dx*Dy/size(Pos,1));
Neighbours=Proximity2D(Pos(:,1),Pos(:,2),D,D);
Done=zeros(size(Pos,1),1);
for i=1:size(Pos,1)
    if Done(i,1)>1
        continue
    else
        Done(i,1)=i;
    end
    for jj=1:size(Neighbours{i,1},1)
        j=Neighbours{i,1}(jj,1);
        if j<=i | Done(j,1)>1
            continue
        end
        if norm([Pos(i,1)-Pos(j,1),Pos(i,2)-Pos(j,2)])<D*epsilon
            Done(j,1)=i;
        end
    end
end
n=1;
for i=1:size(Done,1)
    if Done(i,1)==i
        Done(i,2)=n;
        n=n+1;
    end
end

Cleaned_Tri=zeros(size(Tri,1),3);
for i=1:size(Tri,1)
    Cleaned_Tri(i,1)=Done(Done(Tri(i,1),1),2);
    Cleaned_Tri(i,2)=Done(Done(Tri(i,2),1),2);
    Cleaned_Tri(i,3)=Done(Done(Tri(i,3),1),2);
end

Cleaned_Pos=zeros(max(Done(:,2)),2);
for i=1:size(Done,1)
    if Done(i,2)>0
        Cleaned_Pos(Done(i,2),:)=Pos(i,:);
    end
end


