clear all
Width=10000;
Height=1000;
Grainsize=10;
MinSideLength=1;

%[C,Cup,Cdown]=GenerateMicrostructurePorous('Lattice','Flat','Flat','Periodic',Width*Height/Grainsize^2,0,Width,0,Height,1,MinSideLength,70,0,[]);
%return

Points=Lattice(Width*Height/Grainsize^2,0,Width,0,Height,1);
Npoints=size(Points,1);
%Pointsxp=[Points(:,1)+(Width-0),Points(:,2)];
%Pointsxm=[Points(:,1)-(Width-0),Points(:,2)];
Pointsxm=[-Points(:,1),Points(:,2)];
Pointsxp=[Width+(Width-Points(:,1)),Points(:,2)];
%
Points=[Points;Pointsxp;Pointsxm];
Pointsym=[Points(:,1),0-(Points(:,2)-0)];
Pointsyp=[Points(:,1),Height+(Height-Points(:,2))];
Points=[Points;Pointsyp;Pointsym];
[vertices,cells]=voronoin(Points);
GrainsContours=cell(Npoints,1);
for i=1:Npoints
    GrainsContours{i,1}=[vertices(cells{i,1}',:);vertices(cells{i,1}(1),:)];
end

G=GrainsContours;
minsize=MinSideLength;
Loop=0;
Endloop=0;
while Endloop==0;
    Ncorrect=0;
    Endloop=1;
    Loop=Loop+1;
    for i=1:size(G,1)
        C=G{i,1};
        for j=1:size(C,1)-1
            d=sqrt((C(j,1)-C(j+1,1))^2+(C(j,2)-C(j+1,2))^2);
            if d<minsize
                Ncorrect=Ncorrect+1;
                Endloop=0;
                if j==1
                    Cm=C(size(C,1)-1,:);
                else
                    Cm=C(j-1,:);
                end
                if j==size(C,1)-1
                    Cp=C(2,:);
                else
                    Cp=C(j+2,:);
                end
                C(j,:)=C(j,:)+0.1*(Cm-C(j,:));
                C(j+1,:)=C(j+1,:)+0.1*(Cp-C(j+1,:));
                if j==1
                    C(size(C,1),:)=C(1,:);
                elseif j==size(C,1)-1
                    C(1,:)=C(size(C,1),:);
                end
            end
        end
        G{i,1}=C;
    end
    [Loop,Ncorrect]
end

GrainsContours=G;
for i=1:size(GrainsContours,1)
    for j=1:size(GrainsContours{i,1},1)
        if abs(GrainsContours{i,1}(j,1)<1e-12)
            GrainsContours{i,1}(j,1)=0;
        end
        if abs(GrainsContours{i,1}(j,2)<1e-12)
            GrainsContours{i,1}(j,2)=0;
        end
    end
end
return
figure;
for i=1:size(G,1)
    patch(GrainsContours{i,1}(:,1),GrainsContours{i,1}(:,2),rand(1,3),'linestyle','none')
    hold on
end
axis equal
