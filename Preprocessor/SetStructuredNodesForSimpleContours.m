function [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetStructuredNodesForSimpleContours(Contours,Average_Nodal_Distance)
for i=1:size(Contours,1)
    if strcmp(Contours{i,1},'Simple')==0
        i=i-1;
        break
    end
end
Ncontours=i;
Contour=[];
for n=1:Ncontours
    Contour=cat(1,Contour,Contours{n,2}(1:size(Contours{n,2},1)-1,:));
end
INITIAL_POSITIONS=Contour;
BORDERS=cell(Ncontours,3);

% Boundary Field nodes
Nstart=0;
for n=1:Ncontours
    for i=1:size(Contours{n,2},1)-1
        Nini=i;
        Nfin=i+1;
        Xini=Contours{n,2}(Nini,1);
        Yini=Contours{n,2}(Nini,2);
        Xfin=Contours{n,2}(Nfin,1);
        Yfin=Contours{n,2}(Nfin,2);
        Nnewnodes=floor(norm([Xfin-Xini,Yfin-Yini])/(0.999999*Average_Nodal_Distance))-1;
        if Nnewnodes>0
            if abs(Xini-Xfin)<1e-8*Average_Nodal_Distance
                XNewnodes=Xini*ones(Nnewnodes,1);
            else
                XNewnodes=[Xini+0.999999*(Xfin-Xini)/(Nnewnodes+1):0.999999*(Xfin-Xini)/(Nnewnodes+1):Xfin-0.999999*(Xfin-Xini)/(Nnewnodes+1)]';
            end
            if abs(Yini-Yfin)<1e-8*Average_Nodal_Distance
                YNewnodes=Yini*ones(Nnewnodes,1);
            else
                YNewnodes=[Yini+0.999999*(Yfin-Yini)/(Nnewnodes+1):0.999999*(Yfin-Yini)/(Nnewnodes+1):Yfin-0.999999*(Yfin-Yini)/(Nnewnodes+1)]';
            end
            NewNodes=[size(INITIAL_POSITIONS,1)+1:size(INITIAL_POSITIONS,1)+Nnewnodes]';
            NewSegments=cat(2,NewNodes,NewNodes+1);
            NewSegments=cat(1,[Nini+Nstart,NewSegments(1,1)],NewSegments);
            NewSegments(size(NewSegments,1),2)=Nfin+Nstart;
            INITIAL_POSITIONS=cat(1,INITIAL_POSITIONS,[XNewnodes,YNewnodes]);
            BORDERS{n,3}=cat(1,BORDERS{n,3},NewSegments);
        else
            NewSegments=[Nini,Nfin];
            BORDERS{n,3}=cat(1,BORDERS{n,3},NewSegments);
        end
    end
    Nstart=Nstart+size(Contours{n,2},1)-1;
    BORDERS{n,1}='Simple';
    BORDERS{n,2}=size(BORDERS{n,3},1)+1;
    if n==Ncontours
        BORDERS{n,3}=[BORDERS{n,3}(:,1);1];
    else
        BORDERS{n,3}=[BORDERS{n,3}(:,1);BORDERS{n,3}(size(BORDERS{n,3},1),2)];
    end
    BORDERS{n,4}=Contours{n,3};
end

% Fill with nodes and triangulate
Xmin=min(Contour(:,1));Xmax=max(Contour(:,1));
Ymin=min(Contour(:,2));Ymax=max(Contour(:,2));
Nxinter=round((Xmax-Xmin)/Average_Nodal_Distance);
Nyinter=round((Ymax-Ymin)/Average_Nodal_Distance);
Xinter=(Xmax-Xmin)/Nxinter;
Yinter=(Ymax-Ymin)/Nyinter;
for i=Xmin+Xinter:Xinter:Xmax-Xinter
    for j=Ymin+Yinter:Yinter:Ymax-Yinter
        Xcandidate=i;
        Ycandidate=j;
        if inpolygon(Xcandidate,Ycandidate,Contour(:,1),Contour(:,2))
            add=1;
            for j=1:size(BORDERS,1)
                for k=1:size(BORDERS{j,3},1)
                    xp=INITIAL_POSITIONS(BORDERS{j,3}(k,1),1);
                    yp=INITIAL_POSITIONS(BORDERS{j,3}(k,1),2);
                    if sqrt((Xcandidate-xp)^2+(Ycandidate-yp)^2)<Average_Nodal_Distance/2
                        add=0;
                        break
                    end
                end
                if add==0
                    break
                end
            end
            if add==1
                INITIAL_POSITIONS=cat(1,INITIAL_POSITIONS,[Xcandidate,Ycandidate]);
            end
        end
    end
end
TRIANGULATION=delaunay(INITIAL_POSITIONS(:,1),INITIAL_POSITIONS(:,2));