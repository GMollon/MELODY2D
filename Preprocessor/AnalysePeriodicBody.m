function [INFLUENCE_DOMAINS,MASSES,INVERSE_MASSES,STIFFNESS_MATRIX,DAMPING_MATRIX,BODY_FORCES,INTEGRATION,INTEGRATION_TYPE,DOF_TO_DISP]=AnalysePeriodicBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,NODAL_DISTANCE,Mesh_Ratio,MATERIAL,Materials,Support_Neighbours,Gravity,Periodic_Boundaries,Number_Gauss_Points,Integration_Type,Shape_Type)

% Initialize
INFLUENCE_DOMAINS=zeros(NUMBER_NODES,1);
MASSES=zeros(NUMBER_NODES,2);
INVERSE_MASSES=zeros(NUMBER_NODES,2);
STIFFNESS_MATRIX=cell(2*NUMBER_NODES,1);
DAMPING_MATRIX=cell(2*NUMBER_NODES,1);
BODY_FORCES=zeros(NUMBER_NODES,2);
INTEGRATION=cell(Number_Gauss_Points*NUMBER_CELLS,6);
DOF_TO_DISP=cell(NUMBER_NODES,1);
Period=Periodic_Boundaries(2)-Periodic_Boundaries(1);
param=5;
%Distance=param*NODAL_DISTANCE*max(Mesh_Ratio);
%Nright=NUMBER_NODES;
%ListRight=[1:NUMBER_NODES]';
%Nleft=NUMBER_NODES;
%ListLeft=[1:NUMBER_NODES]';
ListLeft=find(INITIAL_POSITIONS(:,1)<Periodic_Boundaries(1)+param*NODAL_DISTANCE*max(Mesh_Ratio));
ListRight=find(INITIAL_POSITIONS(:,1)>Periodic_Boundaries(2)-param*NODAL_DISTANCE*max(Mesh_Ratio));
Nright=size(ListRight,1);
Nleft=size(ListLeft,1);

% Creating the triple body
Pos=cat(1,INITIAL_POSITIONS,[INITIAL_POSITIONS(ListRight,1)-Period,INITIAL_POSITIONS(ListRight,2)],[INITIAL_POSITIONS(ListLeft,1)+Period,INITIAL_POSITIONS(ListLeft,2)]);
Tri=[];
for i=1:size(TRIANGULATION,1)
    N1=TRIANGULATION(i,1);
    N2=TRIANGULATION(i,2);
    N3=TRIANGULATION(i,3);
    X1=INITIAL_POSITIONS(N1,1);
    X2=INITIAL_POSITIONS(N2,1);
    X3=INITIAL_POSITIONS(N3,1);
    if X1<Periodic_Boundaries(1)+param*NODAL_DISTANCE*max(Mesh_Ratio) & X2<Periodic_Boundaries(1)+param*NODAL_DISTANCE*max(Mesh_Ratio) & X3>Periodic_Boundaries(2)-param*NODAL_DISTANCE*max(Mesh_Ratio)
        Tri=cat(1,Tri,[N1,NUMBER_NODES+find(ListRight==N3),N2],[Nright+NUMBER_NODES+find(ListLeft==N1),N3,Nright+NUMBER_NODES+find(ListLeft==N2)]);
    elseif X1<Periodic_Boundaries(1)+param*NODAL_DISTANCE*max(Mesh_Ratio) & X3<Periodic_Boundaries(1)+param*NODAL_DISTANCE*max(Mesh_Ratio) & X2>Periodic_Boundaries(2)-param*NODAL_DISTANCE*max(Mesh_Ratio)
        Tri=cat(1,Tri,[N1,N3,NUMBER_NODES+find(ListRight==N2)],[Nright+NUMBER_NODES+find(ListLeft==N1),Nright+NUMBER_NODES+find(ListLeft==N3),N2]);
    elseif X3<Periodic_Boundaries(1)+param*NODAL_DISTANCE*max(Mesh_Ratio) & X2<Periodic_Boundaries(1)+param*NODAL_DISTANCE*max(Mesh_Ratio) & X1>Periodic_Boundaries(2)-param*NODAL_DISTANCE*max(Mesh_Ratio)
        Tri=cat(1,Tri,[NUMBER_NODES+find(ListRight==N1),N3,N2],[N1,Nright+NUMBER_NODES+find(ListLeft==N3),Nright+NUMBER_NODES+find(ListLeft==N2)]);
    elseif X1<Periodic_Boundaries(1)+param*NODAL_DISTANCE*max(Mesh_Ratio) & X2>Periodic_Boundaries(2)-param*NODAL_DISTANCE*max(Mesh_Ratio) & X3>Periodic_Boundaries(2)-param*NODAL_DISTANCE*max(Mesh_Ratio)
        Tri=cat(1,Tri,[N1,NUMBER_NODES+find(ListRight==N3),NUMBER_NODES+find(ListRight==N2)],[Nright+NUMBER_NODES+find(ListLeft==N1),N3,N2]);
    elseif X2<Periodic_Boundaries(1)+param*NODAL_DISTANCE*max(Mesh_Ratio) & X1>Periodic_Boundaries(2)-param*NODAL_DISTANCE*max(Mesh_Ratio) & X3>Periodic_Boundaries(2)-param*NODAL_DISTANCE*max(Mesh_Ratio)
        Tri=cat(1,Tri,[NUMBER_NODES+find(ListRight==N1),NUMBER_NODES+find(ListRight==N3),N2],[N1,N3,Nright+NUMBER_NODES+find(ListLeft==N2)]);
    elseif X3<Periodic_Boundaries(1)+param*NODAL_DISTANCE*max(Mesh_Ratio) & X2>Periodic_Boundaries(2)-param*NODAL_DISTANCE*max(Mesh_Ratio) & X1>Periodic_Boundaries(2)-param*NODAL_DISTANCE*max(Mesh_Ratio)
        Tri=cat(1,Tri,[NUMBER_NODES+find(ListRight==N1),N3,NUMBER_NODES+find(ListRight==N2)],[N1,Nright+NUMBER_NODES+find(ListLeft==N3),N2]);
    elseif max([X1,X2,X3])<Periodic_Boundaries(1)+param*NODAL_DISTANCE*max(Mesh_Ratio)
        Tri=cat(1,Tri,[N1,N2,N3],NUMBER_NODES+Nright+[find(ListLeft==N1),find(ListLeft==N2),find(ListLeft==N3)]);
    elseif min([X1,X2,X3])>Periodic_Boundaries(2)-param*NODAL_DISTANCE*max(Mesh_Ratio)
        Tri=cat(1,Tri,[N1,N2,N3],NUMBER_NODES+[find(ListRight==N1),find(ListRight==N2),find(ListRight==N3)]);
    else
        Tri=cat(1,Tri,[N1,N2,N3]);
    end
end

ListNodes=unique([Tri(:,1);Tri(:,2);Tri(:,3)]);
NewTri=zeros(size(Tri,1),3);
for i=1:size(ListNodes,1)
    for j=1:size(Tri,1)
        if Tri(j,1)==ListNodes(i);
            NewTri(j,1)=i;
        end
        if Tri(j,2)==ListNodes(i);
            NewTri(j,2)=i;
        end
        if Tri(j,3)==ListNodes(i);
            NewTri(j,3)=i;
        end
    end 
end
Tri=NewTri;
Pos=Pos(ListNodes,:);
%figure;hold on;triplot(Tri,Pos(:,1),Pos(:,2),'.-b');axis equal;pause

% Analyse the triple body
[Influence_Domains,Masses,Inverse_Masses,Stiffness_Matrix,Damping_Matrix,Body_Forces,Integration,INTEGRATION_TYPE,Dof_To_Disp]=AnalyseSimpleBody(size(Pos,1),Pos,size(Tri,1),Tri,MATERIAL,Materials,Support_Neighbours,Gravity,Number_Gauss_Points,Integration_Type,Shape_Type);

% Retain only the relevant nodes
% for i=1:NUMBER_NODES*2
%     for j=1:size(Stiffness_Matrix{i,1},2)
%         index=Stiffness_Matrix{i,1}(1,j);
%         if index/2==floor(index/2)
%             node=index/2;
%             rest=0;
%         else
%             node=(index+1)/2;
%             rest=-1;
%         end
%         if node>NUMBER_NODES+Nright
%             STIFFNESS_MATRIX{i,1}(1,j)=rest+2*ListLeft(node-NUMBER_NODES-Nright);
%         elseif node>NUMBER_NODES
%             STIFFNESS_MATRIX{i,1}(1,j)=rest+2*ListRight(node-NUMBER_NODES);
%         else
%             STIFFNESS_MATRIX{i,1}(1,j)=index;
%         end
%         STIFFNESS_MATRIX{i,1}(2,j)=Stiffness_Matrix{i,1}(2,j);
%     end
%     for j=1:size(Damping_Matrix{i,1},2)
%         index=Damping_Matrix{i,1}(1,j);
%         if index/2==floor(index/2)
%             node=index/2;
%             rest=0;
%         else
%             node=(index+1)/2;
%             rest=-1;
%         end
%         if node>NUMBER_NODES+Nright
%             DAMPING_MATRIX{i,1}(1,j)=rest+2*ListLeft(node-NUMBER_NODES-Nright);
%         elseif node>NUMBER_NODES
%             DAMPING_MATRIX{i,1}(1,j)=rest+2*ListRight(node-NUMBER_NODES);
%         else
%             DAMPING_MATRIX{i,1}(1,j)=index;
%         end
%         DAMPING_MATRIX{i,1}(2,j)=Damping_Matrix{i,1}(2,j);
%     end
% end
for indexi=1:2*NUMBER_NODES
    DAMPING_MATRIX{indexi,1}=full(Damping_Matrix{indexi,1});
    STIFFNESS_MATRIX{indexi,1}=full(Stiffness_Matrix{indexi,1});
    for j=1:size(Stiffness_Matrix{indexi,1},2)
        indexj=Stiffness_Matrix{indexi,1}(1,j);
        if indexj/2==floor(indexj/2)
            restj=0;
        else
            restj=1;
        end
        nodej=(indexj+restj)/2;
        xj=Pos(nodej,1);
        yj=Pos(nodej,2);
        if xj<Periodic_Boundaries(1)
            newnodej=find(((Pos(:,1)-xj-Period).^2+(Pos(:,2)-yj).^2).^0.5<1e-6*NODAL_DISTANCE);
%             if newnodej>NUMBER_NODES
%                 full([1,indexi,nodej,indexj,restj,xj,yj,newnodej])
%                 pause
%             end
        elseif xj>=Periodic_Boundaries(2)
            newnodej=find(((Pos(:,1)-xj+Period).^2+(Pos(:,2)-yj).^2).^0.5<1e-6*NODAL_DISTANCE);
%             if newnodej>NUMBER_NODES
%                 full([2,indexi,nodej,indexj,restj,xj,yj,newnodej])
%                 pause
%             end
        else
            continue
        end
        newindexj=2*newnodej-restj;
        DAMPING_MATRIX{indexi,1}(1,j)=newindexj;
        STIFFNESS_MATRIX{indexi,1}(1,j)=newindexj;
    end
end
for i=1:NUMBER_NODES
    INFLUENCE_DOMAINS(i,1)=Influence_Domains(i,1);
    MASSES(i,:)=Masses(i,:);
    INVERSE_MASSES(i,:)=Inverse_Masses(i,:);
    BODY_FORCES(i,:)=Body_Forces(i,:);
    for j=1:size(Dof_To_Disp{i,1},1)
        if Pos(Dof_To_Disp{i,1}(j,1),1)<Periodic_Boundaries(1)
            DOF_TO_DISP{i,1}(j,1)=find(((Pos(:,1)-Pos(Dof_To_Disp{i,1}(j,1),1)-Period).^2+(Pos(:,2)-Pos(Dof_To_Disp{i,1}(j,1),2)).^2).^0.5<1e-6*NODAL_DISTANCE);
        elseif Pos(Dof_To_Disp{i,1}(j,1),1)>=Periodic_Boundaries(2)
            DOF_TO_DISP{i,1}(j,1)=find(((Pos(:,1)-Pos(Dof_To_Disp{i,1}(j,1),1)+Period).^2+(Pos(:,2)-Pos(Dof_To_Disp{i,1}(j,1),2)).^2).^0.5<1e-6*NODAL_DISTANCE);
        else
            DOF_TO_DISP{i,1}(j,1)=Dof_To_Disp{i,1}(j,1);
        end
        DOF_TO_DISP{i,1}(j,2:4)=Dof_To_Disp{i,1}(j,2:4);
    end
end
n=0;
for i=1:size(Integration,1)
    if Integration{i,6}(1)<Periodic_Boundaries(1) | Integration{i,6}(1)>Periodic_Boundaries(2)
        continue
    end
    n=n+1;
    for j=1:size(Integration{i,1},1)
        if Pos(Integration{i,1}(j,1),1)<Periodic_Boundaries(1)
            INTEGRATION{n,1}(j,1)=find(((Pos(:,1)-Pos(Integration{i,1}(j,1),1)-Period).^2+(Pos(:,2)-Pos(Integration{i,1}(j,1),2)).^2).^0.5<1e-6*NODAL_DISTANCE);
        elseif Pos(Integration{i,1}(j,1),1)>=Periodic_Boundaries(2)
            INTEGRATION{n,1}(j,1)=find(((Pos(:,1)-Pos(Integration{i,1}(j,1),1)+Period).^2+(Pos(:,2)-Pos(Integration{i,1}(j,1),2)).^2).^0.5<1e-6*NODAL_DISTANCE);
        else
            INTEGRATION{n,1}(j,1)=Integration{i,1}(j,1);
        end
    end
    INTEGRATION{n,2}=CreateListDOF(INTEGRATION{n,1});
    INTEGRATION{n,3}=Integration{i,3};
    INTEGRATION{n,4}=Integration{i,4};
    INTEGRATION{n,5}=Integration{i,5};
    INTEGRATION{n,6}=Integration{i,6};
end