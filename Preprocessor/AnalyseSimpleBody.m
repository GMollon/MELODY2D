function [INFLUENCE_DOMAINS,MASSES,INVERSE_MASSES,STIFFNESS_MATRIX,DAMPING_MATRIX,BODY_FORCES,INTEGRATION,INTEGRATION_TYPE,DOF_TO_DISP]=AnalyseSimpleBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials,Support_Neighbours,Gravity,Number_Gauss_Points,Integration_Type,Shape_Type)

% Initialize
INFLUENCE_DOMAINS=zeros(NUMBER_NODES,1);
MASSES=zeros(NUMBER_NODES,2);
BODY_FORCES=zeros(NUMBER_NODES,2);
Distributed_Forces_Vector=zeros(2*NUMBER_NODES,1);
DOF_TO_DISP=cell(NUMBER_NODES,1);
STIFFNESS_MATRIX=cell(2*NUMBER_NODES,1);
DAMPING_MATRIX=cell(2*NUMBER_NODES,1);
if NUMBER_NODES>100000
    FlagLargeSize=1;
else
    FlagLargeSize=0;
end
if FlagLargeSize==0
    Mass_Matrix=sparse(2*NUMBER_NODES,2*NUMBER_NODES);
    %Mass_Matrix=spalloc(2*NUMBER_NODES,2*NUMBER_NODES,20*Support_Neighbours*2*NUMBER_NODES);
    Stiffness_Matrix=sparse(2*NUMBER_NODES,2*NUMBER_NODES);
    %Stiffness_Matrix=spalloc(2*NUMBER_NODES,2*NUMBER_NODES,3*Support_Neighbours*2*NUMBER_NODES);
else
    Mass_cell=cell(2*NUMBER_NODES,1);
    Stiffness_cell=cell(2*NUMBER_NODES,1);
    for i=1:2*NUMBER_NODES
        Mass_cell{i}=zeros(20*Support_Neighbours,2);
        Stiffness_cell{i}=zeros(20*Support_Neighbours,2);
    end
    Mass_count=zeros(2*NUMBER_NODES,1);
    Stiffness_count=zeros(2*NUMBER_NODES,1);
end

% Define the connectivity of each node
CONNECTIONS=cell(NUMBER_NODES,1);
Segments=zeros(size(TRIANGULATION,1)*3,2);
for i=1:size(TRIANGULATION,1)
    Segments(3*(i-1)+1,:)=TRIANGULATION(i,1:2);
    Segments(3*(i-1)+2,:)=TRIANGULATION(i,2:3);
    Segments(3*(i-1)+3,:)=TRIANGULATION(i,[1,3]);
end
for i=1:NUMBER_NODES
    list1=Segments(find(Segments(:,1)==i),2);
    list2=Segments(find(Segments(:,2)==i),1);
    CONNECTIONS{i,1}=unique([list1;list2]);
end

% Determine the size of the influence domains
for i=1:NUMBER_NODES
    %[1,i,NUMBER_NODES]
%     FlagOut=0;
%     FlagNext=0;
%     Domain=i;
%     while FlagOut==0
%         NewDomain=[];
%         for j=Domain'
%             NewDomain=[NewDomain;CONNECTIONS{j,1}];
%         end
%         Domain=unique(NewDomain);
%         if FlagNext==1
%             FlagOut=1;
%         elseif size(Domain,1)>Support_Neighbours
%             FlagNext=1;
%         end
%     end
    Domain=[1:size(INITIAL_POSITIONS,1)]';
    Distances=sortrows([(INITIAL_POSITIONS(Domain,1)-INITIAL_POSITIONS(i,1)).^2+(INITIAL_POSITIONS(Domain,2)-INITIAL_POSITIONS(i,2)).^2,Domain],[1,2]);
    %
    d=sqrt(Distances(Support_Neighbours+1,1));
    INFLUENCE_DOMAINS(i,1)=d*(1-1e-8);
    %
    %INFLUENCE_DOMAINS(i,1)=1.1*sqrt(Distances(Support_Neighbours+1,1));
end



% Distance_Vector=zeros(NUMBER_NODES,1);
% for i=1:NUMBER_NODES
%     [1,i,NUMBER_NODES]
%     x1=INITIAL_POSITIONS(i,1);
%     y1=INITIAL_POSITIONS(i,2);
%     %for j=1:NUMBER_NODES
%     %    if j==i
%     %        Distance_Vector(j)=10^10;
%     %    else
%     %        x2=INITIAL_POSITIONS(j,1);
%     %        y2=INITIAL_POSITIONS(j,2);
%     %        Distance_Vector(j)=norm([x2-x1,y2-y1]);
%     %    end
%     %end
%     %
%     Distance_Vector=((INITIAL_POSITIONS(:,1)-INITIAL_POSITIONS(i,1)).^2+(INITIAL_POSITIONS(:,2)-INITIAL_POSITIONS(i,2)).^2);
%     Distance_Vector(i)=10^20;
%     for nn=5:100
%         lis=find(Distance_Vector<nn^2*NODAL_DISTANCE^2);
%         if size(lis,1)>Support_Neighbours
%             Distance_Vector=Distance_Vector(lis);
%             break
%         end
%     end
%     %
%     Sorted_Distance_Vector=sort(Distance_Vector);
%     %d1=Sorted_Distance_Vector(Support_Neighbours);
%     %
%     d1=sqrt(Sorted_Distance_Vector(Support_Neighbours));
%     %
%     INFLUENCE_DOMAINS(i,1)=d1-1e-8*NODAL_DISTANCE;
% end
% %INFLUENCE_DOMAINS=0.15*ones(size(INFLUENCE_DOMAINS,1),1);

% Compute integration weights and related quantities
for i=1:size(Materials,1)
    if strcmp(Materials{i,1},MATERIAL)
        Material_Number=i;
        break
    end
end
Material_Type=Materials{Material_Number,2};
Material_Parameters=Materials{Material_Number,3};
Rho=Material_Parameters(1,1);
Alpha=Material_Parameters(1,2);
Beta=Material_Parameters(1,3);
INTEGRATION=cell(1,6);
Number_Integration_Points=0;
%T=zeros(1,5);
if strcmp(Integration_Type,'Gauss')
    INTEGRATION_TYPE='Gauss';
    for i=1:NUMBER_CELLS
        %[2,i,NUMBER_CELLS]
        X1=INITIAL_POSITIONS(TRIANGULATION(i,1),1);
        Y1=INITIAL_POSITIONS(TRIANGULATION(i,1),2);
        X2=INITIAL_POSITIONS(TRIANGULATION(i,2),1);
        Y2=INITIAL_POSITIONS(TRIANGULATION(i,2),2);
        X3=INITIAL_POSITIONS(TRIANGULATION(i,3),1);
        Y3=INITIAL_POSITIONS(TRIANGULATION(i,3),2);
        [GaussPoints,GaussWeights,Jacobian]=TriangleGaussPoints(X1,Y1,X2,Y2,X3,Y3,Number_Gauss_Points);
        for l=1:Number_Gauss_Points
            Poi_Position=[GaussPoints(l,1),GaussPoints(l,2)];
            Weight=GaussWeights(l);
%            tic
            if strcmp(Shape_Type,'MLS')
                List_Domain_Fnodes=SupportDomain(Poi_Position,INITIAL_POSITIONS,INFLUENCE_DOMAINS,Support_Neighbours,CONNECTIONS,TRIANGULATION(i,:)');
%                t1=toc;
                List_DOF=CreateListDOF(List_Domain_Fnodes);
                Shape_Functions=MLS2DWithDerivatives(Poi_Position,INITIAL_POSITIONS,List_Domain_Fnodes,INFLUENCE_DOMAINS,3);
%                t2=toc;
            elseif strcmp(Shape_Type,'T3')
                List_Domain_Fnodes=TRIANGULATION(i,1:3)';
%                t1=toc;
                List_DOF=CreateListDOF(List_Domain_Fnodes);
                Shape_Functions=T3FEM(Poi_Position,INITIAL_POSITIONS,List_Domain_Fnodes);
%                t2=toc;
            end
            Number_Integration_Points=Number_Integration_Points+1;
            INTEGRATION{Number_Integration_Points,1}=List_Domain_Fnodes;
            INTEGRATION{Number_Integration_Points,2}=List_DOF;
            INTEGRATION{Number_Integration_Points,3}=Shape_Functions;
            INTEGRATION{Number_Integration_Points,4}=Weight;
            INTEGRATION{Number_Integration_Points,5}=Jacobian;
            INTEGRATION{Number_Integration_Points,6}=Poi_Position;
            [PK,Cauchy,Elasticity_Matrix]=ApplyConstitutiveModel([1,0;0,1],Material_Type,Material_Parameters);
%            t3=toc;
            Local_Stiffness_Matrix=ComputeGaussPointStiffnessMatrix(Shape_Functions,Weight,Jacobian,Elasticity_Matrix);
            Local_Distributed_Forces=ComputeGaussPointDistributedForces(Shape_Functions,Weight,Jacobian,Gravity,Rho);
            Local_Mass_Matrix=ComputeGaussPointMassMatrix(Shape_Functions,Weight,Jacobian,Rho);
%            t4=toc;
            Distributed_Forces_Vector(List_DOF,1)=Distributed_Forces_Vector(List_DOF,1)+Local_Distributed_Forces;
            if FlagLargeSize==0
                Stiffness_Matrix(List_DOF,List_DOF)=Stiffness_Matrix(List_DOF,List_DOF)+Local_Stiffness_Matrix;
                Mass_Matrix(List_DOF,List_DOF)=Mass_Matrix(List_DOF,List_DOF)+Local_Mass_Matrix;
            else
                Ndof=size(List_DOF,1);
                for j=1:Ndof
                    jj=List_DOF(j);
                    if Mass_count(jj)+Ndof>20*Support_Neighbours
                        a=sortrows(Mass_cell{jj}(1:Mass_count(jj),:),[1,2]);
                        b=a(1,:);
                        for k=2:size(a,1)
                            if a(k,1)==b(end,1)
                                b(end,2)=b(end,2)+a(k,2);
                            else
                                b(end+1,1)=a(k,1);
                                b(end,2)=a(k,2);
                            end
                        end
                        Mass_count(jj)=size(b,1);
                        Mass_cell{jj}=[b;zeros(20*Support_Neighbours-size(b,1),2)];
                    end
                    Mass_cell{jj}(Mass_count(jj)+1:Mass_count(jj)+Ndof,:)=[List_DOF,Local_Mass_Matrix(:,j)];
                    Mass_count(jj)=Mass_count(jj)+Ndof;
                    if Stiffness_count(jj)+Ndof>20*Support_Neighbours
                        a=sortrows(Stiffness_cell{jj}(1:Stiffness_count(jj),:),[1,2]);
                        b=a(1,:);
                        for k=2:size(a,1)
                            if a(k,1)==b(end,1)
                                b(end,2)=b(end,2)+a(k,2);
                            else
                                b(end+1,1)=a(k,1);
                                b(end,2)=a(k,2);
                            end
                        end
                        Stiffness_cell{jj}=[b;zeros(20*Support_Neighbours-size(b,1),2)];
                        Stiffness_count(jj)=size(b,1);
                    end
                    Stiffness_cell{jj}(Stiffness_count(jj)+1:Stiffness_count(jj)+Ndof,:)=[List_DOF,Local_Stiffness_Matrix(:,j)];
                    Stiffness_count(jj)=Stiffness_count(jj)+Ndof;
                end
            end
%            t5=toc;
%            T=T+[t1,t2-t1,t3-t2,t4-t3,t5-t4];
        end
    end
elseif strcmp(Integration_Type,'Nodal')
    INTEGRATION_TYPE='Nodal';
    Integration_Cells=UpdatedDelaunay(TRIANGULATION,INITIAL_POSITIONS(:,1),INITIAL_POSITIONS(:,2));
    %
    figure;hold on;plot(INITIAL_POSITIONS(:,1),INITIAL_POSITIONS(:,2),'.b')
    for i=1:size(Integration_Cells,1)
        plot(Integration_Cells{i,1}(:,1),Integration_Cells{i,1}(:,2),'-r')
    end
    axis equal
    %
    for i=1:NUMBER_NODES
        Poi_Position=INITIAL_POSITIONS(i,:);
        Cell=Integration_Cells{i,1};
        A=polyarea(Cell(:,1),Cell(:,2));
        WeightX=zeros(NUMBER_NODES,1);
        WeightY=zeros(NUMBER_NODES,1);
        for j=1:size(Cell,1)-1
            Poi_Position1=Cell(j,:);
            Poi_Position2=Cell(j+1,:);
            Length=norm(Poi_Position2-Poi_Position1);
            Tangent=(Poi_Position2-Poi_Position1)/Length;
            Normal=[Tangent(2),-Tangent(1)];
            GaussPW=GaussCoefficients(Number_Gauss_Points);
            for k=1:size(GaussPW,2)
                x=(Poi_Position1(1)+Poi_Position2(1))/2+GaussPW(1,k)*(Poi_Position1(1)-Poi_Position2(1))/2;
                y=(Poi_Position1(2)+Poi_Position2(2))/2+GaussPW(1,k)*(Poi_Position1(2)-Poi_Position2(2))/2;
                List_Domain_Fnodes=SupportDomain([x,y],INITIAL_POSITIONS,INFLUENCE_DOMAINS);
                Shape_Functions=MLS2D([x,y],INITIAL_POSITIONS,List_Domain_Fnodes,INFLUENCE_DOMAINS);
                Shapes=zeros(size(INITIAL_POSITIONS,1),1);
                Shapes(List_Domain_Fnodes,1)=Shape_Functions(:,1);
                WeightX=WeightX+Normal(1)*(Length/2)*Shapes*GaussPW(2,k);
                WeightY=WeightY+Normal(2)*(Length/2)*Shapes*GaussPW(2,k);
            end
        end
        List_Domain_Fnodes=find(WeightX~=0 & WeightY~=0);
        List_DOF=CreateListDOF(List_Domain_Fnodes);
        WeightX=WeightX(List_Domain_Fnodes)/A;
        WeightY=WeightY(List_Domain_Fnodes)/A;
        INTEGRATION{i,1}=List_Domain_Fnodes;
        INTEGRATION{i,2}=List_DOF;
        INTEGRATION{i,3}=WeightX;
        INTEGRATION{i,4}=WeightY;
        INTEGRATION{i,5}=A;
        INTEGRATION{i,6}=Poi_Position;
        [PK,Cauchy,Elasticity_Matrix]=ApplyConstitutiveModel([1,0;0,1],Material_Type,Material_Parameters);
        %
        Local_Stiffness_Matrix=PointStiffnessMatrixWithNodalIntegration(WeightX,WeightY,A,Elasticity_Matrix);
        Stiffness_Matrix(List_DOF,List_DOF)=Stiffness_Matrix(List_DOF,List_DOF)+Local_Stiffness_Matrix;
        %Distributed_Forces_Vector(2*i-1:2*i,1)=A*Rho*Gravity';
        %Mass_Matrix(2*i-1,2*i-1)=A*Rho;
        %Mass_Matrix(2*i,2*i)=A*Rho;
        List_Domain_Fnodes=SupportDomain(Poi_Position,INITIAL_POSITIONS,INFLUENCE_DOMAINS);
        List_DOF=CreateListDOF(List_Domain_Fnodes);
        Shape_Functions=MLS2DWithDerivatives(Poi_Position,INITIAL_POSITIONS,List_Domain_Fnodes,INFLUENCE_DOMAINS,3);
        %Local_Stiffness_Matrix=ComputeGaussPointStiffnessMatrix(Shape_Functions,1,A,Elasticity_Matrix);
        Local_Distributed_Forces=ComputeGaussPointDistributedForces(Shape_Functions,1,A,Gravity,Rho);
        Local_Mass_Matrix=ComputeGaussPointMassMatrix(Shape_Functions,1,A,Rho);
        %Stiffness_Matrix(List_DOF,List_DOF)=Stiffness_Matrix(List_DOF,List_DOF)+Local_Stiffness_Matrix;
        Distributed_Forces_Vector(List_DOF,1)=Distributed_Forces_Vector(List_DOF,1)+Local_Distributed_Forces;
        Mass_Matrix(List_DOF,List_DOF)=Mass_Matrix(List_DOF,List_DOF)+Local_Mass_Matrix;
        %
    end
end

if FlagLargeSize==0
    Damping_Matrix=Alpha*Stiffness_Matrix+Beta*Mass_Matrix;
end

for i=1:2*NUMBER_NODES
    %[3,i,2*NUMBER_NODES]
    if FlagLargeSize==0
        list=find(Stiffness_Matrix(i,:)~=0);
        STIFFNESS_MATRIX{i,1}=cat(1,list,Stiffness_Matrix(i,list));
        list=find(Damping_Matrix(i,:)~=0);
        DAMPING_MATRIX{i,1}=cat(1,list,Damping_Matrix(i,list));
    else
        a=sortrows(Mass_cell{i}(1:Mass_count(i),:),[1,2]);
        if isempty(a)
            continue
        end
        b=a(1,:);
        for k=2:size(a,1)
            if a(k,1)==b(end,1)
                b(end,2)=b(end,2)+a(k,2);
            else
                b(end+1,1)=a(k,1);
                b(end,2)=a(k,2);
            end
        end
        MASS_MATRIX{i,1}=b';
        a=sortrows(Stiffness_cell{i}(1:Stiffness_count(i),:),[1,2]);
        b=a(1,:);
        for k=2:size(a,1)
            if a(k,1)==b(end,1)
                b(end,2)=b(end,2)+a(k,2);
            else
                b(end+1,1)=a(k,1);
                b(end,2)=a(k,2);
            end
        end
        STIFFNESS_MATRIX{i,1}=b';
        DAMPING_MATRIX{i,1}=STIFFNESS_MATRIX{i,1}(1,:);
        DAMPING_MATRIX{i,1}(2,:)=Alpha*STIFFNESS_MATRIX{i,1}(2,:)+Beta*MASS_MATRIX{i,1}(2,:);
    end
end
for i=1:NUMBER_NODES
    %[4,i,NUMBER_NODES]
    if FlagLargeSize==0
        MASSES(i,1)=sum(Mass_Matrix(:,2*i-1));
        MASSES(i,2)=sum(Mass_Matrix(:,2*i));
    else
        if isempty(MASS_MATRIX{2*i-1,1}(2,:))
            continue
        end
        MASSES(i,1)=sum(MASS_MATRIX{2*i-1,1}(2,:));
        MASSES(i,2)=sum(MASS_MATRIX{2*i,1}(2,:));
    end
    BODY_FORCES(i,1)=Distributed_Forces_Vector(2*i-1);
    BODY_FORCES(i,2)=Distributed_Forces_Vector(2*i);
end
INVERSE_MASSES=1./MASSES;

% Precompute the change of variable DOF to Displacements
for i=1:NUMBER_NODES
    %[5,i,NUMBER_NODES]
    Poi_Position=INITIAL_POSITIONS(i,:);
    if strcmp(Shape_Type,'MLS')
        %List_Domain_Fnodes=SupportDomain(Poi_Position,INITIAL_POSITIONS,INFLUENCE_DOMAINS);
        List_Domain_Fnodes=SupportDomain(Poi_Position,INITIAL_POSITIONS,INFLUENCE_DOMAINS,Support_Neighbours,CONNECTIONS,i);
        %List_Domain_Fnodes=NODES_DOMAINS{i,1};
        %List_DOF=CreateListDOF(List_Domain_Fnodes);
        Shape_Functions=MLS2DWithDerivatives(Poi_Position,INITIAL_POSITIONS,List_Domain_Fnodes,INFLUENCE_DOMAINS,3);
        DOF_TO_DISP{i,1}=cat(2,List_Domain_Fnodes,Shape_Functions);
    elseif strcmp(Shape_Type,'T3')
        %List_Domain_Fnodes=TRIANGULATION(i,1:3)';
        %List_DOF=CreateListDOF(List_Domain_Fnodes);
        %Shape_Functions=T3FEM(Poi_Position,INITIAL_POSITIONS,List_Domain_Fnodes);
        %DOF_TO_DISP{i,1}=cat(2,List_Domain_Fnodes,Shape_Functions);
        DOF_TO_DISP{i,1}=[i,1.,0,0];
    end
    
    
    
    %List_Domain_Fnodes=SupportDomain(Poi_Position,INITIAL_POSITIONS,INFLUENCE_DOMAINS);
    %Shape_Functions=MLS2DWithDerivatives(Poi_Position,INITIAL_POSITIONS,List_Domain_Fnodes,INFLUENCE_DOMAINS,3);
    %DOF_TO_DISP{i,1}=cat(2,List_Domain_Fnodes,Shape_Functions);
end
%T