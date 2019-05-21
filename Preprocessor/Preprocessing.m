function [NUMBER_BODIES,NUMBER_RIGIDS,BODIES_STATIC,BODIES_DYNAMIC,BODIES_CURRENT,MATERIALS,CONTACT_LAWS,TIME_INI,TIME_STEP,TIME_END,TIME,SAVE_PERIOD,PRINT_PERIOD,CONTACT_UPDATING_PERIOD,PERIODIC_BOUNDARIES,SCHEME,TARGET_ERROR,CONTROL_PARAMETER,ACCEPTED_RATIO,GRAVITY,MONITORINGS,DESACTIVATIONS,SPIES]=Preprocessing(Contours,Distributions,Interpolations,Integrations,Detections,Mesh_Ratios,Materials,Contact_Laws,Bodies_Materials,Imposed_Pressures,Imposed_Velocities,Initial_Velocities,Periodic_Boundaries,Gravity,Time_Stepping_Parameters,Save_Periods,Contact_Updating_Period,Scheme,Scheme_Parameters,Activate_Plot,Monitorings,Desactivations,Spies,Status,Alid)
disp(' ')
disp('Preprocessing')
disp(' ')
addpath('C:\MESHFREE\PREPRO\Distmesh')
%NUMBER_BODIES=size(Contours,1);
NUMBER_BODIES=0;
NUMBER_RIGIDS=0;
for Body=1:size(Contours,1);
    tic
    Distribution_Type=Distributions{Body,1};
    Shape_Type=Interpolations{Body,1};
    Body_Support_Domain=Interpolations{Body,2};
    if strcmp(Distribution_Type,'Structured') | strcmp(Distribution_Type,'Unstructured') | strcmp(Distribution_Type,'Manual')
        NUMBER_BODIES=NUMBER_BODIES+1;
        CLASS='Body';
        disp(['   Body ',int2str(Body),'/',int2str(NUMBER_BODIES)])
    elseif strcmp(Distribution_Type,'Rigid')
        NUMBER_RIGIDS=NUMBER_RIGIDS+1;
        CLASS='Rigid';
        disp(['   Rigid ',int2str(Body),'/',int2str(NUMBER_RIGIDS)])
    end
    NODAL_DISTANCE=Distributions{Body,2};
    NUMBER_REGIONS=Distributions{Body,3};
    MATERIAL=Bodies_Materials{Body,1};
    Mesh_Ratio=Mesh_Ratios(Body,:);
    if strcmp(Contours{Body,1}{1,1},'Closed') & size(Contours{Body,1},1)==1
        % Bodies with closed contours
        disp('       Positionning Nodes')
        Contour=Contours{Body,1}{1,2};
        Interpolant=Contours{Body,1}{1,3};
        if strcmp(Distribution_Type,'Structured')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetStructuredNodesForClosedContour(Contour,Interpolant,NODAL_DISTANCE);
        elseif strcmp(Distribution_Type,'Unstructured')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetUnstructuredNodesForClosedContour(Contour,Interpolant,NODAL_DISTANCE,Mesh_Ratio,Activate_Plot);
        elseif strcmp(Distribution_Type,'Manual')
            INITIAL_POSITIONS=Distributions{Body,4};
            TRIANGULATION=Distributions{Body,5};
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetManualBorderForClosedContours(Contour,Interpolant,INITIAL_POSITIONS,TRIANGULATION,NODAL_DISTANCE);
        elseif strcmp(Distribution_Type,'Rigid')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetRigidNodesForClosedContour(Contour,Interpolant,NODAL_DISTANCE,Activate_Plot);
            %[INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetUnstructuredNodesForClosedContour(Contour,Interpolant,NODAL_DISTANCE,Mesh_Ratio,Activate_Plot);
        end
        NUMBER_NODES=size(INITIAL_POSITIONS,1);
        NUMBER_CELLS=size(TRIANGULATION,1);
        TYPE='Simple';
        if strcmp(Distribution_Type,'Structured') | strcmp(Distribution_Type,'Unstructured') | strcmp(Distribution_Type,'Manual') 
            disp('       Precomputing Numerical Integration')
            [INFLUENCE_DOMAINS,MASSES,INVERSE_MASSES,STIFFNESS_MATRIX,DAMPING_MATRIX,BODY_FORCES,INTEGRATION,INTEGRATION_TYPE,DOF_TO_DISP]=AnalyseSimpleBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials,Body_Support_Domain,Gravity,Integrations{Body,2},Integrations{Body,1},Shape_Type);
            CENTRE_OF_MASS=[];
        elseif strcmp(Distribution_Type,'Rigid')
            [MASSES,INVERSE_MASSES,CENTRE_OF_MASS]=AnalyseRigidBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials);
            INFLUENCE_DOMAINS=[];
            STIFFNESS_MATRIX=[];
            DAMPING_MATRIX=[];
            BODY_FORCES=[];
            INTEGRATION=[];
            INTEGRATION_TYPE=[];
            DOF_TO_DISP=[];
        end
        %
        %
        %
    elseif strcmp(Contours{Body,1}{1,1},'Closed') & size(Contours{Body,1},1)==2
        % Bodies with closed contour and a hole
        disp('       Positionning Nodes')
        Contour1=Contours{Body,1}{1,2};
        Contour2=Contours{Body,1}{2,2};
        Interpolant=Contours{Body,1}{1,3};
        %if strcmp(Distribution_Type,'Structured')
        %    [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetStructuredNodesForClosedContour(Contour,Interpolant,NODAL_DISTANCE);
        %elseif strcmp(Distribution_Type,'Unstructured')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetUnstructuredNodesForHoledContour(Contour2,Contour1,Interpolant,NODAL_DISTANCE,Mesh_Ratio,Activate_Plot);
        %elseif strcmp(Distribution_Type,'Rigid')
        %    [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetRigidNodesForClosedContour(Contour,Interpolant,NODAL_DISTANCE,Activate_Plot);
        %end
        NUMBER_NODES=size(INITIAL_POSITIONS,1);
        NUMBER_CELLS=size(TRIANGULATION,1);
        TYPE='Simple';
        if strcmp(Distribution_Type,'Structured') | strcmp(Distribution_Type,'Unstructured') | strcmp(Distribution_Type,'Manual') 
            disp('       Precomputing Numerical Integration')
            [INFLUENCE_DOMAINS,MASSES,INVERSE_MASSES,STIFFNESS_MATRIX,DAMPING_MATRIX,BODY_FORCES,INTEGRATION,INTEGRATION_TYPE,DOF_TO_DISP]=AnalyseSimpleBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials,Body_Support_Domain,Gravity,Integrations{Body,2},Integrations{Body,1},Shape_Type);
            CENTRE_OF_MASS=[];
        elseif strcmp(Distribution_Type,'Rigid')
            [MASSES,INVERSE_MASSES,CENTRE_OF_MASS]=AnalyseRigidBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials);
            INFLUENCE_DOMAINS=[];
            STIFFNESS_MATRIX=[];
            DAMPING_MATRIX=[];
            BODY_FORCES=[];
            INTEGRATION=[];
            INTEGRATION_TYPE=[];
            DOF_TO_DISP=[];
        end
        %
        %
        %
    elseif strcmp(Contours{Body,1}{1,1},'Periodic') & strcmp(Contours{Body,1}{2,1},'Periodic')
        % Bodies with periodic contours
        disp('       Positionning Nodes')
        Contour1=Contours{Body,1}{1,2};
        Contour2=Contours{Body,1}{2,2};
        Interpolant1=Contours{Body,1}{1,3};
        Interpolant2=Contours{Body,1}{2,3};
        if strcmp(Distribution_Type,'Structured') | strcmp(Distribution_Type,'Unstructured')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetUnstructuredNodesForPeriodicContours(Contour1,Contour2,Interpolant1,Interpolant2,NODAL_DISTANCE,Mesh_Ratio,Periodic_Boundaries,Activate_Plot);
            %save('Temp.mat','INITIAL_POSITIONS','BORDERS','TRIANGULATION');
            %load('Temp.mat');
        elseif strcmp(Distribution_Type,'Rigid')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetRigidNodesForPeriodicContours(Contour1,Contour2,Interpolant1,Interpolant2,NODAL_DISTANCE,Mesh_Ratio,Periodic_Boundaries,Activate_Plot);
        end
        NUMBER_NODES=size(INITIAL_POSITIONS,1);
        NUMBER_CELLS=size(TRIANGULATION,1);
        TYPE='Periodic';
        if strcmp(Distribution_Type,'Structured') | strcmp(Distribution_Type,'Unstructured')
            disp('       Precomputing Numerical Integration')
            [INFLUENCE_DOMAINS,MASSES,INVERSE_MASSES,STIFFNESS_MATRIX,DAMPING_MATRIX,BODY_FORCES,INTEGRATION,INTEGRATION_TYPE,DOF_TO_DISP]=AnalysePeriodicBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,NODAL_DISTANCE,Mesh_Ratio,MATERIAL,Materials,Body_Support_Domain,Gravity,Periodic_Boundaries,Integrations{Body,2},Integrations{Body,1},Shape_Type);
            CENTRE_OF_MASS=[];
        elseif strcmp(Distribution_Type,'Rigid')
            [MASSES,INVERSE_MASSES,CENTRE_OF_MASS]=AnalyseRigidBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials);
            INFLUENCE_DOMAINS=[];
            STIFFNESS_MATRIX=[];
            DAMPING_MATRIX=[];
            BODY_FORCES=[];
            INTEGRATION=[];
            INTEGRATION_TYPE=[];
            DOF_TO_DISP=[];
        end
    elseif strcmp(Contours{Body,1}{1,1},'PeriodicPolygon') & strcmp(Contours{Body,1}{2,1},'PeriodicPolygon')
        % Bodies with periodic contours
        disp('       Positionning Nodes')
        Contour1=Contours{Body,1}{1,2};
        Contour2=Contours{Body,1}{2,2};
        Interpolant1=Contours{Body,1}{1,3};
        Interpolant2=Contours{Body,1}{2,3};
        if strcmp(Distribution_Type,'Structured') | strcmp(Distribution_Type,'Unstructured')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetUnstructuredNodesForPeriodicPolygonContours(Contour1,Contour2,Interpolant1,Interpolant2,NODAL_DISTANCE,Mesh_Ratio,Periodic_Boundaries,Activate_Plot);
        elseif strcmp(Distribution_Type,'Rigid')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetUnstructuredNodesForPeriodicPolygonContours(Contour1,Contour2,Interpolant1,Interpolant2,NODAL_DISTANCE,Mesh_Ratio,Periodic_Boundaries,Activate_Plot);
        end
        NUMBER_NODES=size(INITIAL_POSITIONS,1);
        NUMBER_CELLS=size(TRIANGULATION,1);
        TYPE='Periodic';
        if strcmp(Distribution_Type,'Structured') | strcmp(Distribution_Type,'Unstructured')
            disp('       Precomputing Numerical Integration')
            [INFLUENCE_DOMAINS,MASSES,INVERSE_MASSES,STIFFNESS_MATRIX,DAMPING_MATRIX,BODY_FORCES,INTEGRATION,INTEGRATION_TYPE,DOF_TO_DISP]=AnalysePeriodicBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,NODAL_DISTANCE,Mesh_Ratio,MATERIAL,Materials,Body_Support_Domain,Gravity,Periodic_Boundaries,Integrations{Body,2},Integrations{Body,1},Shape_Type);
            CENTRE_OF_MASS=[];
        elseif strcmp(Distribution_Type,'Rigid')
            [MASSES,INVERSE_MASSES,CENTRE_OF_MASS]=AnalyseRigidBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials);
            INFLUENCE_DOMAINS=[];
            STIFFNESS_MATRIX=[];
            DAMPING_MATRIX=[];
            BODY_FORCES=[];
            INTEGRATION=[];
            INTEGRATION_TYPE=[];
            DOF_TO_DISP=[];
        end
    elseif strcmp(Contours{Body,1}{1,1},'Polygon') & size(Contours{Body,1},1)==1
        % Bodies with polygonal closed contours
        disp('       Positionning Nodes')
        Contour=Contours{Body,1}{1,2};
        Interpolant=Contours{Body,1}{1,3};
        if strcmp(Distribution_Type,'Structured')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetStructuredNodesForPolygonContour(Contour,Interpolant,NODAL_DISTANCE);
        elseif strcmp(Distribution_Type,'Unstructured')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetUnstructuredNodesForPolygonContour(Contour,Interpolant,NODAL_DISTANCE,Mesh_Ratio,Activate_Plot);
        elseif strcmp(Distribution_Type,'Rigid')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetRigidNodesForPolygonContour(Contour,Interpolant,NODAL_DISTANCE,Activate_Plot);
            %[INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetUnstructuredNodesForPolygonContour(Contour,Interpolant,NODAL_DISTANCE,Mesh_Ratio,Activate_Plot);
        end
        NUMBER_NODES=size(INITIAL_POSITIONS,1);
        NUMBER_CELLS=size(TRIANGULATION,1);
        %TYPE='Closed';
        TYPE='Simple';
        if strcmp(Distribution_Type,'Structured') | strcmp(Distribution_Type,'Unstructured')
            disp('       Precomputing Numerical Integration')
            [INFLUENCE_DOMAINS,MASSES,INVERSE_MASSES,STIFFNESS_MATRIX,DAMPING_MATRIX,BODY_FORCES,INTEGRATION,INTEGRATION_TYPE,DOF_TO_DISP]=AnalyseSimpleBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials,Body_Support_Domain,Gravity,Integrations{Body,2},Integrations{Body,1},Shape_Type);
            CENTRE_OF_MASS=[];
        elseif strcmp(Distribution_Type,'Rigid')
            [MASSES,INVERSE_MASSES,CENTRE_OF_MASS]=AnalyseRigidBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials);
            INFLUENCE_DOMAINS=[];
            STIFFNESS_MATRIX=[];
            DAMPING_MATRIX=[];
            BODY_FORCES=[];
            INTEGRATION=[];
            INTEGRATION_TYPE=[];
            DOF_TO_DISP=[];
        end
        %
        %
        %
    elseif strcmp(Contours{Body,1}{1,1},'Polygon') & size(Contours{Body,1},1)==2
        % Bodies with polygonal closed contours and a hole
        disp('       Positionning Nodes')
        Contour1=Contours{Body,1}{1,2};
        Contour2=Contours{Body,1}{2,2};
        Interpolant=Contours{Body,1}{1,3};
        %if strcmp(Distribution_Type,'Structured')
        %    [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetStructuredNodesForPolygonContour(Contour,Interpolant,NODAL_DISTANCE);
        %elseif strcmp(Distribution_Type,'Unstructured')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetUnstructuredNodesForHoledPolygonContour(Contour2,Contour1,Interpolant,NODAL_DISTANCE,Mesh_Ratio,Activate_Plot);
        %elseif strcmp(Distribution_Type,'Rigid')
        %    [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetRigidNodesForPolygonContour(Contour,Interpolant,NODAL_DISTANCE,Activate_Plot);
        %end
        NUMBER_NODES=size(INITIAL_POSITIONS,1);
        NUMBER_CELLS=size(TRIANGULATION,1);
        %TYPE='Closed';
        TYPE='Simple';
        if strcmp(Distribution_Type,'Structured') | strcmp(Distribution_Type,'Unstructured')
            disp('       Precomputing Numerical Integration')
            [INFLUENCE_DOMAINS,MASSES,INVERSE_MASSES,STIFFNESS_MATRIX,DAMPING_MATRIX,BODY_FORCES,INTEGRATION,INTEGRATION_TYPE,DOF_TO_DISP]=AnalyseSimpleBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials,Body_Support_Domain,Gravity,Integrations{Body,2},Integrations{Body,1},Shape_Type);
            CENTRE_OF_MASS=[];
        elseif strcmp(Distribution_Type,'Rigid')
            [MASSES,INVERSE_MASSES,CENTRE_OF_MASS]=AnalyseRigidBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials);
            INFLUENCE_DOMAINS=[];
            STIFFNESS_MATRIX=[];
            DAMPING_MATRIX=[];
            BODY_FORCES=[];
            INTEGRATION=[];
            INTEGRATION_TYPE=[];
            DOF_TO_DISP=[];
        end
        %
        %
        %
    else
        % Bodies with decomposed contours
        disp('       Positionning Nodes')
        if strcmp(Distribution_Type,'Structured')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetStructuredNodesForSimpleContours(Contours{Body,1},NODAL_DISTANCE);
        elseif strcmp(Distribution_Type,'Unstructured')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetUnstructuredNodesForSimpleContours(Contours{Body,1},NODAL_DISTANCE,Activate_Plot);
        elseif strcmp(Distribution_Type,'Manual')
            INITIAL_POSITIONS=Distributions{Body,4};
            TRIANGULATION=Distributions{Body,5};
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetManualBordersForSimpleContours(Contours{Body,1},INITIAL_POSITIONS,TRIANGULATION,NODAL_DISTANCE);
        elseif strcmp(Distribution_Type,'Rigid')
            [INITIAL_POSITIONS,BORDERS,TRIANGULATION]=SetRigidNodesForSimpleContours(Contours{Body,1},NODAL_DISTANCE,Activate_Plot);
        end
        NUMBER_NODES=size(INITIAL_POSITIONS,1);
        NUMBER_CELLS=size(TRIANGULATION,1);
        TYPE='Simple';
        if strcmp(Distribution_Type,'Structured') | strcmp(Distribution_Type,'Unstructured') | strcmp(Distribution_Type,'Manual')
            disp('       Precomputing Numerical Integration')
            [INFLUENCE_DOMAINS,MASSES,INVERSE_MASSES,STIFFNESS_MATRIX,DAMPING_MATRIX,BODY_FORCES,INTEGRATION,INTEGRATION_TYPE,DOF_TO_DISP]=AnalyseSimpleBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials,Body_Support_Domain,Gravity,Integrations{Body,2},Integrations{Body,1},Shape_Type);
            CENTRE_OF_MASS=[];
        elseif strcmp(Distribution_Type,'Rigid')
            [MASSES,INVERSE_MASSES,CENTRE_OF_MASS]=AnalyseRigidBody(NUMBER_NODES,INITIAL_POSITIONS,NUMBER_CELLS,TRIANGULATION,MATERIAL,Materials);
            INFLUENCE_DOMAINS=[];
            STIFFNESS_MATRIX=[];
            DAMPING_MATRIX=[];
            BODY_FORCES=[];
            INTEGRATION=[];
            INTEGRATION_TYPE=[];
            DOF_TO_DISP=[];
        end
    end
    NUMBER_BORDERS=size(BORDERS,1);
    DIRICHLET_BC=Imposed_Velocities{Body,1};
    NEUMANN_BC=Imposed_Pressures{Body,1};
    ALID={Alid{Body,1},Alid{Body,2},Alid{Body,3}};
      
    BODIES_STATIC(Body,1).MATERIAL=MATERIAL;
    BODIES_STATIC(Body,1).CLASS=CLASS;
    BODIES_STATIC(Body,1).TYPE=TYPE;
    BODIES_STATIC(Body,1).CENTRE_OF_MASS=CENTRE_OF_MASS;
    BODIES_STATIC(Body,1).NUMBER_NODES=NUMBER_NODES;
    BODIES_STATIC(Body,1).NUMBER_BORDERS=NUMBER_BORDERS;
    BODIES_STATIC(Body,1).NUMBER_CELLS=NUMBER_CELLS;
    BODIES_STATIC(Body,1).NUMBER_REGIONS=NUMBER_REGIONS;
    BODIES_STATIC(Body,1).NODAL_DISTANCE=NODAL_DISTANCE;
    BODIES_STATIC(Body,1).INITIAL_POSITIONS=INITIAL_POSITIONS;
    BODIES_STATIC(Body,1).BORDERS=BORDERS;
    BODIES_STATIC(Body,1).TRIANGULATION=TRIANGULATION;
    BODIES_STATIC(Body,1).INFLUENCE_DOMAINS=INFLUENCE_DOMAINS;
    BODIES_STATIC(Body,1).MASSES=MASSES;
    BODIES_STATIC(Body,1).INVERSE_MASSES=INVERSE_MASSES;
    BODIES_STATIC(Body,1).STIFFNESS_MATRIX=STIFFNESS_MATRIX;
    BODIES_STATIC(Body,1).DAMPING_MATRIX=DAMPING_MATRIX;
    BODIES_STATIC(Body,1).INTEGRATION=INTEGRATION;
    BODIES_STATIC(Body,1).INTEGRATION_TYPE=INTEGRATION_TYPE;
    BODIES_STATIC(Body,1).DOF_TO_DISP=DOF_TO_DISP;
    BODIES_STATIC(Body,1).DIRICHLET_BC=DIRICHLET_BC;
    BODIES_STATIC(Body,1).NEUMANN_BC=NEUMANN_BC;
    BODIES_STATIC(Body,1).DETECTIONS=Detections(Body,:);
    BODIES_STATIC(Body,1).STATUS=Status{Body,1};
    BODIES_STATIC(Body,1).ALID=ALID;
    
    BODIES_DYNAMIC(Body,1).CURRENT_POSITIONS=INITIAL_POSITIONS;
    BODIES_DYNAMIC(Body,1).DISPLACEMENTS=zeros(NUMBER_NODES,2);
    BODIES_DYNAMIC(Body,1).DISPLACEMENTS_PARAMETERS=zeros(NUMBER_NODES,2);
    BODIES_DYNAMIC(Body,1).VELOCITIES=[Initial_Velocities{Body,1}(1)*ones(NUMBER_NODES,1),Initial_Velocities{Body,1}(2)*ones(NUMBER_NODES,1)];
    BODIES_DYNAMIC(Body,1).VELOCITIES_PARAMETERS=BODIES_DYNAMIC(Body,1).VELOCITIES;
    BODIES_DYNAMIC(Body,1).ACCELERATIONS=zeros(NUMBER_NODES,2);
    BODIES_DYNAMIC(Body,1).ACCELERATIONS_PARAMETERS=zeros(NUMBER_NODES,2);
    BODIES_DYNAMIC(Body,1).NUMBER_NEIGHBOURS=0;
    BODIES_DYNAMIC(Body,1).NEIGHBOURS=[];
    BODIES_DYNAMIC(Body,1).INTERNAL_FORCES=zeros(NUMBER_NODES,2);
    BODIES_DYNAMIC(Body,1).CONTACT_FORCES=zeros(NUMBER_NODES,2);
    BODIES_DYNAMIC(Body,1).CONTACT_FORCES=zeros(NUMBER_NODES,2);
    BODIES_DYNAMIC(Body,1).SELF_CONTACT_FORCES=zeros(NUMBER_NODES,2);
    BODIES_DYNAMIC(Body,1).BODY_FORCES=BODY_FORCES;
    BODIES_DYNAMIC(Body,1).DIRICHLET_FORCES=zeros(NUMBER_NODES,2);
    BODIES_DYNAMIC(Body,1).NEUMANN_FORCES=zeros(NUMBER_NODES,2);
    BODIES_DYNAMIC(Body,1).DAMPING_FORCES=zeros(NUMBER_NODES,2);
    BODIES_DYNAMIC(Body,1).FORCES=zeros(NUMBER_NODES,2);
    BODIES_DYNAMIC(Body,1).CONTACTS=cell(NUMBER_BORDERS,1);
    BODIES_DYNAMIC(Body,1).CONTACT_PRESSURES=cell(NUMBER_BORDERS,1);
    BODIES_DYNAMIC(Body,1).PROXIMITIES=cell(NUMBER_BORDERS,1);
    BODIES_DYNAMIC(Body,1).BC=cell(NUMBER_BORDERS,4);
    BODIES_DYNAMIC(Body,1).BC_PRESSURES=cell(NUMBER_BORDERS,1);
    BODIES_DYNAMIC(Body,1).PREVIOUS_GRADIENT=zeros(size(INTEGRATION,1),4);
    for Border=1:BODIES_STATIC(Body).NUMBER_BORDERS
        BODIES_DYNAMIC(Body).CONTACTS{Border}=zeros(BODIES_STATIC(Body).BORDERS{Border,2},20);
        BODIES_DYNAMIC(Body).PROXIMITIES{Border}=zeros(BODIES_STATIC(Body).BORDERS{Border,2},1);
        BODIES_DYNAMIC(Body).CONTACT_PRESSURES{Border}=zeros(BODIES_STATIC(Body).BORDERS{Border,2},3);
        BODIES_DYNAMIC(Body).BC_PRESSURES{Border}=zeros(BODIES_STATIC(Body).BORDERS{Border,2},3);
    end
    
    if strcmp(Scheme,'Euler') | strcmp(Scheme,'Adaptive_Euler')
        BODIES_CURRENT(Body)=BODIES_DYNAMIC(Body);
    elseif strcmp(Scheme,'BS')
        BODIES_CURRENT(Body,1).CURRENT_POSITIONS=INITIAL_POSITIONS;
        BODIES_CURRENT(Body,1).DISPLACEMENTS=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body,1).DISPLACEMENTS_PARAMETERS=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body,1).VELOCITIES=BODIES_DYNAMIC(Body,1).VELOCITIES;
        BODIES_CURRENT(Body,1).VELOCITIES_PARAMETERS=BODIES_DYNAMIC(Body,1).VELOCITIES_PARAMETERS;
        BODIES_CURRENT(Body,1).ACCELERATIONS=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body,1).ACCELERATIONS_PARAMETERS=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body,1).NUMBER_NEIGHBOURS=0;
        BODIES_CURRENT(Body,1).NEIGHBOURS=[];
        BODIES_CURRENT(Body,1).INTERNAL_FORCES=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body,1).CONTACT_FORCES=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body,1).CONTACT_FORCES=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body,1).SELF_CONTACT_FORCES=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body,1).BODY_FORCES=BODY_FORCES;
        BODIES_CURRENT(Body,1).DIRICHLET_FORCES=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body,1).NEUMANN_FORCES=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body,1).DAMPING_FORCES=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body,1).FORCES=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body,1).CONTACTS=cell(NUMBER_BORDERS,1);
        BODIES_CURRENT(Body,1).CONTACT_PRESSURES=cell(NUMBER_BORDERS,1);
        BODIES_CURRENT(Body,1).PROXIMITIES=cell(NUMBER_BORDERS,1);
        BODIES_CURRENT(Body,1).BC=cell(NUMBER_BORDERS,4);
        BODIES_CURRENT(Body,1).BC_PRESSURES=cell(NUMBER_BORDERS,1);
        BODIES_CURRENT(Body).CONTACTS=BODIES_DYNAMIC(Body).CONTACTS;
        BODIES_CURRENT(Body).PROXIMITIES=BODIES_DYNAMIC(Body).PROXIMITIES;
        BODIES_CURRENT(Body).Up1=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body).Vp1=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body).Up2=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body).Vp2=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body).Up3=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body).Vp3=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body).Up4=zeros(NUMBER_NODES,2);
        BODIES_CURRENT(Body).Vp4=zeros(NUMBER_NODES,2);
        BODIES_DYNAMIC(Body).Up1=zeros(NUMBER_NODES,2);
        BODIES_DYNAMIC(Body).Vp1=zeros(NUMBER_NODES,2);
    end
    disp(['            -> Elapsed time : ',num2str(toc),' s'])
end

MATERIALS=Materials;
CONTACT_LAWS=Contact_Laws;
TIME_INI=Time_Stepping_Parameters(1);
TIME_STEP=Time_Stepping_Parameters(2);
TIME_END=Time_Stepping_Parameters(3);
TIME=TIME_INI;
SAVE_PERIOD=Save_Periods(1);
PRINT_PERIOD=Save_Periods(2);
CONTACT_UPDATING_PERIOD=Contact_Updating_Period;
PERIODIC_BOUNDARIES=Periodic_Boundaries;
SCHEME=Scheme;
TARGET_ERROR=Scheme_Parameters(1);
CONTROL_PARAMETER=Scheme_Parameters(2);
ACCEPTED_RATIO=Scheme_Parameters(3);
GRAVITY=Gravity;
MONITORINGS=Monitorings;
DESACTIVATIONS=Desactivations;
SPIES=Spies;