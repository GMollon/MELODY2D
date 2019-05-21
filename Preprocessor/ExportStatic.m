file_data=[Simulation_Name,'/STATIC_DATA.asc'];
fid_data = fopen(file_data,'w');

fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid_data,'%s\n','%%%%%%%          GENERAL DATA          %%%%%%%');
fprintf(fid_data,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid_data,'%s\n',' ');

%%% SIMULATION NAME %%%
fprintf(fid_data,'%s\n','SIMULATION_NAME');
fprintf(fid_data,'%s\n',Simulation_Name);
fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n',' ');

%%% NUMBER_BODIES %%%
fprintf(fid_data,'%s\n','NUMBER_BODIES');
fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC,1));
fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n',' ');

%%% BODIES AND RIGIDS %%%
for Body=1:size(BODIES_STATIC,1)
    if strcmp(BODIES_STATIC(Body,1).CLASS,'Body')
        fprintf(fid_data,'%s\n',' ');
        fprintf(fid_data,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid_data,'%s\n',['%%%%%%%              BODY ',int2str(Body-1),'            %%%%%%%']);
        fprintf(fid_data,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid_data,'%s\n',' ');

        %%% PROPERTIES %%%
        fprintf(fid_data,'%s\n','DEFORMABLE');
        fprintf(fid_data,'%0.16g\n',Body-1);
        fprintf(fid_data,'%s\n',' ');

        %%% NODES %%%
        fprintf(fid_data,'%s\n','NODES');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_NODES);
        fprintf(fid_data,'%s\n',' ');
        table=cat(2,[1:BODIES_STATIC(Body,1).NUMBER_NODES]'-1,BODIES_STATIC(Body,1).INITIAL_POSITIONS,BODIES_STATIC(Body,1).INFLUENCE_DOMAINS,BODIES_STATIC(Body,1).MASSES,BODIES_STATIC(Body,1).INVERSE_MASSES);
        fprintf(fid_data,'%0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g\n',table');
        fprintf(fid_data,'%s\n',' ');

        %%% BORDERS %%%
        fprintf(fid_data,'%s\n','BORDERS');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_BORDERS);
        fprintf(fid_data,'%s\n',' ');
        for Border=1:BODIES_STATIC(Body,1).NUMBER_BORDERS
            fprintf(fid_data,'%s',BODIES_STATIC(Body,1).BORDERS{Border,1});
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%s\n',BODIES_STATIC(Body,1).BORDERS{Border,4});
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).BORDERS{Border,2});
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).BORDERS{Border,5});
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).BORDERS{Border,6});
            fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).BORDERS{Border,3}-1);
            fprintf(fid_data,'%s\n',' ');
        end

        %%% INTEGRATION %%%
        fprintf(fid_data,'%s\n','INTEGRATION');
        fprintf(fid_data,'%s\n',BODIES_STATIC(Body,1).INTEGRATION_TYPE);
        if strcmp(BODIES_STATIC(Body,1).INTEGRATION_TYPE,'Gauss')
            fprintf(fid_data,'%0.16g',size(BODIES_STATIC(Body,1).INTEGRATION,1));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_REGIONS);
            fprintf(fid_data,'%s\n',' ');
            for gp=1:size(BODIES_STATIC(Body,1).INTEGRATION,1)
                fprintf(fid_data,'%s\n',['Gauss Point ',int2str(gp-1)]);
                fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).INTEGRATION{gp,6}(1));
                fprintf(fid_data,'%s',' ');
                fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).INTEGRATION{gp,6}(2));
                fprintf(fid_data,'%s',' ');
                fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).INTEGRATION{gp,4});
                fprintf(fid_data,'%s',' ');
                fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).INTEGRATION{gp,5});
                table=BODIES_STATIC(Body,1).INTEGRATION{gp,1}-1;
                for i=1:size(BODIES_STATIC(Body,1).INTEGRATION{gp,1},1)
                    table(i,2:3)=BODIES_STATIC(Body,1).INTEGRATION{gp,2}(2*i-1:2*i,1)'-1;
                end
                table=cat(2,table,BODIES_STATIC(Body,1).INTEGRATION{gp,3});
                fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).INTEGRATION{gp,1},1));
                fprintf(fid_data,'%0.16g %0.16g %0.16g %0.16g %0.16g %0.16g\n',table');
                fprintf(fid_data,'%s\n',' ');
            end
        elseif strcmp(BODIES_STATIC(Body,1).INTEGRATION_TYPE,'Nodal')
            fprintf(fid_data,'%0.16g',size(BODIES_STATIC(Body,1).INTEGRATION,1));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_REGIONS);
            fprintf(fid_data,'%s\n',' ');
            for n=1:size(BODIES_STATIC(Body,1).INTEGRATION,1)
                fprintf(fid_data,'%s\n',['Node ',int2str(n-1)]);
                fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).INTEGRATION{n,6}(1));
                fprintf(fid_data,'%s',' ');
                fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).INTEGRATION{n,6}(2));
                fprintf(fid_data,'%s',' ');
                fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).INTEGRATION{n,5});
                table=BODIES_STATIC(Body,1).INTEGRATION{n,1}-1;
                for i=1:size(BODIES_STATIC(Body,1).INTEGRATION{n,1},1)
                    table(i,2:3)=BODIES_STATIC(Body,1).INTEGRATION{n,2}(2*i-1:2*i,1)'-1;
                end
                table=cat(2,table,BODIES_STATIC(Body,1).INTEGRATION{n,3},BODIES_STATIC(Body,1).INTEGRATION{n,4});
                fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).INTEGRATION{n,1},1));
                fprintf(fid_data,'%0.16g %0.16g %0.16g %0.16g %0.16g\n',table');
                fprintf(fid_data,'%s\n',' ');
            end
        end

        %%% TRIANGULATION %%%
        fprintf(fid_data,'%s\n','TRIANGULATION');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_CELLS);
        table=cat(2,[1:BODIES_STATIC(Body,1).NUMBER_CELLS]',BODIES_STATIC(Body,1).TRIANGULATION);
        fprintf(fid_data,'%0.16g %0.16g %0.16g %0.16g\n',table'-1);
        fprintf(fid_data,'%s\n',' ');

        %%% DOF TO DISP %%%
        fprintf(fid_data,'%s\n','DOF_TO_DISP');
        fprintf(fid_data,'%s\n',' ');
        for node=1:BODIES_STATIC(Body,1).NUMBER_NODES
            fprintf(fid_data,'%s\n',['NODE ',int2str(node-1)]);
            fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).DOF_TO_DISP{node,1},1));
            table=BODIES_STATIC(Body,1).DOF_TO_DISP{node,1}(:,1:4);
            table(:,1)=table(:,1)-1;
            fprintf(fid_data,'%0.16g %0.16g %0.16g %0.16g\n',table');
            fprintf(fid_data,'%s\n',' ');
        end
        fprintf(fid_data,'%s\n',' ');

        %%% STIFFNESS AND DAMPING MATRICES %%%
        fprintf(fid_data,'%s\n','MATRICES');
        %table=[];
        table=zeros(1e7,4);
        Nline=1;
        for i=1:size(BODIES_STATIC(Body,1).STIFFNESS_MATRIX,1)
            t1=full(BODIES_STATIC(Body,1).STIFFNESS_MATRIX{i,1})';
            t2=full(BODIES_STATIC(Body,1).DAMPING_MATRIX{i,1})';
            if isempty(t2)
                t2=ones(size(t1,1),2);
            end
            t=[t1,t2(:,2)];
            list=find(t(:,1)>=i);
            table(Nline:Nline+size(list,1)-1,:)=cat(2,i*ones(size(list,1),1),t(list,:));
            Nline=Nline+size(list,1);
        end
        table=table(1:Nline-1,:);
        fprintf(fid_data,'%0.16g\n',size(table,1));
        fprintf(fid_data,'%s\n',' ');
        table(:,1:2)=table(:,1:2)-1;
        fprintf(fid_data,'%0.16g %0.16g %0.16g %0.16g\n',table');
        fprintf(fid_data,'%s\n',' ');
        fprintf(fid_data,'%s\n',' ');
        fprintf(fid_data,'%s\n',' ');
    elseif strcmp(BODIES_STATIC(Body,1).CLASS,'Rigid')
        fprintf(fid_data,'%s\n',' ');
        fprintf(fid_data,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid_data,'%s\n',['%%%%%%%             BODY ',int2str(Body-1),'            %%%%%%%']);
        fprintf(fid_data,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid_data,'%s\n',' ');

        %%% PROPERTIES %%%
        fprintf(fid_data,'%s\n','RIGID');
        fprintf(fid_data,'%0.16g\n',Body-1);
        fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).CENTRE_OF_MASS);
        fprintf(fid_data,'%0.16g %0.16g %0.16g %0.16g\n',cat(2,BODIES_STATIC(Body,1).MASSES,BODIES_STATIC(Body,1).INVERSE_MASSES));
        fprintf(fid_data,'%s\n',' ');
        
        %%% NODES %%%
        fprintf(fid_data,'%s\n','NODES');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_NODES);
        fprintf(fid_data,'%s\n',' ');
        table=cat(2,[1:BODIES_STATIC(Body,1).NUMBER_NODES]'-1,BODIES_STATIC(Body,1).INITIAL_POSITIONS);
        fprintf(fid_data,'%0.16g %0.16g %0.16g\n',table');
        fprintf(fid_data,'%s\n',' ');

        %%% BORDERS %%%
        fprintf(fid_data,'%s\n','BORDERS');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_BORDERS);
        fprintf(fid_data,'%s\n',' ');
        for Border=1:BODIES_STATIC(Body,1).NUMBER_BORDERS
            fprintf(fid_data,'%s',BODIES_STATIC(Body,1).BORDERS{Border,1});
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%s\n',BODIES_STATIC(Body,1).BORDERS{Border,4});
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).BORDERS{Border,2});
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).BORDERS{Border,5});
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).BORDERS{Border,6});
            fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).BORDERS{Border,3}-1);
            fprintf(fid_data,'%s\n',' ');
        end
        
        %%% TRIANGULATION %%%
        fprintf(fid_data,'%s\n','TRIANGULATION');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_CELLS);
        table=cat(2,[1:BODIES_STATIC(Body,1).NUMBER_CELLS]',BODIES_STATIC(Body,1).TRIANGULATION);
        fprintf(fid_data,'%0.16g %0.16g %0.16g %0.16g\n',table'-1);
        fprintf(fid_data,'%s\n',' ');
        
    end
end
fclose(fid_data);


file_control=[Simulation_Name,'/STATIC_CONTROL.asc'];
fid_data = fopen(file_control,'w');

fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid_data,'%s\n','%%%%%%%          GENERAL DATA          %%%%%%%');
fprintf(fid_data,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid_data,'%s\n',' ');

%%% SIMULATION NAME %%%
fprintf(fid_data,'%s\n','SIMULATION_NAME');
fprintf(fid_data,'%s\n',Simulation_Name);
fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n',' ');

%%% MATERIALS %%%
fprintf(fid_data,'%s\n','MATERIALS');
fprintf(fid_data,'%0.16g\n',size(MATERIALS,1));
fprintf(fid_data,'%s\n',' ');
for i=1:size(MATERIALS,1)
    fprintf(fid_data,'%s\n',MATERIALS{i,1});
    fprintf(fid_data,'%s\n',MATERIALS{i,2});
    if isempty(MATERIALS{i,3})==0
        for j=1:size(MATERIALS{i,3},2)
            fprintf(fid_data,'%0.16g',MATERIALS{i,3}(1,j));
            fprintf(fid_data,'%s',' ');
        end
        fprintf(fid_data,'%s\n',' ');
    end
    fprintf(fid_data,'%s\n',' ');
end
fprintf(fid_data,'%s\n',' ');

%%% CONTACT LAWS %%%
fprintf(fid_data,'%s\n','CONTACT_LAWS');
fprintf(fid_data,'%0.16g\n',size(CONTACT_LAWS,1));
fprintf(fid_data,'%s\n',' ');
for i=1:size(CONTACT_LAWS,1)
    fprintf(fid_data,'%s\n',CONTACT_LAWS{i,1});
    fprintf(fid_data,'%s\n',CONTACT_LAWS{i,2});
    fprintf(fid_data,'%s\n',CONTACT_LAWS{i,3});
    fprintf(fid_data,'%s\n',CONTACT_LAWS{i,4});
    %if strcmp(CONTACT_LAWS{i,3},'CZMlinear')==0 & strcmp(CONTACT_LAWS{i,3},'CZMfatigue')==0
    %    fprintf(fid_data,'%0.16g',PENALTY_PARAMETER);
    %    fprintf(fid_data,'%s',' ');
    %    fprintf(fid_data,'%0.16g',PENALTY_PARAMETER);
    %    fprintf(fid_data,'%s',' ');
    %end
    if strcmp(CONTACT_LAWS{i,3},'Frictionless')==0
        for j=1:size(CONTACT_LAWS{i,5},2)
            fprintf(fid_data,'%0.16g',CONTACT_LAWS{i,5}(1,j));
            fprintf(fid_data,'%s',' ');
        end
    end
    fprintf(fid_data,'%s\n',' ');
    fprintf(fid_data,'%s\n',' ');
end
fprintf(fid_data,'%s\n',' ');

%%% SOLVER %%%
fprintf(fid_data,'%s\n','SOLVER');
fprintf(fid_data,'%s\n',SCHEME);
fprintf(fid_data,'%0.16g',TIME_INI);
fprintf(fid_data,'%s',' ');
fprintf(fid_data,'%0.16g',TIME_STEP);
fprintf(fid_data,'%s',' ');
fprintf(fid_data,'%0.16g\n',TIME_END);
fprintf(fid_data,'%0.16g',TARGET_ERROR);
fprintf(fid_data,'%s',' ');
fprintf(fid_data,'%0.16g',CONTROL_PARAMETER);
fprintf(fid_data,'%s',' ');
fprintf(fid_data,'%0.16g\n',ACCEPTED_RATIO);
fprintf(fid_data,'%0.16g',SAVE_PERIOD);
fprintf(fid_data,'%s',' ');
fprintf(fid_data,'%0.16g',PRINT_PERIOD);
fprintf(fid_data,'%s',' ');
fprintf(fid_data,'%0.16g\n',CONTACT_UPDATING_PERIOD);
fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n',' ');

%%% PERIODICITY %%%
fprintf(fid_data,'%s\n','PERIODIC_BOUNDARIES');
fprintf(fid_data,'%0.16g',PERIODIC_BOUNDARIES(1));
fprintf(fid_data,'%s',' ');
fprintf(fid_data,'%0.16g\n',PERIODIC_BOUNDARIES(2));
fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n',' ');

%%% GRAVITY %%%
fprintf(fid_data,'%s\n','GRAVITY');
fprintf(fid_data,'%0.16g',GRAVITY(1));
fprintf(fid_data,'%s',' ');
fprintf(fid_data,'%0.16g\n',GRAVITY(2));
fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n',' ');

%%% PLOTTING OPTIONS %%%
%fprintf(fid_data,'%s\n','PLOTTING');
%fprintf(fid_data,'%0.16g\n',ACTIVATE_PLOT);
%fprintf(fid_data,'%0.16g',GRAPHIC_BOUNDARIES(1));
%fprintf(fid_data,'%s',' ');
%fprintf(fid_data,'%10.10g',GRAPHIC_BOUNDARIES(2));
%fprintf(fid_data,'%s',' ');
%fprintf(fid_data,'%0.16g',GRAPHIC_BOUNDARIES(3));
%fprintf(fid_data,'%s',' ');
%fprintf(fid_data,'%0.16g\n',GRAPHIC_BOUNDARIES(4));
%fprintf(fid_data,'%s\n',' ');
%fprintf(fid_data,'%s\n',' ');

%%% NUMBER_BODIES %%%
fprintf(fid_data,'%s\n','NUMBER_BODIES');
fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC,1));
fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n',' ');

%%% MONITORINGS %%%
fprintf(fid_data,'%s\n','MONITORING');
if isempty(MONITORINGS{1,1})
    fprintf(fid_data,'%s\n','0');
else
    fprintf(fid_data,'%0.16g\n',size(MONITORINGS,1));
    for i=1:size(MONITORINGS,1)
        fprintf(fid_data,'%s\n',MONITORINGS{i,1});
    end
end
fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n',' ');

%%% SPIES %%%
fprintf(fid_data,'%s\n','SPIES');
if isempty(SPIES{1,1})
    fprintf(fid_data,'%s\n','0');
else
    fprintf(fid_data,'%0.16g\n',size(SPIES,1));
    for n=1:size(SPIES,1)
        fprintf(fid_data,'%s\n',' ');
        fprintf(fid_data,'%s',SPIES{n,1});
        fprintf(fid_data,'%s',' ');
        fprintf(fid_data,'%0.16g',SPIES{n,2});
        fprintf(fid_data,'%s',' ');
        fprintf(fid_data,'%0.16g\n',SPIES{n,3});
        for i=1:size(SPIES{n,4},1)
            fprintf(fid_data,'%s\n',SPIES{n,4}{i,1});
        end
    end
end
fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n',' ');

%%% DEACTIVATIONS %%%
fprintf(fid_data,'%s\n','DEACTIVATION');
if isempty(DEACTIVATIONS{1,1})
    fprintf(fid_data,'%s\n','0');
else
    fprintf(fid_data,'%0.16g\n',size(DEACTIVATIONS,1));
    for i=1:size(DEACTIVATIONS,1)
        fprintf(fid_data,'%s\n',DEACTIVATIONS{i,1});
    end
end
fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n',' ');

%%% GRAPHIC %%%
fprintf(fid_data,'%s\n','GRAPHIC');
fprintf(fid_data,'%0.16g',To_Plot(1));fprintf(fid_data,'%s\n',' Body Index');
fprintf(fid_data,'%0.16g',To_Plot(2));fprintf(fid_data,'%s\n',' Initial Position');
fprintf(fid_data,'%0.16g',To_Plot(3));fprintf(fid_data,'%s\n',' Current Position');
fprintf(fid_data,'%0.16g',To_Plot(4));fprintf(fid_data,'%s\n',' Displacement');
fprintf(fid_data,'%0.16g',To_Plot(5));fprintf(fid_data,'%s\n',' Velocity');
fprintf(fid_data,'%0.16g',To_Plot(6));fprintf(fid_data,'%s\n',' Acceleration');
fprintf(fid_data,'%0.16g',To_Plot(7));fprintf(fid_data,'%s\n',' Force');
fprintf(fid_data,'%0.16g',To_Plot(8));fprintf(fid_data,'%s\n',' Internal Force');
fprintf(fid_data,'%0.16g',To_Plot(9));fprintf(fid_data,'%s\n',' Contact Force');
fprintf(fid_data,'%0.16g',To_Plot(10));fprintf(fid_data,'%s\n',' Body Force');
fprintf(fid_data,'%0.16g',To_Plot(11));fprintf(fid_data,'%s\n',' Dirichlet Force');
fprintf(fid_data,'%0.16g',To_Plot(12));fprintf(fid_data,'%s\n',' Neumann Force');
fprintf(fid_data,'%0.16g',To_Plot(13));fprintf(fid_data,'%s\n',' Damping Force');
fprintf(fid_data,'%0.16g',To_Plot(14));fprintf(fid_data,'%s\n',' Alid Force');
fprintf(fid_data,'%0.16g',To_Plot(15));fprintf(fid_data,'%s\n',' Jacobian');
fprintf(fid_data,'%0.16g',To_Plot(16));fprintf(fid_data,'%s\n',' Cauchy XX Stress');
fprintf(fid_data,'%0.16g',To_Plot(17));fprintf(fid_data,'%s\n',' Cauchy YY Stress');
fprintf(fid_data,'%0.16g',To_Plot(18));fprintf(fid_data,'%s\n',' Cauchy XY Stress');
fprintf(fid_data,'%0.16g',To_Plot(19));fprintf(fid_data,'%s\n',' Cauchy ZZ Stress');
fprintf(fid_data,'%0.16g',To_Plot(20));fprintf(fid_data,'%s\n',' Tresca Stress');
fprintf(fid_data,'%0.16g',To_Plot(21));fprintf(fid_data,'%s\n',' Von Mises Stress');
fprintf(fid_data,'%0.16g',To_Plot(22));fprintf(fid_data,'%s\n',' Major Principal Stress');
fprintf(fid_data,'%0.16g',To_Plot(23));fprintf(fid_data,'%s\n',' Intermediate Principal Stress');
fprintf(fid_data,'%0.16g',To_Plot(24));fprintf(fid_data,'%s\n',' Minor Principal Stress');
fprintf(fid_data,'%0.16g',To_Plot(25));fprintf(fid_data,'%s\n',' Spherical Stress');
fprintf(fid_data,'%0.16g',To_Plot(26));fprintf(fid_data,'%s\n',' Green-Lagrange XX strain');
fprintf(fid_data,'%0.16g',To_Plot(27));fprintf(fid_data,'%s\n',' Green-Lagrange YY strain');
fprintf(fid_data,'%0.16g',To_Plot(28));fprintf(fid_data,'%s\n',' Green-Lagrange XY strain');
fprintf(fid_data,'%0.16g',To_Plot(29));fprintf(fid_data,'%s\n',' Norm of the Green-Lagrange strain tensor');
fprintf(fid_data,'%0.16g',To_Plot(30));fprintf(fid_data,'%s\n',' Body Damage');
fprintf(fid_data,'%0.16g',To_Plot(31));fprintf(fid_data,'%s\n',' Normalized Displacement Error');
fprintf(fid_data,'%0.16g',To_Plot(32));fprintf(fid_data,'%s\n',' Displacement Error');
fprintf(fid_data,'%s\n',' ');
fprintf(fid_data,'%s\n',' ');

%%% BODIES %%%
for Body=1:size(BODIES_STATIC,1)
    if strcmp(BODIES_STATIC(Body,1).CLASS,'Body')
        fprintf(fid_data,'%s\n',' ');
        fprintf(fid_data,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid_data,'%s\n',['%%%%%%%              BODY ',int2str(Body-1),'            %%%%%%%']);
        fprintf(fid_data,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid_data,'%s\n',' ');

        %%% PROPERTIES %%%
        fprintf(fid_data,'%s\n','DEFORMABLE');
        fprintf(fid_data,'%0.16g\n',Body-1);
        fprintf(fid_data,'%s\n',BODIES_STATIC(Body,1).MATERIAL);
        fprintf(fid_data,'%s\n',BODIES_STATIC(Body,1).TYPE);
        fprintf(fid_data,'%0.16g %0.16g %0.16g\n',[BODIES_STATIC(Body,1).NODAL_DISTANCE,BODIES_STATIC(Body,1).DETECTIONS]);
        fprintf(fid_data,'%s\n',' ');

        %%% NODES %%%
        fprintf(fid_data,'%s\n','NODES');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_NODES);
        fprintf(fid_data,'%s\n',' ');

        %%% BORDERS %%%
        fprintf(fid_data,'%s\n','BORDERS');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_BORDERS);
        fprintf(fid_data,'%s\n',' ');
        for Border=1:BODIES_STATIC(Body,1).NUMBER_BORDERS
            fprintf(fid_data,'%s',BODIES_STATIC(Body,1).BORDERS{Border,1});
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%s\n',BODIES_STATIC(Body,1).BORDERS{Border,4});
            fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).BORDERS{Border,2});
            if isempty(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1}) & isempty(BODIES_STATIC(Body,1).NEUMANN_BC{Border,1})
                fprintf(fid_data,'%s\n',['X None']);
            elseif isempty(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1})==0
                if strcmp(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,2},'Driven')
                    fprintf(fid_data,'%s\n',['X Dirichlet Driven']);
                elseif strcmp(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,2},'Soft')
                    fprintf(fid_data,'%s',['X Dirichlet Soft ']);
                    fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,3}(1));
                    fprintf(fid_data,'%s',[' ']);
                    fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,3}(2));
                end
                fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1},1));
                fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1}');
            elseif isempty(BODIES_STATIC(Body,1).NEUMANN_BC{Border,1})==0
                if strcmp(BODIES_STATIC(Body,1).NEUMANN_BC{Border,2},'Oriented')
                    fprintf(fid_data,'%s\n',['X Neumann Oriented']);
                elseif strcmp(BODIES_STATIC(Body,1).NEUMANN_BC{Border,2},'Following')
                    fprintf(fid_data,'%s\n',['X Neumann Following']);
                end
                fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).NEUMANN_BC{Border,1},1));
                fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).NEUMANN_BC{Border,1}');
            end
            if isempty(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,4}) & isempty(BODIES_STATIC(Body,1).NEUMANN_BC{Border,3})
                fprintf(fid_data,'%s\n',['Y None']);
            elseif isempty(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,4})==0
                if strcmp(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,5},'Driven')
                    fprintf(fid_data,'%s\n',['Y Dirichlet Driven']);
                elseif strcmp(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,5},'Soft')
                    fprintf(fid_data,'%s',['Y Dirichlet Soft ']);
                    fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,6}(1));
                    fprintf(fid_data,'%s',[' ']);
                    fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,6}(2));
                end
                fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,4},1));
                fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,4}');
            elseif isempty(BODIES_STATIC(Body,1).NEUMANN_BC{Border,3})==0
                if strcmp(BODIES_STATIC(Body,1).NEUMANN_BC{Border,4},'Oriented')
                    fprintf(fid_data,'%s\n',['Y Neumann Oriented']);
                elseif strcmp(BODIES_STATIC(Body,1).NEUMANN_BC{Border,4},'Following')
                    fprintf(fid_data,'%s\n',['Y Neumann Following']);
                end
                fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).NEUMANN_BC{Border,3},1));
                fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).NEUMANN_BC{Border,3}');
            end
            fprintf(fid_data,'%s\n',' '); 
                
                
                
            %fprintf(fid_data,'%s\n',['Dirichlet X']);
            %if isempty(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1})==0
            %    fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1},1));
            %    fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1}');
            %else
            %    fprintf(fid_data,'%0.16g\n',0);
            %end
            %fprintf(fid_data,'%s\n',['Dirichlet Y']);
            %if isempty(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,2})==0
            %    fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,2},1));
            %    fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,2}');
            %else
            %    fprintf(fid_data,'%0.16g\n',0);
            %end
            %fprintf(fid_data,'%s\n',['Neumann X']);
            %if isempty(BODIES_STATIC(Body,1).NEUMANN_BC{Border,1})==0
            %    fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).NEUMANN_BC{Border,1},1));
            %    fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).NEUMANN_BC{Border,1}');
            %else
            %    fprintf(fid_data,'%0.16g\n',0);
            %end
            %fprintf(fid_data,'%s\n',['Neumann Y']);
            %if isempty(BODIES_STATIC(Body,1).NEUMANN_BC{Border,2})==0
            %    fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).NEUMANN_BC{Border,2},1));
            %    fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).NEUMANN_BC{Border,2}');
            %else
            %    fprintf(fid_data,'%0.16g\n',0);
            %end
            %fprintf(fid_data,'%s\n',' ');
        end
        
        %%% ALID %%%
        if size(BODIES_STATIC(Body,1).ALID{1,1},1)>0
            fprintf(fid_data,'%s\n','ALID');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,1}(1));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,1}(2));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,1}(3));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,1}(4));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,1}(5));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,1}(6));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,1}(7));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,1}(8));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,1}(9));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,1}(10));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,1}(11));
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).ALID{1,1}(12));
            if BODIES_STATIC(Body,1).ALID{1,1}(1)>0
                fprintf(fid_data,'%s\n',' ');
                for line=1:BODIES_STATIC(Body,1).ALID{1,1}(1)
                    fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,2}(line,1));
                    fprintf(fid_data,'%s',' ');
                    fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).ALID{1,2}(line,2));
                end
            end
            if BODIES_STATIC(Body,1).ALID{1,1}(2)>0
                fprintf(fid_data,'%s\n',' ');
                for line=1:BODIES_STATIC(Body,1).ALID{1,1}(2)
                    fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).ALID{1,3}(line,1));
                    fprintf(fid_data,'%s',' ');
                    fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).ALID{1,3}(line,2));
                end
            end
            fprintf(fid_data,'%s\n',' ');
        end
        
        %%% INTEGRATION %%%
        fprintf(fid_data,'%s\n','INTEGRATION');
        fprintf(fid_data,'%0.16g',size(BODIES_STATIC(Body,1).INTEGRATION,1));
        fprintf(fid_data,'%s',' ');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_REGIONS);
        fprintf(fid_data,'%s\n',' ');

        %%% TRIANGULATION %%%
        fprintf(fid_data,'%s\n','TRIANGULATION');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_CELLS);
        fprintf(fid_data,'%s\n',' ');
        fprintf(fid_data,'%s\n',' ');
        fprintf(fid_data,'%s\n',' ');
    elseif strcmp(BODIES_STATIC(Body,1).CLASS,'Rigid')
        fprintf(fid_data,'%s\n',' ');
        fprintf(fid_data,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid_data,'%s\n',['%%%%%%%             BODY ',int2str(Body-1),'            %%%%%%%']);
        fprintf(fid_data,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid_data,'%s\n',' ');

        %%% PROPERTIES %%%
        fprintf(fid_data,'%s\n','RIGID');
        fprintf(fid_data,'%0.16g\n',Body-1);
        fprintf(fid_data,'%s\n',BODIES_STATIC(Body,1).MATERIAL);
        fprintf(fid_data,'%s\n',BODIES_STATIC(Body,1).TYPE);
        %fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NODAL_DISTANCE);
        fprintf(fid_data,'%0.16g %0.16g %0.16g\n',[BODIES_STATIC(Body,1).NODAL_DISTANCE,BODIES_STATIC(Body,1).DETECTIONS]);
        fprintf(fid_data,'%s\n',' ');
        
        %%% NODES %%%
        fprintf(fid_data,'%s\n','NODES');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_NODES);
        fprintf(fid_data,'%s\n',' ');

        %%% BORDERS %%%
        fprintf(fid_data,'%s\n','BORDERS');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_BORDERS);
        fprintf(fid_data,'%s\n',' ');
        for Border=1:BODIES_STATIC(Body,1).NUMBER_BORDERS
            fprintf(fid_data,'%s',BODIES_STATIC(Body,1).BORDERS{Border,1});
            fprintf(fid_data,'%s',' ');
            fprintf(fid_data,'%s\n',BODIES_STATIC(Body,1).BORDERS{Border,4});
            fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).BORDERS{Border,2});
            if isempty(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1}) & isempty(BODIES_STATIC(Body,1).NEUMANN_BC{Border,1})
                fprintf(fid_data,'%s\n',['X None']);
            elseif isempty(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1})==0
                if strcmp(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,2},'Driven')
                    fprintf(fid_data,'%s\n',['X Dirichlet Driven']);
                elseif strcmp(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,2},'Soft')
                    fprintf(fid_data,'%s',['X Dirichlet Soft ']);
                    fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,3}(1));
                    fprintf(fid_data,'%s',[' ']);
                    fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,3}(2));
                end
                fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1},1));
                fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1}');
            elseif isempty(BODIES_STATIC(Body,1).NEUMANN_BC{Border,1})==0
                if strcmp(BODIES_STATIC(Body,1).NEUMANN_BC{Border,2},'Oriented')
                    fprintf(fid_data,'%s\n',['X Neumann Oriented']);
                elseif strcmp(BODIES_STATIC(Body,1).NEUMANN_BC{Border,2},'Following')
                    fprintf(fid_data,'%s\n',['X Neumann Following']);
                end
                fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).NEUMANN_BC{Border,1},1));
                fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).NEUMANN_BC{Border,1}');
            end
            if isempty(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,4}) & isempty(BODIES_STATIC(Body,1).NEUMANN_BC{Border,3})
                fprintf(fid_data,'%s\n',['Y None']);
            elseif isempty(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,4})==0
                if strcmp(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,5},'Driven')
                    fprintf(fid_data,'%s\n',['Y Dirichlet Driven']);
                elseif strcmp(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,5},'Soft')
                    fprintf(fid_data,'%s',['Y Dirichlet Soft ']);
                    fprintf(fid_data,'%0.16g',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,6}(1));
                    fprintf(fid_data,'%s',[' ']);
                    fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,6}(2));
                end
                fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,4},1));
                fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,4}');
            elseif isempty(BODIES_STATIC(Body,1).NEUMANN_BC{Border,3})==0
                if strcmp(BODIES_STATIC(Body,1).NEUMANN_BC{Border,4},'Oriented')
                    fprintf(fid_data,'%s\n',['Y Neumann Oriented']);
                elseif strcmp(BODIES_STATIC(Body,1).NEUMANN_BC{Border,4},'Following')
                    fprintf(fid_data,'%s\n',['Y Neumann Following']);
                end
                fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).NEUMANN_BC{Border,3},1));
                fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).NEUMANN_BC{Border,3}');
            end
            fprintf(fid_data,'%s\n',' ');
        end
        %fprintf(fid_data,'%s\n','BORDERS');
        %fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_BORDERS);
        %fprintf(fid_data,'%s\n',' ');
        %for Border=1:BODIES_STATIC(Body,1).NUMBER_BORDERS
        %    fprintf(fid_data,'%s',BODIES_STATIC(Body,1).BORDERS{Border,1});
        %    fprintf(fid_data,'%s',' ');
        %    fprintf(fid_data,'%s\n',BODIES_STATIC(Body,1).BORDERS{Border,4});
        %    fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).BORDERS{Border,2});
        %    fprintf(fid_data,'%s\n',['Dirichlet X']);
        %    if isempty(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1})==0
        %        fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1},1));
        %        fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,1}');
        %    else
        %        fprintf(fid_data,'%0.16g\n',0);
        %    end
        %    fprintf(fid_data,'%s\n',['Dirichlet Y']);
        %    if isempty(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,2})==0
        %        fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).DIRICHLET_BC{Border,2},1));
        %        fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).DIRICHLET_BC{Border,2}');
        %    else
        %        fprintf(fid_data,'%0.16g\n',0);
        %    end
        %    fprintf(fid_data,'%s\n',['Neumann X']);
        %    if isempty(BODIES_STATIC(Body,1).NEUMANN_BC{Border,1})==0
        %        fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).NEUMANN_BC{Border,1},1));
        %        fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).NEUMANN_BC{Border,1}');
        %    else
        %        fprintf(fid_data,'%0.16g\n',0);
        %    end
        %    fprintf(fid_data,'%s\n',['Neumann Y']);
        %    if isempty(BODIES_STATIC(Body,1).NEUMANN_BC{Border,2})==0
        %        fprintf(fid_data,'%0.16g\n',size(BODIES_STATIC(Body,1).NEUMANN_BC{Border,2},1));
        %        fprintf(fid_data,'%0.16g %0.16g\n',BODIES_STATIC(Body,1).NEUMANN_BC{Border,2}');
        %    else
        %        fprintf(fid_data,'%0.16g\n',0);
        %    end
        %    fprintf(fid_data,'%s\n',' ');
        %end
        
        %%% TRIANGULATION %%%
        fprintf(fid_data,'%s\n','TRIANGULATION');
        fprintf(fid_data,'%0.16g\n',BODIES_STATIC(Body,1).NUMBER_CELLS);
        fprintf(fid_data,'%s\n',' ');
        fprintf(fid_data,'%s\n',' ');
        fprintf(fid_data,'%s\n',' ');
    end
end
fclose(fid_data);