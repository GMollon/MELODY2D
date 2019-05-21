if NUMBER_SAVE<10
    file=[Simulation_Name,'/DYNAMIC_0000',int2str(NUMBER_SAVE),'.asc'];
elseif NUMBER_SAVE<100
    file=[Simulation_Name,'/DYNAMIC_000',int2str(NUMBER_SAVE),'.asc'];
elseif NUMBER_SAVE<1000
    file=[Simulation_Name,'/DYNAMIC_00',int2str(NUMBER_SAVE),'.asc'];Initialize_CZM=0;
elseif NUMBER_SAVE<10000
    file=[Simulation_Name,'/DYNAMIC_0',int2str(NUMBER_SAVE),'.asc'];
elseif NUMBER_SAVE<100000
    file=[Simulation_Name,'/DYNAMIC_',int2str(NUMBER_SAVE),'.asc'];
end
fid = fopen(file,'w');

fprintf(fid,'%s\n',' ');
fprintf(fid,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid,'%s\n','%%%%%%%          GENERAL DATA          %%%%%%%');
fprintf(fid,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid,'%s\n',' ');

%%% SIMULATION NAME %%%
fprintf(fid,'%s\n','SIMULATION_NAME');
fprintf(fid,'%s\n',Simulation_Name);
fprintf(fid,'%s\n',' ');

%%% TIME %%%
fprintf(fid,'%s\n','SOLVER');
fprintf(fid,'%s\n',SCHEME);
fprintf(fid,'%0.16g',NUMBER_SAVE);
fprintf(fid,'%s',' ');
fprintf(fid,'%0.16g',NUMBER_ITERATIONS);
fprintf(fid,'%s',' ');
fprintf(fid,'%0.16g',TIME);
fprintf(fid,'%s',' ');
fprintf(fid,'%0.16g',NUMBER_PRINT);
fprintf(fid,'%s\n',' ');
fprintf(fid,'%0.16g\n',TIME_STEP);
fprintf(fid,'%0.16g',NEXT_SAVE);
fprintf(fid,'%s',' ');
fprintf(fid,'%0.16g',NEXT_PRINT);
fprintf(fid,'%s',' ');
fprintf(fid,'%0.16g',NEXT_UPDATE);
fprintf(fid,'%s\n',' ');
if Initialize_CZM==1
    fprintf(fid,'%s\n','INITIALIZE_CZM');
end
fprintf(fid,'%s\n','UPDATE_MASS_MATRIX');
fprintf(fid,'%s\n','UPDATE_DAMPING_MATRIX');
fprintf(fid,'%s\n',' ');

%%% BODIES %%%
for Body=1:size(BODIES_DYNAMIC,1)
    if strcmp(BODIES_STATIC(Body,1).CLASS,'Body')
        fprintf(fid,'%s\n',' ');
        fprintf(fid,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid,'%s\n',['%%%%%%%              BODY ',int2str(Body-1),'            %%%%%%%']);
        fprintf(fid,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid,'%s\n',' ');
        fprintf(fid,'%s\n','DEFORMABLE');
        fprintf(fid,'%s\n',BODIES_STATIC(Body,1).STATUS);
        fprintf(fid,'%0.16g\n',Body-1);
        fprintf(fid,'%s\n',' ');

        %%% KINEMATICS %%%
        fprintf(fid,'%s\n','KINEMATICS');
        table=cat(2,BODIES_DYNAMIC(Body,1).CURRENT_POSITIONS,BODIES_DYNAMIC(Body,1).DISPLACEMENTS,BODIES_DYNAMIC(Body,1).VELOCITIES,BODIES_DYNAMIC(Body,1).ACCELERATIONS);
        table=cat(2,[1:size(table,1)]',table,BODIES_DYNAMIC(Body,1).DISPLACEMENTS_PARAMETERS,BODIES_DYNAMIC(Body,1).VELOCITIES_PARAMETERS,BODIES_DYNAMIC(Body,1).ACCELERATIONS_PARAMETERS);
        fprintf(fid,'%0.16g\n',size(table,1));
        table(:,1)=table(:,1)-1;
        fprintf(fid,'%0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g\n',table');
        fprintf(fid,'%s\n',' ');

        %%% FORCES %%%
        fprintf(fid,'%s\n','FORCES');
        table=cat(2,BODIES_DYNAMIC(Body,1).INTERNAL_FORCES,BODIES_DYNAMIC(Body,1).CONTACT_FORCES,BODIES_DYNAMIC(Body,1).SELF_CONTACT_FORCES,BODIES_DYNAMIC(Body,1).BODY_FORCES);
        table=cat(2,[1:size(table,1)]',table,BODIES_DYNAMIC(Body,1).DIRICHLET_FORCES,BODIES_DYNAMIC(Body,1).NEUMANN_FORCES,BODIES_DYNAMIC(Body,1).DAMPING_FORCES,BODIES_DYNAMIC(Body,1).FORCES);
        table(:,1)=table(:,1)-1;
        fprintf(fid,'%0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g\n',table');
        fprintf(fid,'%s\n',' ');

        %%% NEIGHBOURS %%%
        fprintf(fid,'%s\n','NEIGHBOURS');
        fprintf(fid,'%0.16g\n',BODIES_DYNAMIC(Body,1).NUMBER_NEIGHBOURS);
        if BODIES_DYNAMIC(Body,1).NUMBER_NEIGHBOURS>0
            table=BODIES_DYNAMIC(Body,1).NEIGHBOURS;
            table(:,1:2)=table(:,1:2)-1;
            fprintf(fid,'%0.16g %0.16g %0.16g\n',table');
        end
        fprintf(fid,'%s\n',' ');

        %%% BORDERS %%%
        fprintf(fid,'%s\n','BORDERS');
        fprintf(fid,'%s\n',' ');
        NUMBER_BORDERS=size(BODIES_DYNAMIC(Body,1).BC,1);
        for Border=1:NUMBER_BORDERS
            fprintf(fid,'%0.16g\n',Border);
            table=BODIES_DYNAMIC(Body,1).CONTACTS{Border,1};
            table=cat(2,table(:,5),zeros(size(table,1),1),table(:,1:4),table(:,6:20));
            for i=1:size(table,1)
                if table(i,3)>0
                    table(i,2)=1;
                end
            end
            for i=1:size(table,1)
                if sum(table(i,2:21))~=0
                    table(i,3:4)=table(i,3:4)-1;
                    table(i,6)=table(i,6)-1;
                    table(i,14)=table(i,14)-1;
                    table(i,16)=table(i,16)-1;
                    table(i,18)=table(i,18)-1;
                    table(i,20)=table(i,20)-1;
                end
                if table(i,1)==0
                    if BODIES_DYNAMIC(Body,1).CONTACT_PRESSURES{Border,1}(i,1)~=0
                        table(i,1)=BODIES_DYNAMIC(Body,1).CONTACT_PRESSURES{Border,1}(i,1);
                    elseif isempty(BODIES_DYNAMIC(Body,1).BC_PRESSURES{Border,1})==0
                        if BODIES_DYNAMIC(Body,1).BC_PRESSURES{Border,1}(i,1)~=0
                            table(i,1)=BODIES_DYNAMIC(Body,1).BC_PRESSURES{Border,1}(i,1);
                        end
                    end
                end
            end
            table=table(:,[1,13,14,15,16,17,18,19,20]);
            fprintf(fid,'%0.16g\n',size(table,1));
            fprintf(fid,'%0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g %0.16g\n',table');
            fprintf(fid,'%s\n',' ');
        end

        %%% CONTACTS %%%
        fprintf(fid,'%s\n','CONTACTS');
        fprintf(fid,'%s\n','0');
        fprintf(fid,'%s\n',' ');

        %%% CONTACT PRESSURES %%%
        fprintf(fid,'%s\n','CONTACT_PRESSURES');
        fprintf(fid,'%s\n',' ');
        for Border=1:NUMBER_BORDERS
            fprintf(fid,'%0.16g\n',Border);
            table=BODIES_DYNAMIC(Body,1).CONTACT_PRESSURES{Border,1};
            fprintf(fid,'%0.16g\n',size(table,1));
            fprintf(fid,'%0.16g %0.16g\n',table(:,2:3)');
            fprintf(fid,'%s\n',' ');
        end

        %%% BC PRESSURES %%%
        fprintf(fid,'%s\n','BC_PRESSURES');
        fprintf(fid,'%s\n',' ');
        for Border=1:NUMBER_BORDERS
            fprintf(fid,'%0.16g\n',Border);
            if isempty(BODIES_DYNAMIC(Body).BC_PRESSURES{Border,1})
                fprintf(fid,'%0.16g\n',size(BODIES_DYNAMIC(Body,1).PROXIMITIES{Border,1},1));
                table=zeros(size(BODIES_DYNAMIC(Body,1).PROXIMITIES{Border,1},1),3);
                fprintf(fid,'%0.16g %0.16g\n',table(:,2:3)');
            else
                table=BODIES_DYNAMIC(Body).BC_PRESSURES{Border,1};
                fprintf(fid,'%0.16g\n',size(table,1));
                fprintf(fid,'%0.16g %0.16g\n',table(:,2:3)');
            end
            fprintf(fid,'%s\n',' ');
        end
    elseif strcmp(BODIES_STATIC(Body,1).CLASS,'Rigid')
        fprintf(fid,'%s\n',' ');
        fprintf(fid,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid,'%s\n',['%%%%%%%             BODY ',int2str(Body-1),'            %%%%%%%']);
        fprintf(fid,'%s\n','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        fprintf(fid,'%s\n',' ');
        fprintf(fid,'%s\n','RIGID');
        fprintf(fid,'%s\n','active');
        fprintf(fid,'%0.16g\n',Body-1);
        fprintf(fid,'%s\n',' ');

        %%% KINEMATICS %%%
        fprintf(fid,'%s\n','KINEMATICS');
        fprintf(fid,'%0.16g %0.16g %0.16g\n',cat(2,BODIES_STATIC(Body,1).CENTRE_OF_MASS,0));
        fprintf(fid,'%0.16g %0.16g %0.16g\n',[0,0,0]);
        fprintf(fid,'%0.16g %0.16g %0.16g\n',cat(2,BODIES_DYNAMIC(Body,1).VELOCITIES(1,:),0));
        fprintf(fid,'%s\n',' ');

        %%% FORCES %%%
        fprintf(fid,'%s\n','FORCES');
        fprintf(fid,'%0.16g %0.16g %0.16g\n',[0,0,0]);
        fprintf(fid,'%0.16g %0.16g %0.16g\n',[0,0,0]);
        fprintf(fid,'%0.16g %0.16g %0.16g\n',[0,0,0]);
        fprintf(fid,'%0.16g %0.16g %0.16g\n',[0,0,0]);
        fprintf(fid,'%0.16g %0.16g %0.16g\n',[0,0,0]);
        fprintf(fid,'%0.16g %0.16g %0.16g\n',[0,0,0]);
        fprintf(fid,'%s\n',' ');

        %%% NEIGHBOURS %%%
        fprintf(fid,'%s\n','NEIGHBOURS');
        fprintf(fid,'%0.16g\n',BODIES_DYNAMIC(Body,1).NUMBER_NEIGHBOURS);
        if BODIES_DYNAMIC(Body,1).NUMBER_NEIGHBOURS>0
            table=BODIES_DYNAMIC(Body,1).NEIGHBOURS;
            table(:,1:2)=table(:,1:2)-1;
            fprintf(fid,'%0.16g %0.16g %0.16g\n',table');
        end
        fprintf(fid,'%s\n',' ');

        %%% CONTACTS %%%
        fprintf(fid,'%s\n','CONTACTS');
        fprintf(fid,'%s\n','0');
        fprintf(fid,'%s\n',' ');

        %%% CONTACT PRESSURES %%%
        fprintf(fid,'%s\n','CONTACT_PRESSURES');
        fprintf(fid,'%s\n',' ');
        for Border=1:BODIES_STATIC(Body,1).NUMBER_BORDERS
            fprintf(fid,'%0.16g\n',Border);
            table=BODIES_DYNAMIC(Body,1).CONTACT_PRESSURES{Border,1};
            fprintf(fid,'%0.16g\n',size(table,1));
            fprintf(fid,'%0.16g %0.16g\n',table(:,2:3)');
            fprintf(fid,'%s\n',' ');
        end
        fprintf(fid,'%s\n',' ');

        %%% BC PRESSURES %%%
        fprintf(fid,'%s\n','BC_PRESSURES');
        fprintf(fid,'%s\n',' ');
        for Border=1:BODIES_STATIC(Body,1).NUMBER_BORDERS
            fprintf(fid,'%0.16g\n',Border);
            if isempty(BODIES_DYNAMIC(Body).BC_PRESSURES{Border,1})
                fprintf(fid,'%0.16g\n',size(BODIES_DYNAMIC(Body,1).PROXIMITIES{Border,1},1));
                table=zeros(size(BODIES_DYNAMIC(Body,1).PROXIMITIES{Border,1},1),3);
                fprintf(fid,'%0.16g %0.16g\n',table(:,2:3)');
            else
                table=BODIES_DYNAMIC(Body).BC_PRESSURES{Border,1};
                fprintf(fid,'%0.16g\n',size(table,1));
                fprintf(fid,'%0.16g %0.16g\n',table(:,2:3)');
            end
            fprintf(fid,'%s\n',' ');
        end
    end
end
fclose(fid);