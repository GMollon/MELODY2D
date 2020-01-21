%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                       %%
%%                  MELODY PREPROCESSOR                  %%
%%            Multibody ELement-free Open code           %%
%%                 for DYnamic simulation                %%
%%                                                       %%
%%                     Main Program                      %%
%%               Version 3.11 ; May 2019                 %%
%%                                                       %%
%%                Author: Guilhem Mollon                 %%
%%                                                       %%
%%            Developed at LaMCoS, INSA Lyon,            %%
%%                     Lyon, France                      %%
%%                      Year 2015                        %%
%%                                                       %%
%%            Please read enclosed .pdf file             %%
%%                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%%% Preprocessing %%%
MELODYLoadData_Example;
global NUMBER_BODIES BODIES_STATIC BODIES_DYNAMIC BODIES_CURRENT MATERIALS CONTACT_LAWS GRAPHIC_BOUNDARIES PERIODIC_BOUNDARIES GRAVITY
mkdir(Simulation_Name);
disp(' ')
[NUMBER_BODIES,NUMBER_RIGIDS,BODIES_STATIC,BODIES_DYNAMIC,BODIES_CURRENT,MATERIALS,CONTACT_LAWS,TIME_INI,TIME_STEP,TIME_END,TIME,SAVE_PERIOD,PRINT_PERIOD,CONTACT_UPDATING_PERIOD,PERIODIC_BOUNDARIES,SCHEME,TARGET_ERROR,CONTROL_PARAMETER,ACCEPTED_RATIO,GRAVITY,MONITORINGS,DEACTIVATIONS,SPIES]=Preprocessing(Contours,Distributions,Interpolations,Integrations,Detections,Mesh_Ratios,Materials,Contact_Laws,Bodies_Materials,Imposed_Pressures,Imposed_Velocities,Initial_Velocities,Periodic_Boundaries,Gravity,Time_Stepping_Parameters,Save_Periods,Contact_Updating_Period,Scheme,Scheme_Parameters,Activate_Plot,Monitorings,Deactivations,Spies,Status,Alid);
NUMBER_ITERATIONS=0;
NUMBER_SAVE=0;
NUMBER_PRINT=0;
NEXT_PRINT=PRINT_PERIOD;
NEXT_SAVE=SAVE_PERIOD;
NEXT_UPDATE=CONTACT_UPDATING_PERIOD;   
MONITORING=[];
   
save([Simulation_Name,'/STATIC.mat'],'NUMBER_BODIES','NUMBER_RIGIDS','BODIES_STATIC','MATERIALS','CONTACT_LAWS','TIME_INI','TIME_STEP','TIME_END','SAVE_PERIOD','PRINT_PERIOD','CONTACT_UPDATING_PERIOD','PERIODIC_BOUNDARIES','SCHEME','TARGET_ERROR','CONTROL_PARAMETER','ACCEPTED_RATIO','GRAVITY','MONITORINGS','DEACTIVATIONS','SPIES','-v7.3');
save([Simulation_Name,'/',Simulation_Name,'_Simulation_',int2str(NUMBER_SAVE),'.mat'],'Simulation_Name','BODIES_DYNAMIC','TIME','TIME_STEP','NUMBER_ITERATIONS','MONITORING','NUMBER_SAVE','NUMBER_PRINT','NEXT_PRINT','NEXT_SAVE','NEXT_UPDATE','-v7.3');
Number_Save=0;
ExportStatic;
ExportDynamic;