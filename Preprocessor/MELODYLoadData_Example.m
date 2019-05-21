% 1. Simulation name
Simulation_Name='Example';


% 2. Initilization
Contours=cell(1,1);
Distributions=cell(1,1);
Interpolations=cell(1,1);
Integrations=cell(1,1);
Bodies_Materials=cell(1,1);
Imposed_Pressures=cell(1,1);
Imposed_Velocities=cell(1,1);
Alid=cell(1,1);
Deactivations=cell(1,1);
To_Plot=zeros(27,1);


% 3. Bodies
Contours{1,1}={'Closed',cat(2,-1+0.5*cos([0:1:360]'/180*pi),1+0.5*sin([0:1:360]'/180*pi)),'Linear'};
Distributions{1,1}={'Unstructured'};
Distributions{1,2}=0.05;
Distributions{1,3}=8;
Interpolations{1,1}='MLS';
Interpolations{1,2}=10;
Integrations{1,1}='Gauss';
Integrations{1,2}=3;
Detections(1,1)=0.01;
Detections(1,2)=0.01;
Bodies_Materials{1,1}='Material1';
Imposed_Pressures{1,1}={[],'None',[],'None';...
                        };
Imposed_Velocities{1,1}={[],'None',[],[],'None',[];...
                         };
Initial_Velocities{1,1}=4*rand(1,2)-[2,2];
Mesh_Ratios(1,1:2)=[3,1];
Status{1,1}='active';
Alid{1,1}=[];
Alid{1,2}=[];
Alid{1,3}=[];

Contours{2,1}={'Closed',cat(2,1+0.5*cos([0:1:360]'/180*pi),1+0.5*sin([0:1:360]'/180*pi)),'Linear'};
Distributions{2,1}={'Rigid'};
Distributions{2,2}=0.05;
Distributions{2,3}=8;
Interpolations{2,1}='None';
Interpolations{2,2}=10;
Integrations{2,1}='Gauss';
Integrations{2,2}=3;
Detections(2,1)=0.01;
Detections(2,2)=0.01;
Bodies_Materials{2,1}='Material2';
Imposed_Pressures{2,1}={[],'None',[],'None';...
                        };
Imposed_Velocities{2,1}={[],'None',[],[],'None',[];...
                         };
Initial_Velocities{2,1}=4*rand(1,2)-[2,2];
Mesh_Ratios(2,1:2)=[1,1];
Status{2,1}='active';
Alid{2,1}=[];
Alid{2,2}=[];
Alid{2,3}=[];

Contours{3,1}={'Periodic',[-2,-1;2,-1],'Linear';'Periodic',[2,0;-2,0],'Linear'};
Distributions{3,1}={'Unstructured'};
Distributions{3,2}=0.05;
Distributions{3,3}=8;
Interpolations{3,1}='MLS';
Interpolations{3,2}=10;
Integrations{3,1}='Gauss';
Integrations{3,2}=3;
Detections(3,1)=0.01;
Detections(3,2)=0.01;
Bodies_Materials{3,1}='Material1';
Imposed_Pressures{3,1}={[],'None',[],'None';...
                        [],'None',[],'None';...
                        };
Imposed_Velocities{3,1}={[0,0;1e6,0],'Soft',[1e6,0],[0,0;1e6,0],'Soft',[1e6,0];...
                         [],'None',[],[],'None',[];...
                         };
Initial_Velocities{3,1}=[0,0];
Mesh_Ratios(3,1:2)=[4,1];
Status{3,1}='active';
Alid{3,1}=[];
Alid{3,2}=[];
Alid{3,3}=[];


% 4. Materials and contact laws
Materials={'Material1','NeoHookean',[1,0.001,0.001,500,0.3];...
           'Material2','NeoHookean',[1,0,0,0,0];...
           };
Contact_Laws={'Material1','Material1','MohrCoulomb','Evolutive',[1e6,1e6,0.5,0,0];...
              'Material1','Material2','MohrCoulomb','Evolutive',[1e6,1e6,0.5,0,0];...
              };


% 5. General boundary conditions
Periodic_Boundaries=[-2,2];
Gravity=[0;-10];


% 6. Text outputs
Monitorings=cell(1,1);
Spies=cell(1,1);


% 7. Graphic outputs
To_Plot(1)=1;   % Body Index
To_Plot(2)=1;   % Initial Position
To_Plot(3)=0;   % Current Position
To_Plot(4)=1;   % Displacement
To_Plot(5)=1;   % Velocity
To_Plot(6)=0;   % Acceleration
To_Plot(7)=0;   % Force
To_Plot(8)=0;   % Internal Force
To_Plot(9)=0;   % Contact Force
To_Plot(10)=0;  % Body Force
To_Plot(11)=0;  % Dirichlet Force
To_Plot(12)=0;  % Neumann Force
To_Plot(13)=0;  % Damping Force
To_Plot(14)=0;  % Alid Force
To_Plot(15)=0;  % Jacobian
To_Plot(16)=1;  % Cauchy XX Stress
To_Plot(17)=1;  % Cauchy YY Stress
To_Plot(18)=1;  % Cauchy XY Stress
To_Plot(19)=0;  % Cauchy ZZ Stress
To_Plot(20)=0;  % Tresca Stress
To_Plot(21)=1;  % Von Mises Stress
To_Plot(22)=0;  % Major Principal Stress
To_Plot(23)=0;  % Intermediate Principal Stress
To_Plot(24)=0;  % Minor Principal Stress
To_Plot(25)=0;  % Spherical Stress
To_Plot(26)=0;  % Green-Lagrange XX strain
To_Plot(27)=0;  % Green-Lagrange YY strain
To_Plot(28)=0;  % Green-Lagrange XY strain
To_Plot(29)=0;  % Norm of the Green-Lagrange strain tensor
To_Plot(30)=1;  % Body Damage
To_Plot(31)=1;  % Normalized Displacement Error
To_Plot(32)=0;  % Displacement Error


% 8. Numerical parameters
Scheme='Adaptive_Euler';
Scheme_Parameters=[1e-4,0.2,2];
Contact_Updating_Period=0.001;
Time_Stepping_Parameters=[0,1e-4,5];
Save_Periods=[0.01,0.01];


% 9. Flags
Activate_Plot=1;
Initialize_CZM=0;