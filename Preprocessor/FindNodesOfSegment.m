function [X0,Y0,Node0,X1,Y1,Node1,X2,Y2,Node2,X3,Y3,Node3]=FindNodesOfSegment(Body,Border,Shift,Segment)
global BODIES_STATIC BODIES_CURRENT PERIODIC_BOUNDARIES
Border_Type=BODIES_STATIC(Body).BORDERS{Border,1};
Number_Nodes=BODIES_STATIC(Body).BORDERS{Border,2};
Period=PERIODIC_BOUNDARIES(2)-PERIODIC_BOUNDARIES(1);
if Segment==1
    if strcmp(Border_Type,'Closed')
        B0=Border;
        N0=Number_Nodes;
        S0=0;
    elseif strcmp(Border_Type,'Simple')
        if Border==1
            B0=BODIES_STATIC(Body).NUMBER_BORDERS;
            N0=BODIES_STATIC(Body).BORDERS{B0,2};
        else
            B0=Border-1;
            N0=BODIES_STATIC(Body).BORDERS{B0,2};
        end
        S0=0;
    elseif strcmp(Border_Type,'Periodic')
        B0=Border;
        N0=Number_Nodes;
        if Border==1
            S0=-1;
        elseif Border==2
            S0=1;
        end
    end
    B1=Border;
    N1=Segment;
    S1=0;
    B2=Border;
    N2=Segment+1;
    S2=0;
    B3=Border;
    N3=Segment+2;
    S3=0;
elseif Segment==Number_Nodes-1
    B0=Border;
    N0=Segment-1;
    S0=0;
    B1=Border;
    N1=Segment;
    S1=0;
    B2=Border;
    N2=Segment+1;
    S2=0;
    if strcmp(Border_Type,'Closed')
        B3=Border;
        N3=1;
        S3=0;
    elseif strcmp(Border_Type,'Simple')
        if Border==BODIES_STATIC(Body).NUMBER_BORDERS
            B3=1;
        else
            B3=Border+1;
        end
        N3=1;
        S3=0;
    elseif strcmp(Border_Type,'Periodic')
        B3=Border;
        N3=1;
        if Border==1
            S3=1;
        elseif Border==2
            S3=-1;
        end
    end
elseif Segment==Number_Nodes
    B0=Border;
    N0=Segment-1;
    S0=0;
    B1=Border;
    N1=Segment;
    S1=0;
    if strcmp(Border_Type,'Closed')
        B2=Border;
        N2=1;
        S2=0;
        B3=Border;
        N3=2;
        S3=0;
    elseif strcmp(Border_Type,'Simple')
        if Border==BODIES_STATIC(Body).NUMBER_BORDERS
            B2=1;
            B3=1;
        else
            B2=Border+1;
            B3=Border+1;
        end
        N2=1;
        N3=2;
        S2=0;
        S3=0;
    elseif strcmp(Border_Type,'Periodic')
        B2=Border;
        N2=1;
        B3=Border;
        N3=2;
        if Border==1
            S2=1;
            S3=1;
        elseif Border==2
            S2=-1;
            S3=-1;
        end
    end
else
    B0=Border;
    N0=Segment-1;
    S0=0;
    B1=Border;
    N1=Segment;
    S1=0;
    B2=Border;
    N2=Segment+1;
    S2=0;
    B3=Border;
    N3=Segment+2;
    S3=0;
end
Node0=BODIES_STATIC(Body).BORDERS{B0,3}(N0);
X0=BODIES_CURRENT(Body).CURRENT_POSITIONS(Node0,1)+(S0+Shift)*Period;
Y0=BODIES_CURRENT(Body).CURRENT_POSITIONS(Node0,2);
Node1=BODIES_STATIC(Body).BORDERS{B1,3}(N1);
X1=BODIES_CURRENT(Body).CURRENT_POSITIONS(Node1,1)+(S1+Shift)*Period;
Y1=BODIES_CURRENT(Body).CURRENT_POSITIONS(Node1,2);
Node2=BODIES_STATIC(Body).BORDERS{B2,3}(N2);
X2=BODIES_CURRENT(Body).CURRENT_POSITIONS(Node2,1)+(S2+Shift)*Period;
Y2=BODIES_CURRENT(Body).CURRENT_POSITIONS(Node2,2);
Node3=BODIES_STATIC(Body).BORDERS{B3,3}(N3);
X3=BODIES_CURRENT(Body).CURRENT_POSITIONS(Node3,1)+(S3+Shift)*Period;
Y3=BODIES_CURRENT(Body).CURRENT_POSITIONS(Node3,2);
