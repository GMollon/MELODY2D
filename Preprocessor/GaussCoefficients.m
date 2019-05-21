function GaussPW=GaussCoefficients(Nperdimension)
if Nperdimension==1
    GaussPW=[0;2];
elseif Nperdimension==2
    GaussPW=[-0.57735,0.57735;1,1];
elseif Nperdimension==3
    GaussPW=[-0.77459,0,0.77459;0.55555,0.88888,0.55555];
elseif Nperdimension==4
    GaussPW=[-0.86113,-0.33998,0.33998,0.86113;0.34785,0.65214,0.65214,0.34785];
elseif Nperdimension==5
    GaussPW=[-0.93246,-0.66120,-0.23861,0.23861,0.66120,0.93246;...
        0.17132,0.36076,0.46791,0.46791,0.36076,0.17132];
elseif Nperdimension>=6
    GaussPW=[-0.96028,-0.79666,-0.52553,-0.18343,0.18343,0.52553,0.79666,0.96028;...
    0.10122,0.22238,0.31370,0.36268,0.36268,0.31370,0.22238,0.10122];
end