clc
close all
clear

%% ========================INPUT PARAMETERS=========================
%===Geometrical Parameters
    % Average Length of filaments
    L_Filament = 0.4;
    % Proportion of each filament = total length of the filements / cell area 
    Area_Ratio_Major = 4;
    Area_Ratio_Minor = 4;
    %Average Size of mesh
    Avg_Mesh_size = 1.5;
    %Radius to search for elements with Max or Min strains
    Search_Rad_Maj=9;
    Search_Rad_Min=9;
%=========================

%====Mechanical Parameters
    %Young's modulus of cytoplasm
    E_Cell=0.5; %kPa
    Nu_Cell=0.3;
    %Young's modulus of nucleus
    E_Nucl=1.5; %kPa
    Nu_Nucl=0.3;
    %Properties of Fibers
    E_Fiber=1.5;%kPa
    A_Fiber=0.03;%um^2
    %Membrane Properties
    E_Membrane = 100000;%kPa
    A_Membrane=0.004;%um^2
    %Nuclear Membrane Properties
    E_Nuc_Mem=100000;%kPa
    A_Nuc_Mem=0.004;%um^2
    %Spring Suport Stiffness
    K_Support=20000;%pN/um
%=========================

%====Filament Dynamics Parameters
    %Max stress in filaments
    Sigma_Max=2.5; %kPa (nN/um^2) 
    %decay constant of signal
    Theta=300;%sec
    %constant controlling rate of formation
    kf=10;
    %constant controlling rate of dissociation
    kb=1;
    %Shortening Rate
    Eps0_dot=2.8e-3; %1/sec
    %Constant controlling fractional redution in fiber stress
    kv=10;
    %Signalling Process
    c_signal=0;
    %Parameters adjusting the rate of filament dislocation/reorientation
    Disp_Factor_Maj = 0.2;
    Rotat_Factor_Maj=0.2;
    Disp_Factor_Min = 0.2;
    Rotat_Factor_Min=0.2;
    %Time step
    Time_Step=10;%sec
    %Number of steps
    N_Step=100;
%=================================================================
save('Parameters.mat',"Avg_Mesh_size","Search_Rad_Maj","L_Filament","Area_Ratio_Major","Area_Ratio_Minor",...
    "Disp_Factor_Maj","Rotat_Factor_Maj","Disp_Factor_Min","Rotat_Factor_Min","E_Cell","Nu_Cell","E_Nucl",...
    "Nu_Nucl","E_Fiber","A_Fiber","E_Membrane","A_Membrane","E_Nuc_Mem","A_Nuc_Mem","K_Support",...
    "Sigma_Max","Theta","kf","kb","Eps0_dot","kv","c_signal","Time_Step","N_Step","Search_Rad_Min");

%% To read coordinates of boundary, nucleus, and supports from text files
[Boundary_Coord,Nucleus_Coord,Support_Coord] = Read_Coord('Boundary_Coorrdinates.txt','Nucleus_Coordinates.txt','Support_Coordinates.txt');

%To generate mesh in the cell (cytoplasmic) and its nucleus.
model=Mesh_Generation(Boundary_Coord,Nucleus_Coord, Avg_Mesh_size);
save('FE_Data.mat','model','Boundary_Coord','Nucleus_Coord','Support_Coord');

% To generate Major and Minor fibres
[Major_Points_X,Major_Points_Y,Minor_Points_X,Minor_Points_Y] = Fibre_Generation(Boundary_Coord,Nucleus_Coord);
save('Fillement_Data.mat','Major_Points_X','Major_Points_Y','Minor_Points_X','Minor_Points_Y');

%% Finite Element Analysis
FE_Model();


