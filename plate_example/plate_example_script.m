%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors :
% Baptiste Chomette, Ecole Centrale de Lyon, LTDS, baptiste.chomette@ec-lyon.fr
% Jean-Loic Le Carrou : Sorbonne Université, d'Alembert, jean-loic.le_carrou@sorbonne-universite.fr
% Sami Karkar : Sorbonne Université, d'Alembert, sami.karkar@free.fr
% François Fabre : Sorbonne Université, d'Alembert, fabrefrancois8@gmail.com
%
% Aknowledgements :
% This work, part of the project Ngombi, was funded by the Agence Nationale de la Recherche
% (French National research agency), Grant No. ANR-19-CE27-0013-01.
%
% License :
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%
% Realease : v0 ... 2025
%
% Reference : 
% B. Chomette, J-L. Le Carrou, S. Karkar, F. Fabre, ..., JOSS, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% % make sure the relative path is correct
addpath('../MPIFD_MIMO')

% define the type of measurement: here, accelerance
frfType = 'acc' ;

% Do we want to make the mode computation using only real-valued residues?
useRealResidues = false ;

% Do we want to use the corrective terms LR/UR for out of band modes?
useLRUR = true ;

% LOAD MEASURES
% LAM rectabgular plate, measured by François Fabre
% see picture for DOF numbering
% excitation in point 35
load('FRF_plaque_moyenne.mat')
w = 2*pi*FRF_plaque_moyenne.freq ;
frf=FRF_plaque_moyenne.frf_moy;

% choose the frequency range of interest for modelling
[w1, frf1] = select_frf(w, frf, 10, 2000) ;


% List of orders to try for drawing the stabilization chart
order = 10:5:100 ;

% Description of the system : number of the DOFs with actuators and sensors
dof_a = 35 ;
dof_s = 1:45 ;
nPts = 45;

% FRF(s) to display
% DOF for actuation: 35, DOF for sensor: 35
% We choose to display the colocalised response function H[35,35]
listIOVisu = [35 35];

% Coordinates of the DOFs in a frame of reference (used for mode shape plots)
load('coord_points_plaque.mat') ;
x=reshape(x_y_plaque(:,1),5,9)';
x=x(:);
y=reshape(x_y_plaque(:,2),5,9)';
y=y(:);

COORD = [x,y] ;

% % option 1 : 2D visu
% meshStructure.dimension=2;
% meshStructure.COORD=COORD;

% option 2 : generating the structure for 3D mode shape visualization
sizeX = size(x') ;
meshStructure.model3D.vertices = [x' ; y' ; zeros(sizeX)] ;
meshStructure.model3D.vertexNormals = [zeros(sizeX) ; zeros(sizeX) ; ones(sizeX)] ;
% Manually creating faces for the 3D model using rectangles
nY=numel(unique(y));
nX=numel(unique(x));
faces=zeros(4,(nX-1)*(nY-1));
for ii=1:nY-1
    for jj=1:nX-1
        faces(:,jj+(ii-1)*(nX-1)) = [(ii-1)*nX+jj ; (ii-1)*nX+jj+1 ; ii*nX+jj+1 ; ii*nX+jj] ;
    end
end
meshStructure.model3D.faces = faces ;
meshStructure.model3D.frfModelLink = [(1:nPts)', ones(nPts,1)] ;
meshStructure.dimension = 3 ;


%% launch the GUI and start the identification

fig1 = mpifd_mimo(w1, frf1, dof_a, dof_s, order, useLRUR, useRealResidues, frfType, listIOVisu, meshStructure) ;


%% uncomment and run after using the GUI to identify and compute modes
% % % get data from the ui
% st=guidata(fig1);
% phi=st.id.phi;
% phir=st.id.phir;
% fst=st.id.fst;
% 
% %% MAC matrix barplot
% 
% MAC = mac_matrix(phir, phir) ;
