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
% GNU CC BY-NC-SA
%
% Realease : v0 ... 2025
%
% Reference : 
% B. Chomette, J-L. Le Carrou, S. Karkar, F. Fabre, ..., JOSS, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

addpath('../MPIFD_MIMO')

% file test using 4 dofs system with proportional damping
% for MPIFD_MIMO algorithm
%
% modal parameters of the analytical 4 dofs system
% m1=3, m2=1, m3=2, m4=1
% k1=1, k2=2, k3=1, k4=1
% c1=c2=c3=c4=0.1
%
% 4 FRF with excitation located on dof 1 and 4 FRF with excitation located en dof 3 
%
%  -- f --   - xi -
%   0.0494   0.0177
%   0.1203   0.0104
%   0.2037   0.0103
%   0.3059   0.0122

% phi using norm2 so that norm(phi(:,i)) = 1
% 0.2992    0.6152   -0.1589   -0.2352
% 0.4055    0.3959    0.1521    0.9507
% 0.5791   -0.2688    0.5249   -0.1895
% 0.6409   -0.6265   -0.8223    0.0703

% scaled modes using mass matrix (unitary modal mass)
phi_an = [0.2431 0.4549 -0.1380 -0.2197;
 0.3295 0.2927 0.1321 0.8878 ;
 0.4705 -0.1987 0.4558 -0.1770 ;
 0.5208  -0.4632 -0.7141 0.0657] ;

% data type: disp, vel ou acc
test = 'disp' ;

% Use complex or Real-valued residues
useRealResidues = false ;

% use residual terms for LF or HF contributions of unmodelled modes?
useLRUR = false ;

% load frequency response function in matrix form [3001x(NoxNi)]
if strcmp(test, 'disp')
    load('REC_tot.mat') ;
    frf = rec ;
    frfType = 'disp' ;
elseif strcmp(test, 'vel')
    load('MOB_tot.mat') ;
    frf = mob ;
    frfType = 'vel' ;
elseif strcmp(test, 'acc')
    load('ACC_tot.mat') ;
    frf = acc ;
    frfType = 'acc' ;
end

% load angular frequency in rad/s and in vector form [3001x1]
load('omg_tot.mat')

% selection of the FRF range using frequency betwen 0.02 and 0.38 Hz
% [w1, frf1] = select_frf(w, frf, 0.02,  0.38) ;
[w1, frf1] = select_frf(w, frf, 0.086,  0.25) ;

% description of the system inupts and outputs
vectAct = [1 3] ; % driving points: actuators placed on ddl num 1 and 3
vectSens = [1 2 3 4] ; % sensors placed on ddl 1,2,3 and 4

% identification using stabilization chart between order 2 and 20
vectOrder = 2:20 ;

% listIOVisu = [1 1 ; 3 1 ; 1 3 ; 2 3] ; % we want to display the FRF of transfer 
% functions H11, H31, H13, and H23

%% start the Identification UI
fig1 = mpifd_mimo(w1, frf1, vectAct, vectSens, vectOrder, useLRUR, ...
    useRealResidues, frfType) ;



%% after completion of steps in the GUI, uncomment and run the following
% %% mass verification (M=diag(3,1,2,1))
% st=guidata(fig1);
% phir=st.id.phir;
% M = (phir*phir')^-1 ;
% figure(2)
% bar3((M))
% 
% %% mac matrix verification
% phi=st.id.phi;
% MAC = mac_matrix(phi, phi_an) ;
