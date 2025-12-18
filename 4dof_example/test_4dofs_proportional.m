% Mass-Spring-Damper system with 4 DOFs, clamped-free
% 
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

%digits(16)

m1 = 3 ;
m2 = 1 ;
m3 = 2 ;
m4 = 1 ;

k1 = 1 ;
k2 = 2 ;
k3 = 1 ;
k4 = 1 ;

M = [m1 0 0 0; 
0 m2 0 0; 
0 0 m3 0; 
0 0 0 m4] ;

K = [k1+k2 -k2 0 0;
-k2 k2+k3 -k3 0;
0 -k3 k3+k4 -k4;
0  0 -k4 k4];

% proportional Rayleigh damping
alpha = 1.e-2; beta = 1.e-2 ;
C = alpha*K + beta*M ;
rec = []; mob = []; acc = [] ;
for iF=1:2
    % load
    F = zeros(4, 1) ;
    if iF == 1
        F(1)=1;
    elseif iF == 2
        F(3)=1;
    end

    w = 0:0.001:3 ;
    nw = length(w) ;

    nw = length(w) ;
    REC = zeros(4, nw) ;
    ACC = zeros(4, nw) ;
    MOB = zeros(4, nw) ;

    for iw = 1:1:nw
        REC(:, iw) = linsolve((K+1j*w(iw)*C-w(iw)^2*M), F) ;
        ACC(:, iw) = linsolve((K+1j*w(iw)*C-w(iw)^2*M), -w(iw)^2*F) ;
        MOB(:, iw) = linsolve((K+1j*w(iw)*C-w(iw)^2*M), 1j*w(iw)*F) ;    
    end

    figure(1)

    REClog = 20*log10(abs(REC)) ;
    ACClog = 20*log10(abs(ACC)) ;
    MOBlog = 20*log10(abs(MOB)) ;

    subplot(1,3,1)
    plot(w/2/pi, REClog) ;
    title('Receptance')

    subplot(1,3,2)
    plot(w/2/pi, ACClog) ;
    title('Accelerance')

    subplot(1,3,3)
    plot(w/2/pi, MOBlog) ;
    title('Mobility')

    w = w.' ;
    rec = [rec, REC.'] ;
    acc = [acc, ACC.'] ;
    mob = [mob, MOB.'] ;

end



save('omg_tot.mat', 'w') ;
save('REC_tot.mat', 'rec') ;
save('ACC_tot.mat', 'acc') ;
save('MOB_tot.mat', 'mob') ;

%[V,D] = eig(K,M);
[V,D] = eig(M^-1*K);
[lambda,k] = sort(diag(D));
V = V(:,k) ;
Factor = diag(V'*M*V) ;
Vnorm = V*inv(sqrt(diag(Factor))) ;
k = Vnorm' * K*  Vnorm ;
f = diag(sqrt(k)) / (2*pi) ;
c = Vnorm' * C * Vnorm ;
xi = diag(c) ./ (2*2*pi*f) ;

