function [poles] = companion2poles(B, n, dt)
% poles = companion2poles(B, n, dt) solves eigenvalue problem
% using companion matrix
%
% B       denominator coefficient
% n       at order n
% dt      sampling period
%
% poles   complex poles
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


% discrete poles
MC = [zeros(n-1, 1), eye(n-1); -B(1:n)'] ;
%MC = [zeros(n-1, 1), eye(n-1); -B(1:n)'*B(n+1)] ;
valp = eig(MC) ;

% poles
poles = log(valp)/dt ;

% stable poles with negative real part
poles = poles(real(poles) < 0) ;

% from conjugate to poles with positive imaginary part
poles = poles(imag(poles) > 0) ;
