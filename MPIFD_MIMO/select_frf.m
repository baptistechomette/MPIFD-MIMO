function [w1, frf1] = select_frf(w0, frf0, f_start, f_end)
% [w1, frf1] = select_frf(w0, frf0, f_start, f_end)
% selects FRF between f_start and f_end
%
% w0       natural frequency vector in rad/s
% frf0     complex frequency response function matrix
% f_start  start frequency
% f_end    end frequency
%
% w1      natural frequency vector in rad/s between w_start = 2*pi*f_start and
%         w_end = 2*pi*f_end
% frf1    complex frequency response function matrix between f_start and f_end
%         each column corresponds to one FRF between one sensor and one actuator
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
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%
% Realease : v0 ... 2025
%
% Reference : 
% B. Chomette, J-L. Le Carrou, S. Karkar, F. Fabre, ..., JOSS, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f0 = w0/(2*pi) ;
[diff_start, i_start] = min(abs(f0-f_start)) ;
[diff_end, i_end] = min(abs(f0-f_end)) ;
w1 = w0(i_start:i_end) ;
frf1 = frf0(i_start:i_end, :) ;

