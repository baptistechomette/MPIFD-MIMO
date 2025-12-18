function [w, frf] = time2frf(file_name)
% [w, frf] = time2frf(file_name) 
% returns natural frequency vector w in rad/s  and complex frequency response
% function frf
%
% file_name: 'data.txt'
% no title line
% -----------------------------------------------------------------
%   first column |    second column     |        third column
% -----------------------------------------------------------------
%     time in s  ,         input        ,            output
% -----------------------------------------------------------------
%
% delimiter: ','
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


fid = fopen(file_name) ;
format = '%f %f %f' ;
delimiter = ',' ;

A = textscan(fid, format, 'delimiter', delimiter) ;

time = A{:, 1} ;
input = A{:, 2} ;
output = A{:, 3} ;

ts = time(2)-time(1) ;
fs = 1/ts ;

[freq, fft_input] = spectrum(input, fs) ;
[freq, fft_output] = spectrum(output, fs) ;

w = 2*pi*freq ;
w = w.' ;
frf = fft_output./fft_input ;
