function [freq, fftdata] = spectrum(data, fs)
% [freq, fftdata] = spectrum(data, fs)
% Returns frequency vector (freq) and complex spectrum vector (fftdata)
%
% data   time data vector
% fs     sample frequency in Hz
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

ts = 1 / fs ;

L = length(data) ;

nfft = 2^nextpow2(L) ;

fftdata = fft(data, nfft) / L ; % normalization

fftdata = fftdata(1: nfft / 2) ;

% multiply by 2 to take into account the fact that we threw out second half
fftdata = fftdata * 2 ;
freq = fs / 2. * linspace(0., 1., nfft / 2) ;

