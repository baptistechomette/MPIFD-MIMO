function [fst, xist] = stabchart(fp, xip, fmin1, ximin1, f, fst, xist)
% [fst, xist] = stabchart(fp, xip, fmin1, ximin1, yfrf, ip, ...
% f, fst, xist, ff, xixi, mathp, istab, p)
% Stabilization chart to extract stable poles in frequency and damping
%
% fp and xip        frequency and damping identified before poles selection
%                   at iteration p
% fmin1 and ximin1  frequency and damping identified at iteration p-1
% yfrf              amplitude in dB, yfrf = 20*log10(abs(frf))
% ip                iteration index
% f                 frequency vector in Hz
% fst               frequency of stables poles in frequency and damping (s)
% xist              damping of stable poles in frequency and damping (s)
% ff                frequency of stable poles in frequency (f)
% xixi              frequency of stable poles in damping (d)
% mathp             frequency of mathematical poles (o)
% p                 identification order
%
% fst and xist      frequency in Hz and modal damping factor of stable poles in
%                   frequency and damping
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

tolf = 0.01 ;
told = 0.05 ;

fmax = max(f) ;
fmin = min(f) ;
nfp = length(fp) ;
if length(fmin1) > 0
    for ifp = 1:nfp
        if (fp(ifp) < fmax) && (fp(ifp) > fmin)
            [diff_fp, idiff_fp] = min(abs(fp(ifp)-fmin1) / fp(ifp)) ;
            diff_xip = abs(xip(ifp)-ximin1(idiff_fp)) / xip(ifp) ;
            
            if (diff_fp <= tolf) && (diff_xip <= told)
                fst = [fst; fp(ifp)] ;
                xist = [xist xip(ifp)] ;
            end
        end
    end
end
end

