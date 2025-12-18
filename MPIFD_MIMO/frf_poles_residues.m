function [Hnum] = frf_poles_residues(PSI, fp, xip, L, w, LR, UR)
%FRF calculation using complex modes and modal participation factor
% PSI    complex mode before normalization, No x Nm x Ni
% fp     eigen frequency, Nm x 1
% xip    modal damping factor, Nm x 1
% L      modal participation factor, Ni x Nm
% w      frequency vector (rad/s), Nf x 1
%
% Hnum   complex FRF matrix, Nf x No x Ni
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

    [No, Nm, Ni] = size(PSI) ;
    Nf = length(w) ;
    wp = 2*pi*fp ;
    poles = -xip.*wp + 1j*wp.*sqrt(1-xip.^2) ;
    poles_c  = conj(poles) ;
    L_c = conj(L) ;
    s = 1j*w ;
    PSI_c = conj(PSI) ;
    Hnum = zeros(Nf, No, Ni) ;
    for ki=1:Ni
        for ko=1:No
            Hm = zeros(Nf, 1) ;
            for kr=1:Nm
                Hm = Hm + PSI(ko, kr, ki)*L(ki, kr)./(s-poles(kr)) + ...
                     PSI_c(ko, kr, ki)*L_c(ki, kr)./(s-poles_c(kr)) ;
            end
            Hnum(:, ko, ki) = Hm + LR(ko, ki)./(s.^2) + UR(ko, ki) ;
        end
    end
end
