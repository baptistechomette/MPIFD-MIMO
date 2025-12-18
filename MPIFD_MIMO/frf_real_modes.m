function [Hnum] = frf_real_modes(phi, fp, xip, w, LR, UR, vectAct, frfType )
%FRF calculation using real modes normalized using unitary mass matrix
% vectAct   vector of actuator ddl, 1 x Ni
% phi       normalized (unitary mass matrix) real mode, No x Ni
% fp        eigen frequency, Nm x 1
% xip       modal damping factor, Nm x 1
% LR        lower residual matrix, real, No x Ni
% UR        upper residual matrix, real, No x Ni
% w         frequency vector (rad/s), Nf x 1
% frfType  'acc', 'vel' or 'disp' depends on the FRF measurements
%
% Hnum      complex FRF matrix, Nf x (No x Ni)
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

    Ni = length(vectAct) ;
    [No, Nm] = size(phi) ;
    Nf = length(w) ;
    wp = 2*pi*fp ;
    s = 1j*w ;

    pk = -xip.*wp+1j*wp.*(1-xip.^2).^(1/2) ;
    
    Hnum = zeros(Nf, No*Ni) ;
    for ki = 1:Ni
        for ko = 1:No
            Hm = zeros(Nf, 1) ;
            for kr = 1:Nm
                Hm = Hm + ...
                     phi(ko, kr) * phi(vectAct(ki), kr) ./ ...
                     (wp(kr)^2 - w.^2 + 1j*2*xip(kr)*wp(kr)*w) ;
            end
            % FRF in displacement / force formulation
            if strcmp(frfType, 'disp')
                Hnum(:, (ki-1)*No+ko) = Hm + LR((ki-1)*No+ko)./(s.^2) + UR((ki-1)*No+ko) ;
            elseif strcmp(frfType, 'vel')
                Hnum(:, (ki-1)*No+ko) = s.*Hm + LR((ki-1)*No+ko)./s + s.*UR((ki-1)*No+ko) ;
            elseif strcmp(frfType, 'acc')
                Hnum(:, (ki-1)*No+ko) = s.^2.*Hm + LR((ki-1)*No+ko) + s.^2.*UR((ki-1)*No+ko) ;
            end   
        end
    end
end
