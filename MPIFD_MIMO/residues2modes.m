function[phir, phi, residues] = residues2modes(rntot, vectAct, vectSens, fn, xin, useRealResidues)
% returns mode shapes and modal participation factor using complex residues and
% decomposition in singular values
%
% vectAct     vector: ddl for actuator [ddl1, ddl2, ...]
% vectSens     vector: ddl for sensor [ddl1, ddl2, ...]
% rntot    complex residues: nb modes x nb frf, one column = one FRF
% fns      stable frequency
% xins     stable damping
%
% phir     normalized real mode shapes using mass normalization
% phi      complex mode using mass normalization
% residues 3D matrix [nb sensors x nb actuators x nb modes]
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

phi = [] ;
phir = [] ;

[nb_modes, noni] = size(rntot) ;
nba = length(vectAct) ;
nbs = length(vectSens) ;

wn = 2*pi*fn ;

residues = zeros(nbs, nba) ;

if ~useRealResidues
    for k = 1:nb_modes
        for ia = 1:nba
            residues(:, ia, k) = rntot(k, (ia-1)*nbs+1:ia*nbs) ;
            [u, s, v] = svd(residues(:, :, k)) ;
            phik = u(:, 1) ;
        end
        % normalization
        pp = (phik*phik.') ;
        % reshape in column
        pp = reshape(pp(:, vectAct), [noni, 1]) ;
        rr = reshape(residues(:, :, k), [noni, 1]) ;
        muk = lsqminnorm(rr, pp) ;
        im_pk = wn(k)*sqrt(1-xin(k)^2) ;
        % normalization using unitary modal mass
        phink = phik / sqrt(muk) ;
        % real mode simplification (phink*phink must be purely imaginary)
        phink2real = phink*sqrt(2j*im_pk) ;
        if imag(muk) < 0
            phirk = real(phink2real) ;
        elseif imag(muk) > 0
            phirk = imag(phink2real) ;
        end
        % save complex mode with unitary modal mass
        phi = [phi, phink] ;
        % save real mode with unitary modal mass
        phir = [phir, phirk] ;
    end

else
    for k = 1:nb_modes
        for ia = 1:nba
            residues(:, ia, k) = rntot(k, (ia-1)*nbs+1:ia*nbs) ;
            [u, s, v] = svd(residues(:, :, k)) ;
            phik = imag(u(:, 1)) ;
        end
        % normalization
        pp = (phik*phik.') ;
        % reshape in column
        pp = reshape(pp(:, vectAct), [noni, 1]) ;
        rr = reshape(residues(:, :, k), [noni, 1]) ;
        im_pk = wn(k)*sqrt(1-xin(k)^2) ;
        muk = lsqminnorm(-2*imag(rr)*im_pk, pp) ;
        % normalization using unitary modal mass
        phink = phik / sqrt(muk) ;
        if muk > 0
            phirk = real(phink) ;
        elseif muk < 0
            phirk = imag(phink) ;
        end
        % save complex mode with unitary modal mass
        phi = [phi, phink] ;
        % save real mode with unitary modal mass
        phir = [phir, phirk] ;
    end
end
