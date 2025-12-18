function [MAC] = mac_matrix(phi1, phi2, varargin)
% [MAC] = mac_matrix(phi1,phi2)
% 
% Computes and displays the MAC Matrix between two sets of modes as a 3D
% bar plot
% 
% -- Input Arguments --
% phi1: first set of column vectors [Nddl x Nmodes]
% phi2: second set of column vectors [Nddl x Nmodes]
% [hfig] : (optional) handle to a Matlab axes where to plot
% 
% -- Output Arguments --
% [MAC] : (optional) MAC matrix
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

    [ndof, nmodes] = size(phi1) ;
    MAC = zeros(nmodes, nmodes) ;
    for i=1:nmodes
        for ii=1:nmodes
            MAC(i, ii) = mac_mpifd(phi1(:, i), phi2(:, ii)) ;
        end
    end
    
    if ~isempty(varargin) && ishandle(varargin{1})
        myax = varargin{1};
    else
        hfig= figure();
        myax=axes();
    end

    b = bar3(myax, MAC) ;
    % proportionnal color to heigth
    for k = 1:nmodes
        zdata = b(k).ZData ;
        b(k).CData = zdata ;
        b(k).FaceColor = 'interp' ;
    end
    title(myax, 'MAC matrix')
    xlabel(myax, 'phi 1')
    ylabel(myax, 'phi 2')
end

