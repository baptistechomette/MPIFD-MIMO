function [rntot, LRtot, URtot, frfmtot, X] = lsfd(w, frf, fn, xin, res, form, type_frf)
% Least-Squares Frequency-Domain (LSFD) estimator for the residues using
% poles-residues form so that
%
% if form = 'fcomplex', rk is complex
%
%             k=n
%            ----- 
%            \           [  r_k       r_k*  ]
%    H(s) =   .     a(s) [------- + --------] + b(s) LR + c(s) UR
%            /           [ s-p_k     s-p_k* ]
%            -----
%             k=1
%
% if form = 'fnormal' (in the case of proportional damping), rk is purely imaginary
%
%             k=n
%            ----- 
%            \           [       -2 r_k imag(p_k)     ]
%    H(s) =   .     a(s) [----------------------------] + b(s) LR + c(s) UR
%            /           [ s^2 + 2 xi_k w_k s + w_k^2 ]
%            -----
%             k=1
% where
%
% a(s) |  b(s) | c(s) | transfer function
%-------------------------------------------
%  1   | 1/s^2 |  1   | receptance function
%  s   |  1/s  |  s   | mobility function
% s^2  |   1   | s^2  | accelerance function
%
% and p_k = -xi_k w_k + j w_k sqrt(1-xi_k^2), (*) denotes the conjugate
%
% w               angular frequency vector in rad/s
% frf             complex frequency response function matrix between
%                 f_start and f_end
%                 each column corresponds to one FRF between one sensor and one
%                 actuator
% fn              eigen frequency in Hz obtained using lsfc.m
% xin             modal damping factor obtained using lsfc.m
% res             =1 with residual terms
%                 =0 without residual terms
% form            'fnormal' for identification using real residues
%                 'fcomplex' for identification using complex residues
% type_FRF        'acc', 'vel' or 'disp' depends on the FRF measurements
%
% rntot           matrix with complex residues (one column = one FRF)
% LRtot & URtot   matrix with Lower and Upper Residual terms
% frfmtot         matrix, one column for one estimated modal frf
%
% References:
% Sandro Amador, Mahmoud El-Kafafy et al, A New Maximum Likelihood Estimator 
% Formulated in Pole-Residue Modal Model, Applied Sciences 9, 3120, 2019,
% doi:10.3390/app9153120
%
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


if ~isempty(fn)
    [L, nfrf] = size(frf) ;
    s = 1j*w ;
    nmodes = length(fn) ;
    wn = 2*pi*fn ;
    pn = -xin.*wn+1j*wn.*sqrt(1-xin.^2) ;
    pnc = conj(pn) ;

    rntot = zeros(nmodes, nfrf) ;
    LRtot = zeros(1, nfrf) ;
    URtot = zeros(1, nfrf) ;
    frfmtot = zeros(L, nfrf) ;

    % least-square problem formulated to guarantee the realness of X
    %           [X]
    % [real(C)] [real(FRF)]
    % [imag(C)] [imag(FRF)]
    %
    % X=[real(r1) imag(r1) ... real(rn) imag(rn) real(LR) imag(LR) real(UR) imag(UR)]^t
    % to obtain rk and rk* in conjugate form

    if strcmp(form, 'fcomplex')
        C = zeros(L, 2*nmodes+4) ;
        % residual terms formulated using [real(LR) imag(LR) real(UR) imag(UR)]
        if strcmp(type_frf, 'disp')
            coeff = ones(L,1) ;
            C(:,2*nmodes+1:2*nmodes+4) = ...
                res*[1./(s.^2) 1j./(s.^2) ones(L,1) 1j*ones(L,1)] ;
        elseif strcmp(type_frf, 'vel')
            coeff = s ;
            C(:,2*nmodes+1:2*nmodes+4) = res*[1./s 1j./s s 1j*s] ;
        elseif strcmp(type_frf, 'acc')
            coeff = s.^2 ;
            C(:,2*nmodes+1:2*nmodes+4) = res*[ones(L,1) 1j*ones(L,1) s.^2 1j*s.^2] ;
        end


        for k=1:2:2*nmodes-1
            C(:,k:k+1) = coeff.*[1./(s-pn((k+1)/2))+1./(s-pnc((k+1)/2)) ...
                                 1j./(s-pn((k+1)/2))-1j./(s-pnc((k+1)/2))] ;        
        end

        C = [real(C); imag(C)] ;

        for ifrf=1:nfrf
            FRF = [real(frf(:, ifrf)); imag(frf(:, ifrf))] ;
            X = lsqminnorm(C, FRF) ;
            rntot(:, ifrf) = X(1:2:2*nmodes-1) + 1j*X(2:2:2*nmodes)  ;
            LRtot(1, ifrf) = X(2*nmodes+1) + 1j*X(2*nmodes+2) ;
            URtot(1, ifrf) = X(2*nmodes+3) + 1j*X(2*nmodes+4) ;
            frfmtot(:, ifrf) = C(1:L,:)*X + 1j*C(L+1:2*L,:)*X ;
        end

    elseif strcmp(form, 'fnormal')
        C = zeros(L, nmodes+4) ;
        % residual terms formulated using [real(LR) imag(LR) real(UR) imag(UR)]
        if strcmp(type_frf, 'disp')
            coeff = ones(L,1) ;
            C(:,nmodes+1:nmodes+4) = ...
                res*[1./(s.^2) 1j./(s.^2) ones(L,1) 1j*ones(L,1)] ;
        elseif strcmp(type_frf, 'vel')
            coeff = s ;
            C(:,nmodes+1:nmodes+4) = res*[1./s 1j./s s 1j*s] ;
        elseif strcmp(type_frf, 'acc')
            coeff = s.^2 ;
            C(:,nmodes+1:nmodes+4) = res*[ones(L,1) 1j*ones(L,1) s.^2 1j*s.^2] ;
        end    

        for k=1:nmodes
            C(:,k) = coeff.*(1./(s.^2 + 2*xin(k)*wn(k)*s + wn(k)^2)) ;       
        end

        C = [real(C); imag(C)] ;

        for ifrf=1:nfrf
            FRF = [real(frf(:, ifrf)); imag(frf(:, ifrf))] ;
            X = lsqminnorm(C, FRF) ;
            rntot(:, ifrf) = 1j*X(1:nmodes)./(-2*imag(pn))  ; 
            %rntot(:, ifrf) = X(1:nmodes)./(-2*imag(pn))  ;        
            LRtot(1, ifrf) = X(nmodes+1) + 1j*X(nmodes+2) ;
            URtot(1, ifrf) = X(nmodes+3) + 1j*X(nmodes+4) ;
            frfmtot(:, ifrf) = C(1:L,:)*X + 1j*C(L+1:2*L,:)*X ;
        end
    end
else
    disp('no stable poles')
end
