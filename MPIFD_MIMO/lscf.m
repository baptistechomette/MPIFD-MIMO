function [POLES] = lscf(w, frf, n)
% Linear Square Complex Frequency estimator using discrete-time z-model
% version: 31/01/2024
%    
%             k=n
%            ----- 
%            \         k  
%            .     A z
%            /      
%            -----
%             k=0                              j k w dt
% H(z) = -------------------------,    z^k = e
%             k=n
%            ----- 
%            \         k    
%            .     B z
%            /      
%            -----
%             k=0
%
% the frequency axis between f0 and fend is shifted using
%
%             1
% dt = --------------
%       2 (fend - f0)
%
% w = 2 pi (f - f0)
%
% w          natural frequency vector in rad/s
% frf        complex frequency response function matrix
%            each column corresponds to one FRF between one sensor
%            and one actuator
% n          [nmin : nmax] identification using stabilization chart
%
% fn         eigen frequency in Hz
% xin        modal damping factor
% frfnumtot  matrix with numerical identified complex frequency response
%            function using discrete-time z-model
% FST        cell array with frequency of stable poles in frequency and damping
% FF         cell array with frequency of stable poles in frequency
% XIXI       cell array with frequency of stable poles in damping
% MATHP      cell array with frequency of mathematical poles
%
% References:
% H. Van der Auweraer, P. Guillaume, P. Verboven and S. Vanlanduit, Application 
% of a  Fast-Stabilizing Frequency Domain Parameter Estimation Method, 
%Journal of Dynamics Systems, Measurement and Control, 123, pp 651-658, 2001.
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


%% Parameters and frequency shift
w0 = w ; % natural frequency in rad/s
L = length(w) ; % length of the signal
f = w/(2*pi) ; % frequency in Hz
yfrftot = 20*log10(abs(frf)) ;
f0 = min(f) ;
fend = max(f) ;
w = 2*pi*(f-f0) ; % shifted natural frequency
dt = 1/(2*(fend-f0)) ; % sample time
[L , noni] = size(frf) ; % noni: number output * number input

frfnumtot = zeros(L, noni) ;

disp('-----------------------------------------')
disp(['Identification between order ',num2str(n(1)),' and ',num2str(n(end))])
disp('-----------------------------------------')
n_min = min(n) ;
n_max = max(n) ;
Ptot = zeros(L, n_max+1) ;
for k=1:n_max+1
    % matrix P for normal equations using z-domain model
    % dimension L * n_max+1
    Ptot(:, k) = zdomain(w, dt, k-1) ;
end

ip = 0 ;
iFST = 1 ;

%% fast implementation of normal equations using fft
% fft are padded with zeros because 2*L > L and dt=1/(2(fend-f0))
% warning: for z = exp(j w dt)
X2 = fft(conj(frf), 2*L) ;
X1 = fft(frf, 2*L) ;

% correction to avoid warning on the first term of the toeplitz matrices
X2(1,:) = X1(1,:) ;
Y1 = fft(ones(L,1), 2*L) ;
Z1 = fft(abs(frf).^2, 2*L) ;

%% fast implementation using toeplitz matrix and constant denominator order
% at maximal value
Xtot = {} ;
Ytot = {} ;
Ztot = {} ;
Atot = {} ;
M = zeros(n_max+1, n_max+1) ;
for ifrf = 1:noni
    X = toeplitz(-real(X1(1:n_max+1, ifrf)), -real(X2(1:n_max+1, ifrf).')) ;
    Y = toeplitz(real(Y1(1:n_max+1))) ;               
    Z = toeplitz(real(Z1(1:n_max+1, ifrf))) ;
    Xtot{ifrf} = X ; Ytot{ifrf} = Y ; Ztot{ifrf} = Z ;
    Mifrf = Z - X'*Y^(-1)*X ;
    M = M + Mifrf ;
end
M = 2*M ;

%% calcul between n_min and n_max
P = Ptot(:, 1:n_min) ;
h = waitbar(0, 'Identification in process') ;
for p = n_min:n_max
    waitbar(p/n_max, h) ;
    fn = [] ;
    xin = [] ;
    % concatenate in column (2)
    P = cat(2, P, Ptot(:,p+1)) ;

    %B = lsqminnorm(-M(1:p, 1:p), M(1:p, p+1)) ;
    B = lsqminnorm(M(1:p, 1:p), -M(1:p, p+1)) ;
    
    B = [B; 1] ;
    % calcul of frfnumtot only for the last iteration
    if p == n_max
        for ifrf = 1:1:noni
            % YA = -XB
            Atot{ifrf} = -lsqminnorm(Ytot{ifrf},(Xtot{ifrf}*B)) ;
            frfnumtot(:,ifrf) = (P*Atot{ifrf})./(P*B) ;
        end
    end
    % solve eigenvalue problem using companion matrix
    poles = companion2poles(B, p, dt) ;
    % modal parameters
    wp = abs(poles) ;
    % correction of the shift frequency
    wp = wp+2*pi*f0 ;
    xip = -real(poles)./wp ;
    % reorganize modal parameters
    [wp, iwp] = sort(wp) ;
    xip = xip(iwp) ;
    fp = wp/(2*pi) ;
    % stabilization chart
    if ip == 0
        fmin1 = fp ;
        ximin1 = xip ;
        
    elseif ip > 0
        [fn, xin] = stabchart(fp, xip, fmin1, ximin1, f, fn, xin) ;
        xin = xin' ;
        fmin1 = fp ;
        ximin1 = xip ;
        
        % sauvegarde sous la forme de structure
        if ~isempty(fn)
            POLES(iFST).order = p ;
            POLES(iFST).freq = fn ;
            POLES(iFST).xi = xin ;
            iFST = iFST+1 ;
        end
    end
    ip = ip+1 ;
end
close(h) ;


