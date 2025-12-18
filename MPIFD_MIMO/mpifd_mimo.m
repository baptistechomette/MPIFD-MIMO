function [uihandle] = mpifd_mimo(w, frf, vectAct, vectSens, vectOrder, useLRUR, useRealResidues, frfType, varargin)
% MPIFD_MIMO v1.0 beta
% Modal identification tool for dynamical structures.
% Based on LSCF and LSFD algorithms.
% The user inputs experimental (or simulated) data from single or multiple
% input, and multiple output frequency-domain vibration measurements 
% (Frequency Response Functions, 'frf', sampled at a list of angular 
% frequencies 'w'), then specifies a frequency range  and a list of orders 
% to try for modelling, then selects the relevant poles on the graphical 
% interface from the stabilization diagram, to then check how the FRFs are
% rendered by the modal approximation.
% 
% Usage: [h] = mpifd_mimo(w, frf, vectAct, vectSens, vectOrder, useLRUR, useRealResidues, frfType, varargin)
% w      vector of natural frequency in rad/s
% 
% frf    matrix of the complex frf, one column for one frf
% 
% vectAct       vector: ddl for actuator [ddl1, ddl2, ...]
%
% vectSens      vector: ddl for sensor [ddl1, ddl2, ...]
%
% vectOrder     vector: list of model orders to try for identification
% 
% useLRUR         = 1 with LR/UR residual terms
%                 = 0 without residual terms
%
% useRealResidues = 1 for identification using real residues
%                 (proportional damping)
%                 = 0 for identification using complex residues
% 
% FRFtype        'acc', 'vel' or 'disp' depending on the FRF measurements
% 
%
% OPTIONAL ARGUMENTS
% listIOVisu      2-column matrix containing Output/Input pair of dof indices
%                 of the FRF to be displayed. EX: [2 1; 3 2] for H21 and H32
% COORD           2D: N x M x 2 array with coordinates of the DOFs on the
%                   mesh : x = COORD(:,:,1) and y = COORD(:,:,2)
%                 3D: not implemented yet
%
% hfig            a figure handle to restart the GUI in the same window
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
% GNU CC BY-NC-SA
%
% Realease : v0 ... 2025
%
% Reference : 
% B. Chomette, J-L. Le Carrou, S. Karkar, F. Fabre, ..., JOSS, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('----------------------------------------------------------------------------------------------------')
disp(['MPIFD-MIMO: Modal Parameters Identification in the Frequency Domain -- Multi-Input - Multi-Output'])
disp('----------------------------------------------------------------------------------------------------')
fId=fopen('cartouche.txt','r');
credits = textscan(fId,'%s','Delimiter','');
credits = credits{1};
nLines = numel(credits);
for n=1:nLines
    disp(credits{n})
end
fclose(fId);

% disp('This research software is distributed under CC-XX-XX license.')
% disp('This is a research software, no guarantee included.')
fprintf('\n -- \n -- \n')
disp('Parsing input arguments...')
% checking input arguments

if size(w,2)~=1
    error('frequency vector w should be a row vector')
else
    state.w = w;
end

if size(frf,1)~=size(w,1)
    error('frf columns should have the same length as w')
end

if ~isvector(vectAct) 
    error('vectAct should be a *vector* containing the dof numbers of actuators')
end

if ~isvector(vectSens)
    error('vectSens should be a *vector* containing the dof numbers of sensors')
end
state.inputArgs.vectAct = vectAct;
state.inputArgs.vectSens = vectSens;

if size(frf,2)~=numel(vectAct)*numel(vectSens)
    error('frf should have as many columns as length(vectAct)*length(vectSens)')
else
    state.frf = frf;
end

if ~isvector(vectOrder)
    error('vectOrder should be a vector containing all the orders of successive identification')
else
    state.inputArgs.vectOrder = vectOrder;
end

if (useLRUR==1 || useLRUR==0)
    state.inputArgs.useLRUR = useLRUR;
else
    error('option useLRUR should be 1 or 0')
end

if isstring(useRealResidues)
    if strcmp(useRealResidues,'fnormal')
        state.inputArgs.useRealResidues=1;
    elseif strcmp(useRealResidues,'fcomplex')
        state.inputArgs.useRealResidues=0;
    else
        error('error parsing argument ''useRealResidues''. It should be a logical 1/0, or a string containing ''fnormal''/''fcomplex''')
    end
elseif islogical(useRealResidues)
    state.inputArgs.useRealResidues = useRealResidues;
else
    error('error parsing argument ''useRealResidues''. It should be a logical 1/0, or a string containing ''fnormal''/''fcomplex''')
end


if strcmp(frfType,'disp')
    state.frfTypeNum=1;
    disp(['+ Measured FRF provided are of Displacement type'])
elseif strcmp(frfType,'vel')
    state.frfTypeNum=2;
    disp(['+ Measured FRF provided are of Velocity type'])
elseif strcmp(frfType,'acc')
    state.frfTypeNum=3;
    disp(['+ Measured FRF provided are of Acceleration type'])
else
    error('Invalid type of FRF as input argument: must be a string containing ''disp'', ''vel'' or ''acc''.')
end
state.inputArgs.frfType = frfType;



%%% Managing optional user input: 
% - list of I/O pair of FRF to display
% - existing UIhandle
% - DOF coordinates for mode shape visu

% default if no option given
createFig=1;
createListIOVisu=1;
coordKnown = 0;

if ~isempty(varargin)
    for n=1:numel(varargin)
        arg=varargin{n};
        switch class(arg)
            case 'matlab.ui.Figure'
                uihandle = arg;
                createFig = 0 ;
            case 'double'
                if numel(size(arg))==2 && size(arg,2)==2
                    listIOVisu = arg;
                    createListIOVisu = 0;
                end
            case 'struct'
                if isfield(arg,'dimension')
                    meshStructure = arg ;
                    coordKnown=1;
                end
        end
    end
end



if createListIOVisu
    % when no user input is given concerning the FRFs to display, we
    % display all the possible ones
    nI=numel(vectAct);
    nO=numel(vectSens);
    for i=1:nI
        for o=1:nO
            listIOVisu(o+(i-1)*nO,:) = [vectSens(o) vectAct(i)] ;
        end
    end
    disp(['+ Displaying FRFs for all I/O couples. It might be a lot. '...
        'You can use the ''listIOVisu'' ' ...
        'extra argument to specify which FRF to display.'])
else
    if size(listIOVisu,2)~=2
        error('listIOVisu should contain 2 columns: [o i], to display H_oi')
    end
    disp(['+ Displaying FRFs based on provided ''listIOVisu'' argument.'])
end

if coordKnown
    state.meshDim = meshStructure.dimension ;
    fprintf('+ %dD mesh detected: ',state.meshDim)
    if state.meshDim<3
        state.meshCoord = meshStructure.COORD ;
        fprintf(['using provided coordinates for mode shape visualization.\n'])
    elseif state.meshDim==3
        % addpath('loadawobj2016') ;
        state.meshCoord = [] ;
        % hObj = createMeshObj(meshStructure.modelPath, meshStructure.modelOpt); % uses F.Fabre's external script to build the necessary 3D object
        state.model3D = meshStructure.model3D ;
        fprintf(['using provided 3D object desciption for mode shape visualization.\n'])
    end
else
    state.meshCoord = [] ;
    state.meshDim = 1 ; % default value: 1D object / standard plot function
    disp('+ No mesh provided. Falling back to default 1D mode shape visualization.')
end

state.MPCThresholdLow = .1 ; % threshold for "bad" MPC
state.MPCThresholdHigh = .6 ; % threshold for a "good" MPC



lft=25;
btm=50;
wdth=1000;
hght=650;
pos = [lft btm wdth hght];
if createFig
    disp('+ No figure handle provided, creating new UI')
    % Create new figure
    uihandle = uifigure('Name','MPIFD_MIMO: Modal Parameter Identification in the Frequency Domain','Position',[lft btm wdth hght]);
else
    % use existing figure handle
    % pos = uihandle.Position;
    % lft = pos(1) ;
    % btm = pos(2) ;
    % wdth = pos(3) ;
    % hght = pos(4) ;
    disp('+ Figure handle provided: using the same window.')
    set(uihandle,'Position',pos)
    for child=uihandle.Children
        delete(child)
    end
end


state.listIOVisu = listIOVisu ; % user input list of IO index pair for FRF visu (nFRFs x 2 matrix)
state.indexFRFVisu = [] ; % list of indices of plotted FRFs (used and defined in plotFRFCurves function)
state.nFRFs = [] ; % number of plotted FRF (used and defined in plotFRFCurves function)

disp('End of parsing input arguments.')
fprintf('\n')


%%% Start designing the window with components
fprintf('Designing the GUI...')
wPlots = .75*wdth ; % width of the plotting areas

tabGroupTop = uitabgroup(uihandle,'Position',[10 hght/2+5 wPlots hght/2-20]) ;

% plot 1 : Stabilization Chart
tabStabChart = uitab(tabGroupTop,"Title",'Stabilization Chart') ;
axStabChart = uiaxes(tabStabChart,'Position',[1 1 wPlots-10 hght/2-50],'HitTest', 'on');
axStabChart.Interactions = [regionZoomInteraction];
% ax.Toolbar.Visible = 'off';
colororder(axStabChart,{'k','b'})
axStabChart.set('box', 'on');


% Create a new tab for mode shapes visualization
tabModeShapes = uitab(tabGroupTop,"Title",'Mode Shape');
axModeShapes = uiaxes(tabModeShapes, ...
    'Position', [0 0 tabModeShapes.InnerPosition(3) tabModeShapes.InnerPosition(4)]);
if exist('meshStructure','var') 
    if meshStructure.dimension>1
        view(axModeShapes,3);
    end
end


% buttons and table - upper part

lButton = wPlots+15 ; % buttons and table left offset
wButton = wdth-lButton-10 ; % buttons and table width

% "Select poles" button
btnSelect = uibutton(uihandle,'state', 'Text','(De)Select Poles','Position',[lButton hght-50 wButton 35],...
    'ValueChangedFcn',@(s,e) switchSelectMode(s));
% "Clear selection" button
btnClear = uibutton(uihandle,'Text','Clear Selection','Position',[lButton hght-85 wButton 30],...
    'ButtonPushedFcn',@(s,e) clearSelections());
txtInfo = uilabel(uihandle,'Text','Dbl-click a row to display mode shape',...
    'Position', [lButton hght-115 wButton 25]);

% Table to show selected points
tblSelectedPts = uitable(uihandle,...
    'Position',[lButton hght/2+10 wButton hght-120-(hght/2+10)],...
    'ColumnWidth',{'3x', '3x', '2x'},...
    'DoubleClickedFcn', @(s,e) showModeShape(s),...
    'TooltipString', ['Double-click to show a mode shape,', char(10), ...
    '(only after you have computed the modes with the "Compute" button).'] );
tblSelectedPts.ColumnSortable = false ;
tblSelectedPts.ColumnName = {'F (Hz)', 'Damp. %', 'MPC'};

state.styleRedBackground = uistyle("BackgroundColor",[1 .6 .6]); % red background style
state.styleGreenBackground = uistyle("BackgroundColor",[.6 1 .6]); % green background style




% Lower part
% Plot 2 : graphics for measured/modelled FRFs
wPlots=.75*wdth; % keep the same for alignment formerly .71
% FRF mag plot
axFRFs = uiaxes(uihandle,'Position',[10 hght/8+10 wPlots hght/2*.75-10]);
title(axFRFs,'Measured and Modelled FRF');
ylabel(axFRFs,'|H| (dB)')
axFRFs.set('XTickLabel',[])
axFRFs.set('box', 'on')
axFRFs.Interactions = [regionZoomInteraction];

% FRF phase plot
axFRFsPhase = uiaxes(uihandle,'Position',[10 5 wPlots hght/8+10],'HitTest', 'on');
ylabel(axFRFsPhase,'<H (\pi rad)')
xlabel(axFRFsPhase,'F (Hz)')
axFRFsPhase.set('box','on')


% buttons checkboxes and tables
lButton = wPlots+15 ; % buttons and table left offset
wButton = wdth-lButton-10 ; % buttons and table width

% button to compute or update the simulated FRF with the selected modes
btnComputeModes = uibutton(uihandle,'Text','Compute Modes',...
    'Position',[lButton+.5*wButton hght/2-40 .5*wButton 40],...
    'ButtonPushedFcn',@(s,e) computeModes());

% add/remove LR term in computation
boxLRUR = uicheckbox(uihandle,'Text','Use LR/UR','Value',useLRUR,...
    'Position',[lButton hght/2-20 wButton*.45 20],...
    'ValueChangedFcn',@(s,e) switchLRUR(s));

boxUseRealResidues = uicheckbox(uihandle,'Text','Use Real Modes','Value',useRealResidues,...
    'Position',[lButton hght/2-45 wButton*.5 20],...
    'ValueChangedFcn',@(s,e) switchRealResidues(s));

btnDisplayRealModes = uibutton(uihandle,'Text','FRF with approx. Real Modes',...
    'Position',[lButton+15 hght/2-80 wButton-30 30],...
    'ButtonPushedFcn',@(s,e) displayRealModes(s));
if useRealResidues
    btnDisplayRealModes.Enable=0;
end

panelDisplayFRF = uipanel(uihandle,'Title','FRF to display',...
    'Position',[ lButton 10 wButton hght/2-100]);


%%% add a table with the list of O/I for FRF visualisation
tblFRFVisu = uitable(panelDisplayFRF,...
    'Position',[ 0 0 wButton hght/2-125],...
    'ColumnWidth',{'1x', '1x', 'fit'},...
    'TooltipString', ['Right click to add/delete a line.']);

% initialize table content
FRFInputList = listIOVisu(:,2);
FRFOutputList  =listIOVisu(:,1);
nlines = size(listIOVisu,1);
FRFVisibleTag = true([nlines,1]); % display all requested FRFs by default

tblFRFVisu.Data=table();
tblFRFVisu.Data.Output = FRFOutputList(:);
tblFRFVisu.Data.Input = FRFInputList(:);
tblFRFVisu.Data.Display = FRFVisibleTag(:);
tblFRFVisu.ColumnEditable = [true true true];
% tblFRFVisu.DisplayDataChangedFcn = @(s,e) plotFRFCurves();
tblFRFVisu.DisplayDataChangedFcn = @(s,e) updateFRFDisplayTable(s,e);

    function updateFRFDisplayTable(s,e)
        st=guidata(s.Parent);
        lengthTable = numel(st.tblFRFVisu.Data(:,1));
        nbAct = numel(st.inputArgs.vectAct);
        nbSens = numel(st.inputArgs.vectSens);
        if any(~sum(repmat(st.tblFRFVisu.Data.Output,1,nbSens)==repmat(st.inputArgs.vectSens(:)',lengthTable,1),2))
            errordlg('Wrong sensor (output) number, please correct.')
        elseif any(~sum(repmat(st.tblFRFVisu.Data.Input,1,nbAct)==repmat(st.inputArgs.vectAct(:)',lengthTable,1),2))
            errordlg('Wrong actuator (input) number, please correct.')
        else
            plotFRFCurves()
        end
    end


cm = uicontextmenu(uihandle); % add a context menu for right-click on the table
tblFRFVisu.ContextMenu = cm;

mDeleteRow = uimenu(cm); % add an option in the context menu to delete a row
mDeleteRow.Text = "Delete FRF from list";
mDeleteRow.MenuSelectedFcn = @deleteRow;

mAddRow = uimenu(cm) ; % add an option in the context menu to add a row
mAddRow.Text = "Add FRF to list";
mAddRow.MenuSelectedFcn = @(src,event) addRowToFRFList(src,event);

%cm.ContextMenuOpeningFcn = @(src,event) toggleVisibility(src,event,[mDeleteRow,mAddRow]);


%%% State properties of the UI
% These parameters and results can be accessed from the command line through
% 'st=guidata(fig)', as long as you provided the output argument 'fig' when
% launching the mpifd_mimo command, to get the handle of the UI.


state.fig = uihandle; % UI figure handle

state.tabGroupTop = tabGroupTop ; % tab group for top plot
state.tabStabChart = tabStabChart ; % tab containing Stab Chart axis
state.tabModeShapes = tabModeShapes ; % tab containing Mode Shape axis

state.axStabChart = axStabChart; % first axes handle (Stab Chart)
state.axFRFs = axFRFs; % second axes handle (FRFs, measured and modelled)
state.axFRFsPhase = axFRFsPhase; % second axes handle (FRFs, measured and modelled)
state.axModeShapes = axModeShapes ; % axes for showing mode shapes

state.tblSelectedPts = tblSelectedPts; % selected poles table handle
state.tblFRFVisu = tblFRFVisu ; % table with list of displayable FRFs

state.hFRFLines = [] ; % handles for measured FRF mag. plot (input from user)
state.hModFRFLines = []; % handle for simulated FRF mag. plot
state.hFRFPhaseLines = [] ; % handles for measured FRF phase plot
state.hModFRFPhaseLines = []; % handle for simulated FRF phase plot
state.hPlotModeShape = [] ; % handle for mode shape visualization

state.hboxLRUR = boxLRUR;
state.hboxUseRealResidues = boxUseRealResidues;
state.hbtnDisplayRealModes = btnDisplayRealModes;

state.selectedPoles = table('Size',[0 3],...
    'VariableTypes',{'double','double', 'matlab.graphics.Graphics'},...
    'VariableNames',{'Freq', 'Order', 'Marker'}); % table with fields x,y,marker

state.hPolesMarkers = []; % plotted markers handles (all identified poles)
state.n = []; % number of model orders processed
state.POLES = []; % list of poles for each model order

state.id.fst = [] ; % frequencies of currently selected poles
state.id.xist = [] ; % damping of currently selected poles
state.id.phir = []; % computed real modes based on selected poles
state.id.phi = []; % computed complex modes based on selected poles
state.id.rntot=[]; % computed residues based on selected poles
state.id.LRtot=[]; % computed low frequency additional correction term
state.id.URtot=[]; % computed high frequency additional correction term
state.id.frfmtot=[] ; % computed modal FRF without re-computing real modes

state.dFactor = 3*1e-2 ; % multiplier in front of the mode shape to show deformed geometry (3D object visu)

fprintf('........... done.\n')
guidata(uihandle,state);


fprintf('\n -- \n')
fprintf('Computing stabilization chart:\n')
plotStabChart();
fprintf('........................... done.\n')

fprintf('Plotting measured FRF provided...')
plotFRFCurves();
fprintf('........................................... done.\n')
fprintf('\n -- \n -- \n')


% link x axes of all 3 plots
linkaxes([state.axStabChart,state.axFRFs,state.axFRFsPhase],'x')

disp('All set. You can know use the GUI.')

%%% Nested functions -----------------------------------------------------

    function plotStabChart()
        st=guidata(uihandle);
        st.axStabChart.NextPlot = 'add';
        % compute and plot the FRF kind-of-sum
        Hsum = frf_sum(st.frf) ;
        Hsum_dB = 20*log10(abs(Hsum)) ;
        yyaxis(st.axStabChart,'left')
        plot(st.axStabChart,w/2/pi, Hsum_dB,'DisplayName','FRF moduli sum','HitTest','on','PickableParts','none') ;
        xlabel(st.axStabChart,'Frequency (Hz)')
        ylabel(st.axStabChart,'FRF sum')

        % identify poles for each model order in the list of orders
        POLES = lscf(st.w, st.frf, st.inputArgs.vectOrder) ;
        n = length(POLES) ; % number of model order cases identified
        
        % plot the (stable) poles that have been computed
        yyaxis(st.axStabChart,'right');
        ylabel(st.axStabChart,'model order')
        set(st.axStabChart,'YLim',[min(vectOrder) max(vectOrder)*1.1]);
        hPolesMarkers = gobjects(n,1);
        for kp=1:n
            hPolesMarkers(kp)=plot(st.axStabChart,[POLES(kp).freq]', [repmat(POLES(kp).order,size(POLES(kp).freq))], 'bo','MarkerSize', 3,'HitTest','on','PickableParts','all') ;
        end
        ylabel(st.axStabChart,'model order')


        st.hPolesMarkers = hPolesMarkers; % plotted line handles
        st.n = n; % number of model orders processed
        st.POLES = POLES; % list of poles for each model order
        
        yyaxis(st.axStabChart,'left')

        guidata(uihandle,st);
    end

    
    function plotFRFCurves()
        %%% plot FRFs on the lower graph
        % get the data from the gui
        st=guidata(uihandle);
        axesPos = st.axStabChart.InnerPosition;

        % clear FRF plot
        myax = st.axFRFs;
        % clear the graph
        for obj=myax.Children
            delete(obj);
        end

        % reset some properties
        myax.InnerPosition(1)=axesPos(1)+st.tabGroupTop.InnerPosition(1)-1;
        myax.InnerPosition(3)=axesPos(3);
        set(myax,'ColorOrder','default')
        set(myax, 'ColorOrderIndex', 1)


        % local variables

        visuData = st.tblFRFVisu.Data;
        FRFOut = visuData.Output;
        FRFIn = visuData.Input;
        FRFVisi = visuData.Display;
        listAct = st.inputArgs.vectAct;
        listSens = st.inputArgs.vectSens;

        listIOVisuMasked = [FRFOut(FRFVisi) FRFIn(FRFVisi)];
        nFRFs = size(listIOVisuMasked,1);

        % case when there is nothing to display
        if nFRFs==0
            st.axFRFs = myax;
            st.listIOVisuMasked = [];
            st.indexFRFVisu = [];
            st.nFRFs = nFRFs;
            st.leg = {};
            myax=st.axFRFsPhase;
            for obj=myax.Children
                delete(obj);
            end
            st.axFRFsPhase=myax;
            guidata(uihandle,st);
            plotModFRFs()
            return
        end
        % Otherwise, we proceed with plotting FRFs Mag and Phase

        % compute the list of indices to use for the 2nd dimension of 'frf'
        nOutputs=numel(listSens);
        indexFRFVisu = zeros(nFRFs,1);
        for l=1:nFRFs
            indexFRFVisu(l) = (find(listAct==listIOVisuMasked(l,2))-1)*nOutputs+find(listSens==listIOVisuMasked(l,1));
        end
        

        myfrf = st.frf(:,indexFRFVisu); % only load what we need to display
        frf1log = 20*log10(abs(myfrf)) ;

        hFRFLines = plot(myax,st.w/2/pi, frf1log,'-','LineWidth', 1) ;
        % compile and add the legend
        leg=cell(nFRFs,1);
        for l=1:nFRFs
            leg{l}=['Hmes[' num2str(listIOVisuMasked(l,1)) ',' num2str(listIOVisuMasked(l,2)) ']'];
        end
        legend(myax, leg, 'Orientation' , 'vertical');
        st.axFRFs = myax;

        % plot the phase in the lower plot
        myax=st.axFRFsPhase;
        for obj=myax.Children
            delete(obj);
        end
        
        myax.InnerPosition(1)=axesPos(1)+st.tabGroupTop.InnerPosition(1)-1;
        myax.InnerPosition(3)=axesPos(3);

        set(myax,'ColorOrder','default')
        set(myax, 'ColorOrderIndex', 1)
        switch st.frfTypeNum
            case 1
                FRFPhase = mod(angle(myfrf),pi);
                % ylim(myax,[0,pi])
                % yticks(myax,[0 pi/4 pi/2 .75*pi pi])
                % yticklabels(myax,["0", "\pi/4", "\pi/2", "3\pi/4", "pi"])
            case 2
                FRFPhase = mod(angle(myfrf)+pi/2,pi)-pi/2 ;
                % ylim(myax,[-pi/2,pi/2])
                % yticks(myax,[-pi/2 -pi/4 0 pi/4 pi/2])
                % yticklabels(myax,["-\pi/2", "-\pi/4", "0", "\pi/4", "\pi/2"])
            case 3
                FRFPhase = mod(angle(myfrf)+pi,pi)-pi ;
                % ylim(myax,[-pi,0])
                % yticks(myax,[-pi -.75*pi -pi/2 -pi/4 0])
                % yticklabels(myax,["-\pi", "-3\pi/4","-\pi/2", "-\pi/4", "0"])
        end
        hFRFPhaseLines = plot(myax,st.w/2/pi, FRFPhase/pi,...
            '-','LineWidth', 1, 'HandleVisibility','on') ;
        
        st.axFRFsPhase = myax ;
        % store new data
        st.listIOVisuMasked = listIOVisuMasked;
        st.indexFRFVisu = indexFRFVisu;
        st.hFRFLines = hFRFLines;
        st.hFRFPhaseLines = hFRFPhaseLines;
        st.nFRFs = nFRFs;
        st.leg = leg;
       
        guidata(uihandle,st);

        % if there was some modal FRF displayed, replot that too
        if ~isempty(st.hModFRFLines)
            plotModFRFs()
        end
        
    end

    function plotModFRFs()
        st=guidata(uihandle);
        % plotting the FRFs based on the identified modal model
        myax = st.axFRFs;
        set(myax,'ColorOrder','default')
        set(myax, 'ColorOrderIndex', 1)
        hold(myax,'on')
        
        % clear (possibly) previously computed FRFs
        if ~isempty(st.hModFRFLines)
            for line=st.hModFRFLines
                delete(line);
            end
            for line=st.hModFRFPhaseLines
                delete(line);
            end
        end

        
        
        % plot currently computed modal FRFs
        
        % some local variables
        nFRFs = st.nFRFs ; % number of FRFs to display
        if nFRFs==0
            return;
        end
        
        
        indexFRFVisu = st.indexFRFVisu; % index of the corresponding column
        FRFOut = st.listIOVisuMasked(:,1); % list of corresponding Sensor dofs
        FRFIn = st.listIOVisuMasked(:,2); % list of corresponding Actuator dofs

        % load the identified model's numerically simulated FRFs
        Hnum=st.id.Hnum(:,indexFRFVisu); % load only needed columns for display
        Hnumlog = 20*log10(abs(Hnum)) ; % compute mag in dB

        % plot the magnitude
        hModFRFLines=plot(myax,st.w/2/pi, Hnumlog ,...
            ':', 'LineWidth', 1,'DisplayName','FRFmod');
        % compile the legend
        leg2=cell(nFRFs,1);
        for ll=1:nFRFs
            leg2{ll}=['Hmod[' num2str(FRFOut(ll)) ',' num2str(FRFIn(ll)) ']'];
        end
        legend(myax,[st.leg, leg2], 'Orientation' , 'vertical', 'NumColumns', 2);
        hold(myax,'off')
        st.axFRFs = myax;

        % plot the phase
        myax=st.axFRFsPhase;
        set(myax, 'ColorOrderIndex', 1)
        set(myax,'ColorOrder','default')
        hold(myax,'on')
        switch st.frfTypeNum
            case 1
                modFRFPhase = mod(angle(Hnum),pi) ;
            case 2
                modFRFPhase = mod(angle(Hnum)+pi/2,pi)-pi/2 ;
            case 3
                modFRFPhase = mod(angle(Hnum)+pi,pi)-pi ;
        end                
        hModFRFPhaseLines = plot(myax,st.w/2/pi, modFRFPhase/pi,...
            ':','LineWidth', 1, 'HandleVisibility','on') ;
        hold(myax,'off')

        st.hModFRFLines = hModFRFLines;
        st.hModFRFPhaseLines = hModFRFPhaseLines ;
        st.axFRFsPhase = myax;

        guidata(uihandle,st);
    end



    function switchSelectMode(btn)
        st = guidata(uihandle);
        if btn.Value==1
            % after clicking on "Select Poles"
            st.tabGroupTop.set('SelectedTab',st.tabStabChart);
            st.fig.Pointer = 'crosshair';
            disp('Entering pole selection mode. Click the same button to validate your selection.')
            % Capture clicking on lines
            for L = st.hPolesMarkers'
                L.ButtonDownFcn = @lineClick;
            end
            if ~isempty(st.selectedPoles.Marker)
                for P = st.selectedPoles.Marker'
                    P.ButtonDownFcn = @(s,e) removePoint(s);
                end
            end
            btn.Text = 'Validate Selection';
            fprintf('Points selected: ')
            textPlus='';
            for iPts=1:numel(st.selectedPoles.Freq)
                textPlus = [textPlus, '+'];
            end
            fprintf('%s', textPlus)
        else
            % after clicking on 'Validate Selection'
            st = guidata(uihandle);
            st.fig.Pointer = 'arrow';
            st.axStabChart.ButtonDownFcn = '';
            for L = st.hPolesMarkers'
                L.ButtonDownFcn = '';
            end
            if ~isempty(st.selectedPoles.Marker)
                for P = st.selectedPoles.Marker'
                    P.ButtonDownFcn = '';
                end
            end
            btn.Text = '(De)Select Poles';
            fprintf(' done.\n')
            disp(['',num2str(numel(st.selectedPoles.Freq)), ' poles selected.'])
        end
        guidata(uihandle,st);
    end


    function lineClick(src,event)
        % Clicked on a plotted line: snap to nearest data point
        st = guidata(uihandle);
        % identify which line
        lineHandle = src;
        xd = lineHandle.XData;
        yd = lineHandle.YData;
        cp = st.axStabChart.CurrentPoint;
        xclick = cp(1,1);
        % find nearest x index
        [~, idx] = min(abs(xd - xclick));
        xpt = xd(idx);
        ypt = yd(idx);
        addSelectedPoint(xpt,ypt);
    end

    function addSelectedPoint(xp,yp)
        % after detecting a point that was clicked on, add it to the list
        % of selected points

        st = guidata(uihandle);
        % Create a marker on the axes
        hMarker = plot(st.axStabChart,xp,yp,'rs','MarkerSize',8, ...
            'ButtonDownFcn',@(s,e) removePoint(s)); % click marker to remove
        % nPoint = size(st.selectedPoles,1)+1;
        t=table();
        t.Freq = [st.selectedPoles.Freq ; xp];
        t.Order = [st.selectedPoles.Order ; yp];
        t.Marker = [st.selectedPoles.Marker ; hMarker];
        st.selectedPoles = sortrows(t,'Freq'); % sort poles by frequency
        fprintf('+')
        guidata(uihandle,st);
        updateSelectedPolesTable();
    end

    function removePoint(markerHandle)
        % if clicked on an already selected point, remove that point from
        % the list of selected points

        st = guidata(uihandle);
        % find and remove the point associated with markerHandle
        idx = find(st.selectedPoles.Marker(:) == markerHandle);
        if ~isempty(idx)
            delete(st.selectedPoles.Marker(idx));
            st.selectedPoles(idx,:) = [];
            fprintf('\b')
            guidata(uihandle,st);
            updateSelectedPolesTable();
        end
    end

    function clearSelections()
        % remove all selected points
        st = guidata(uihandle);
        if ~isempty(st.selectedPoles)
            for mark=st.selectedPoles.Marker(:)
                if isgraphics(mark)
                    delete(mark);
                end
            end
        end
        st.selectedPoles = table('Size',[0 3],...
    'VariableTypes',{'double','double', 'matlab.graphics.Graphics'},...
    'VariableNames',{'Freq', 'Order', 'Marker'}); % empty the table
        guidata(uihandle,st);
        updateSelectedPolesTable();
        disp("Pole selection cleared.");
    end


    function updateSelectedPolesTable()
        % display the content of the table showing selected poles

        st=guidata(uihandle);
        POLES=st.POLES;

        % Prepare data: Nx2 array of [freq damping] of selected poles
        numPointsSelected = size(st.selectedPoles,1) ;
        poleFreq = zeros(numPointsSelected,1) ;
        poleDamping = zeros(numPointsSelected,1) ;
        % poleOrder = zeros(numPointsSelected,1,'int8') ;
        
        
        for iSel=1:numPointsSelected
            [~, iorder] = min(abs(st.selectedPoles.Order(iSel) - [POLES.order])) ;
            freqThisOrder = [POLES(iorder).freq] ;
            dampThisOrder = [POLES(iorder).xi] ;
            [~, idx] = min(abs(st.selectedPoles.Freq(iSel) - freqThisOrder)) ;
            poleFreq(iSel) = freqThisOrder(idx) ;
            poleDamping(iSel) = dampThisOrder(idx) ;
            % poleOrder(i) = int8(order);
        end

        data=table();
        data.Freq = poleFreq;
        data.Damp = poleDamping*100;
        %data.Ord = poleOrder;
        data.MPC = nan(size(poleFreq)); % MPC=NaN before computing modes
        % data = sortrows(data,1); % sort poles by frequency
        % store data into GUI
        st.tblSelectedPts.Data = data;
        addStyle(st.tblSelectedPts,uistyle("BackgroundColor",[1 1 1]))
        st.id.fst = data.Freq;
        st.id.xist = data.Damp/100;
        guidata(uihandle,st);
    end

    function deleteRow(src,event)
        % function for removing a row from the list of FRFs to display
        fig=src.Parent;
        st=guidata(fig);
        % hpltMS=st.hPlotModeShape;
        tbl = event.ContextObject;
        row = event.InteractionInformation.Row;
        tbl.Data(row,:) = [];
        disp('FRF removed from display.')
        plotFRFCurves()
        % update markers on mode shape if needed
        if ~isempty(st.hPlotModeShape)
            st.tblSelectedPts.Selection=[st.modeNum 1];
            showModeShape(st.tblSelectedPts)
        end
    end


    function addRowToFRFList(src,event)
        % Add a new row to the FRF visualization table
        fig=src.Parent;
        st=guidata(fig);
        hpltMS=st.hPlotModeShape;
        tbl = event.ContextObject;
        newRow = {vectSens(1) , vectAct(1) , true}; % Initialize with default values and visible
        tbl.Data = [tbl.Data; newRow]; % Append new row
        disp('Added a new FRF to display. You can edit the I/O indices in the table.')
        plotFRFCurves()
        % update markers on mode shape if needed
        if ~isempty(hpltMS)
            st.tblSelectedPts.Selection=[st.modeNum 1];
            showModeShape(st.tblSelectedPts)
        end
        guidata(fig,st);
    end


    function toggleVisibility(src,event,menus)
        % allow visibility of the menu only on the table
        row = event.InteractionInformation.Row;
        rowClicked = ~isempty(row);
        for m=menus
            m.Visible = rowClicked;
        end
    end


    function switchLRUR(box)
        % only used to display a message in the command window when
        % clicking the "Use LR/UR" check box
        if box.Value
            disp('''LR/UR'' activated - Hit "Compute Modes" again.')
        else
            disp('''LR/UR'' removed - Hit "Compute Modes" again.')
        end
    end

    function switchRealResidues(box)
        % disable the button to recompute FRF with real modes if we set the
        % computation to use real modes already
        st=guidata(uihandle);
        if box.Value
            disp('Switching to Real residues - Hit "Compute Modes" again.')
            st.hbtnDisplayRealModes.Enable=0;
        else
            disp('Switching to Complex residues - Hit "Compute Modes" again.')
            st.hbtnDisplayRealModes.Enable=1;            
        end
        guidata(uihandle,st);
    end

    function computeModes()
        % Main mode compute routine using LSFD

        st = guidata(uihandle);
        if isempty(st.selectedPoles.Freq)
            uialert(st.fig,'No points selected','Warning');
            return;
        else
            nPointsSelected = numel(st.selectedPoles.Freq);
        end
        % Prepare data: Nx2 array of [freq, damping] of selected poles
        fprintf('Computing residues for the selected %d poles',nPointsSelected)
        if st.hboxUseRealResidues.Value
            fprintf('\n+ Using Real-valued residues')
        else
            fprintf('\n+ Using Complex-valued residues')
        end
        if  st.hboxLRUR.Value==1
            switch st.frfTypeNum
                case 1
                    fprintf('\n+ LR/UR residual terms: Displacement type\n...')
                case 2
                    fprintf('\n+ LR/UR residual terms: Velocity type\n...')
                case 3
                    fprintf('\n+ LR/UR residual terms: Acceleration type\n...')
            end
        end
            

        % compute residues
        [rntot,LRtot, URtot, frfmtot] = lsfd(st.w, st.frf,...
            st.id.fst, st.id.xist, st.hboxLRUR.Value, ...
            st.hboxUseRealResidues.Value, st.inputArgs.frfType) ;
        Hnum = frfmtot;
        fprintf(' done.\n')

        % from the residues compute de the mode shapes 
        fprintf('Computing corresponding mode shapes...')
        [st.id.phir, st.id.phi, st.id.residues] = residues2modes(rntot,...
            st.inputArgs.vectAct, st.inputArgs.vectSens, st.id.fst, st.id.xist,...
            st.hboxUseRealResidues.Value) ;
        fprintf(' done.\n')
        % store results in state variable
        st.id.rntot=rntot; st.id.LRtot=LRtot; st.id.URtot=URtot;
        st.id.frfmtot=frfmtot;
        st.id.Hnum = Hnum;
        
        % compute MPC of the modes and add them to the table
        if ~st.hboxUseRealResidues.Value
            fprintf('Computing Modal Phase Colinearity (MPC) for each mode...')
            MPC = mpc_mpifd(st.id.phi);
            fprintf(' done.\n')
        else
            MPC = nan(size(st.tblSelectedPts.Data.MPC)) ;
        end
        st.tblSelectedPts.Data.MPC = MPC ;

        % set the background in red for modes with low MPC
        removeStyle(st.tblSelectedPts);
        if ~st.hboxUseRealResidues.Value
            rowBadMPC = find(MPC<st.MPCThresholdLow);
            rowGoodMPC = find(MPC>st.MPCThresholdHigh);
            if ~isempty(rowBadMPC)
                addStyle(st.tblSelectedPts,st.styleRedBackground,"row",rowBadMPC);
            end
            if ~isempty(rowGoodMPC)
                addStyle(st.tblSelectedPts,st.styleGreenBackground,"row",rowGoodMPC);
            end
        end

        % save state variable and update plot
        guidata(uihandle,st)
        fprintf('Updating FRF plots...')
        plotModFRFs()
        fprintf(' done.\n')
    end

    function displayRealModes(s)
        % Recompute modal FRFs using real modes (when using LSFD with 
        % complex residues)
        st=guidata(s.Parent);
        fprintf(['Re-computing real modes from complex residues...']);
        st.id.Hnum = frf_real_modes(st.id.phir, st.id.fst, st.id.xist,...
            st.w, st.id.LRtot, st.id.URtot, st.inputArgs.vectAct,...
            st.inputArgs.frfType) ;
        fprintf(' done.\n')
        guidata(uihandle,st)
        fprintf('Updating FRF plots...')
        plotFRFCurves()
        fprintf(' done.\n')
    end



    function showModeShape(s)
        % when dbl click on a cell of the table of selected poles
        % show the associated mode shape
        modeNum = s.Selection(:,1) ; % row of clicked cell
        st=guidata(s.Parent);
        
        st.modeNum = modeNum;
        guidata(s.Parent, st);

        % switch function if we are dealing with a 2D or 3D object
        switch st.meshDim
            case 1
                showModeShape1D(s);
            case 2
                showModeShape2D(s);
            case 3
                showModeShape3D(s);
        end
    end % showModeShape


    function showModeShape1D(s)
        modeNum = s.Selection(:,1) ; % row of clicked cell
        st=guidata(s.Parent);
                disp(['Showing mode(s) number: ',int2str(modeNum(:)'), '...'])

        if isempty(st.id.rntot)
            errordlg('No mode computed yet')
            return
        elseif isempty(st.id.phi)
            [st.id.phir, st.id.phi, st.id.residues] = residues2modes(st.id.rntot,...
            st.inputArgs.vectAct, st.inputArgs.vectSens, st.id.fst, st.id.xist,...
            st.hboxUseRealResidues.Value) ;
        end
        st.tabGroupTop.SelectedTab = st.tabModeShapes;
        cla(st.axModeShapes);
        hpltMS = plot(st.axModeShapes,st.inputArgs.vectSens,st.id.phir(:,modeNum),'--d') ;
        nM=numel(modeNum);
        legMS = cell(nM,1);
        for im=1:nM
            legMS{im} = ['Mode #', int2str(modeNum(im))];
        end
        legend(hpltMS, legMS)
        xlabel(st.axModeShapes,'DOF');
        ylabel(st.axModeShapes,'Magnitude');
        st.axModeShapes.XTick = st.inputArgs.vectSens ;

        fprintf('\b done.\n')
        st.hPlotModeShape = hpltMS ;
        guidata(uihandle,st) ;
    end % showModeShape1D

    function showModeShape2D(s)
        % when dbl click on a cell of the table of selected poles
        % show the associated mode shape in 3D, for a 2D structure
        modeNum = s.Selection(1) ; % row of clicked cell (only one possible here)
        fprintf('Showing mode number: %d...',modeNum)
        st=guidata(s.Parent);

        if isempty(st.id.rntot)
            fprintf('error: no mode computed yet. Hit Compute.\n')
            errordlg('No mode computed yet. Hit Compute.\n')
            return
        elseif isempty(st.id.phi)
            [st.id.phir, st.id.phi, st.id.residues] = residues2modes(st.id.rntot,...
            st.inputArgs.vectAct, st.inputArgs.vectSens, st.id.fst, st.id.xist,...
            st.hboxUseRealResidues.Value) ;
        end

        % using real-valued modes (product of residue2modes)
        PHI = st.id.phir(:,modeNum);
        PHI = PHI./max(abs(PHI),[],1);

        meshCoord = st.meshCoord ;
        x = meshCoord(:,1) ; 
        xUnique = unique(x);
        y = meshCoord(:,2) ;
        yUnique = unique(y) ;
        [X, Y] = meshgrid(xUnique, yUnique);
        PHIGrid = griddata(x, y, PHI, X, Y);

        st.tabGroupTop.SelectedTab = st.tabModeShapes; % set the focus to the Mode Shape tab
        
        % save viewing angles
        if isgraphics(st.axModeShapes) && isvalid(st.axModeShapes)
            viewModeShape = st.axModeShapes.View;
            cla(st.axModeShapes)
            delete(st.axModeShapes)
        else
            viewModeShape = 3;
        end

        % clear the previous graph
        if ~isempty(st.tabModeShapes.Children)
            for chld=st.tabModeShapes.Children
                cla(chld)
                delete(chld)
            end
        end

        
        axMS = uiaxes(st.tabModeShapes, ...
            'Position', [0 0 st.tabModeShapes.InnerPosition(3) st.tabModeShapes.InnerPosition(4)],...
            'XTickLabel', [],...
            'YTickLabel', [],...
            'ZTickLabel', [],...
            'Color', 'none',...
            'TickLength', [0 0],...
            'Box', 'on', 'HitTest', 'on');
        % 'Projection', 'perspective') ;
        % 'DataAspectRatio', [1 1 1],...
        hold(axMS,'on') ;
        colormap(axMS, 'jet');
        
        % axMS = st.axModeShapes ;

        %%% build the visu
        % plot grid points on the z=0 plane
        % hpltMS{1} = mesh(axMS, X', Y', PHIGrid', 'EdgeColor','k', 'LineStyle','--', 'Marker','.', 'FaceColor','none') ; 
        % add DOF number to datatips
        row = dataTipTextRow("DOF#",(1:numel(x)));
        % hpltMS{1}.DataTipTemplate.DataTipRows(end+1) = row;
        
        % add markers for actuators and sensors selected for FRF visu
        actuatorPointNumber = unique(st.tblFRFVisu.Data.Input);
        hpltMS{2}=plot3(axMS, x(actuatorPointNumber), y(actuatorPointNumber), PHI(actuatorPointNumber), 'rs', 'MarkerSize',12, 'MarkerFaceColor','r') ;
        sensorPointNumber = unique(st.tblFRFVisu.Data.Output);
        hpltMS{3}=plot3(axMS, x(sensorPointNumber), y(sensorPointNumber), PHI(sensorPointNumber), 'bo', 'MarkerFaceColor','b') ;

        % add mesh with faces (surf)
        % hpltMS{4} = mesh(axModeShapes, X, Y, PHIGrid, 'FaceAlpha','0') ;
        hpltMS{4} = mesh(axMS, X', Y', PHIGrid', 'EdgeColor', 'k',...
            'EdgeAlpha', '0.2', 'FaceColor', 'interp', 'FaceAlpha','1',...
            'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none') ;
        hpltMS{4}.DataTipTemplate.DataTipRows(end+1) = row; 

        % make everything pretty
        view(axMS, viewModeShape); 
        daspect(axMS,[1 1 1/st.dFactor])
        title(axMS,['Mode #', int2str(modeNum), ': f = ', num2str(st.id.fst(modeNum), '%.1f'), 'Hz']) ;
        legend(axMS, 'off') ;
        % colormap(st.axModeShapes, 'jet') ;
        % axis(st.axModeShapes, 'off');
        % box(st.axModeShapes, 'on')
        % st.axModeShapes.TickLength=[0 0];
        % st.axModeShapes.set('XTickLabel', [],...
        %                                 'YTickLabel', [],...
        %                                 'ZTickLabel', []);
        % st.axModeShapes.Projection='perspective';
        % st.axModeShapes.BoxStyle = "full";

        % add small axes to recognise x,y,z
        lineLen = max(max(xUnique)-min(xUnique),max(yUnique)-min(yUnique))/20;
        dAsp = daspect(axMS);
        plot3(axMS, [0 lineLen*dAsp(1)], [0 0], [0 0], 'b', 'LineWidth',2)        
        plot3(axMS, [0 0], [0 lineLen*dAsp(2)], [0 0], 'g', 'LineWidth',2)        
        plot3(axMS, [0 0],[0 0], [0 lineLen*dAsp(3)], 'r', 'LineWidth',2)        

        fprintf(' done.\n')
        
        % store handle to the plot elements into the GUI
        st.hPlotModeShape = hpltMS ;
        st.axModeShapes = axMS ;
        guidata(s.Parent,st) ;
    end % showModeShape2D

    function showModeShape3D(s)
        % when dbl click on a cell of the table of selected poles
        % show the associated mode shape for the 3 2D-projections of the
        % struture along axes Z, Y, and X.
        modeNum = s.Selection(1) ; % row of clicked cell
        fprintf('Showing mode number: %d...',modeNum)
        st=guidata(s.Parent);
        if isempty(st.id.rntot)
            fprintf('error: no mode computed yet. Hit Compute.\n')
            return
        elseif isempty(st.id.phi)
            [st.id.phir, st.id.phi, st.id.residues] = residues2modes(st.id.rntot,...
            st.inputArgs.vectAct, st.inputArgs.vectSens, st.id.fst, st.id.xist,...
            st.hboxUseRealResidues.Value) ;
        end

        st.tabGroupTop.SelectedTab = st.tabModeShapes; % set the focus to the Mode Shape tab
        
        % clear the previous graph, but keep view angles
        if ishandle(st.axModeShapes(1)) && isvalid(st.axModeShapes(1)) 
            viewModeShape = st.axModeShapes(1).View;
        else
            viewModeShape = 3;
        end
        for ax=st.tabModeShapes.Children
            cla(ax)
            delete(ax)
        end

        % mode shape to plot
        PHI = st.id.phir(:,modeNum);
        
        % 3D model of the structure (user provided, see input arguments)
        model3D = st.model3D ;

        nbMeasureDirections = max(model3D.frfModelLink(:,2)) ;
        
        kk=1; 
        for ll=nbMeasureDirections
        st.axModeShapes(ll, kk) = uiaxes(st.tabModeShapes, 'Position', [(ll-1)/max(model3D.frfModelLink(:,2))*st.tabModeShapes.InnerPosition(3) 0 1/max(model3D.frfModelLink(:,2))*st.tabModeShapes.InnerPosition(3) 1*st.tabModeShapes.InnerPosition(4)],...
                                       'DataAspectRatio', [1 1 1],...
                                       'XTickLabel', [],...
                                       'YTickLabel', [],...
                                       'ZTickLabel', [],...
                                       'Color', 'none',...
                                       'TickLength', [0 0],...
                                       'Box', 'on', ...
                                       'Projection', 'perspective', ...
                                       'UserData',[ll kk]) ;
        hold(st.axModeShapes(ll,kk),'on') ;
        colormap(st.axModeShapes(ll,kk), 'jet');
        end

        
        %%% build the 3D visu
        % using pieces of code from Francois Fabre 
        % (c) 2025 - LAM - Sorbonne University
        % Adapted by Sami Karkar - LAM - Sorbonne University - 11-2025

        dFactor = st.dFactor;
        phi_resh = reshape(dFactor*PHI./max(abs(PHI),[],1), [1, size(PHI)]) ;
        phi_rgp = NaN(3, size(PHI,2), max(model3D.frfModelLink(:,1)), nbMeasureDirections) ;
        linIdx = sub2ind(size(phi_rgp, [3, 4]), model3D.frfModelLink(:,1), model3D.frfModelLink(:,2)) ;

        phi_rgp(:,:, linIdx) = permute(phi_resh.*model3D.vertexNormals, [1 3 2]) ;
        phi_rgp = permute(phi_rgp, [1 3 2 4]) ;

        vertCoord = permute(model3D.vertices + phi_rgp, [2 1 3 4]) ; % coordinates of deformed geometry

        vertColor = zeros(size(PHI,2), max(model3D.frfModelLink(:,1)), max(model3D.frfModelLink(:,2))) ;
        vertColor(:, linIdx) = reshape(phi_resh, size(PHI)).' ;
        vertColor = permute(vertColor, [2 1 3]) ;

        
        % NOTE: for showing multiple modes at a time, one could add tabs, or open a new window
        % because it is going to be very crowded

        % used later to draw x,y,z small color axes
        x=vertCoord(:,1);
        y=vertCoord(:,2);
        z=vertCoord(:,3);
        lineLen = max([max(x)-min(x), max(y)-min(y), max(z)-min(z)])/20;

        % drawing starts here
        kk=1; % for now we're showing 1 mode at a time - Sami - 25.10.2025
        for ll = 1:nbMeasureDirections
                            
            if ~isempty(model3D.faces)
                hpltMS{4,ll,kk}=[patch(st.axModeShapes(ll, kk), 'Vertices', vertCoord(:,:, kk, ll),...
                                             'Faces', model3D.faces.',...
                                             'FaceVertexCData', vertColor(:, kk, ll),...
                                             'FaceColor', 'interp', ...
                                             'EdgeColor', 'k',...
                                             'EdgeAlpha', 0.2,...
                                             'Marker', 'none',...
                                             'LineWidth', 0.5)];
                hpltMS{4,ll,kk}.HitTest = 'off' ;

            end
            
            if isfield(model3D, 'lines')
                hpltMS{4,ll,kk}(2)=patch(st.axModeShapes(ll, kk), 'Vertices', vertCoord(:,:, kk, ll),...
                                             'Faces', model3D.lines.',...
                                             'FaceVertexCData', vertColor(:, kk, ll),...
                                             'FaceColor', 'interp', ...
                                             'EdgeColor', 'interp',...
                                             'EdgeAlpha', 1,...
                                             'LineWidth', 2,...
                                             'Marker','none') ;
                hpltMS{4,ll,kk}(2).HitTest = 'off' ;
            end
            
            [limit_min, limit_max] = bounds(vertColor(:, kk, ll)) ;
            clim(st.axModeShapes(ll,kk),0.9*[-1 1]*mean(abs([limit_min, limit_max]))) ;

            title(st.axModeShapes(ll, kk), ['Mode #', int2str(modeNum), ' : f = ',...
            num2str(st.id.fst(modeNum), '%.1f'), 'Hz , damp. = ', ...
            num2str(st.id.xist(modeNum)*100, '%.1f'), '%'])


            % add markers for actuators and sensors selected for FRF visu
            actuatorPointNumber = unique(st.tblFRFVisu.Data.Input);
            hpltMS{2,ll,kk}=plot3(st.axModeShapes(ll,kk), ...
                x(actuatorPointNumber), y(actuatorPointNumber), ...
                z(actuatorPointNumber), 'rs', ...
                'MarkerSize',12, 'MarkerFaceColor','r', ...
                'HitTest', 'off') ;
            sensorPointNumber = unique(st.tblFRFVisu.Data.Output);
            hpltMS{3,ll,kk}=plot3(st.axModeShapes(ll,kk), ...
                x(sensorPointNumber), y(sensorPointNumber), ...
                z(sensorPointNumber), 'bo', 'MarkerFaceColor','b', ...
                'HitTest','off') ;

            % make everything pretty
            view(st.axModeShapes(ll,kk), viewModeShape);
            legend(st.axModeShapes(ll,kk), 'off') ;

            % add small axes to recognise x,y,z
            daspect(st.axModeShapes(ll,kk),[1 1 1]);
            dAsp = daspect(st.axModeShapes(ll,kk));
            plot3(st.axModeShapes(ll,kk), [0 lineLen*dAsp(1)], [0 0], [0 0], 'b', 'LineWidth',2)
            plot3(st.axModeShapes(ll,kk), [0 0], [0 lineLen*dAsp(2)], [0 0], 'g', 'LineWidth',2)
            plot3(st.axModeShapes(ll,kk), [0 0],[0 0], [0 lineLen*dAsp(3)], 'r', 'LineWidth',2)

            % add grid points
            hpltMS{1,ll,kk} = plot3(st.axModeShapes(1,1),x,y,z,'ok', ...
                'MarkerSize', 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
            % add DOF number to datatips
            row = dataTipTextRow("DOF#",(1:numel(x)));
            hpltMS{1,ll,kk}.DataTipTemplate.DataTipRows(end+1) = row;
            
        end

        



        fprintf(' done.\n')
        
        % store handle to the plot elements into the GUI
        st.hPlotModeShape = hpltMS ;

        guidata(uihandle,st) ;
    end % showModeShape3D


end % main function

