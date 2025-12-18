function model3D = createMeshObj(modelPath, Opt)
% Input:
% modelPath : path to the wavefont object containing the mesh
% Opt : structure including fields
% selectedPoints : original mesh vertices included in the
%                  analysis
% frfModelLink : array linking modeshapes amplitude to vextex normals
% Faces2lines : allows to modify the mesh by removing unwanted faces
%               and edges
% Output:
% hGui : structure containing the mesh informations such as the vertices,
%        vertex normals, edges, faces and the field frfModelLink
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


    % load the mesh of the structure (file extension .obj (wavefont object))
    Opt.Vertex_normals = 'all' ;
    % addpath('loadawobj2016/')
    model3D = loadawobj(modelPath, Opt) ;
    
    [~, reindexer] = ismember(1:size(model3D.v, 2), Opt.selectedPoints) ;
        
    % find fields of model3D listing faces
    fields = fieldnames(model3D) ;
    fPattern = "f" + digitsPattern ;
    fFields = extract([fields{:}], fPattern) ;
    if isfield(model3D, 'l') ; fFields{end+1} = 'l' ; end
    
    % convert faces and edges including unused measurement points
    for iii=1:length(fFields)
        if fFields{iii} == 'l'
            VertCount = 2 ;
        else
            VertCount = str2double(extract(fFields{iii}, digitsPattern)) ;
        end
        
        VrtxCtnPerFace = sum(ismember(model3D.(fFields{iii}), Opt.selectedPoints), 1) ;
        [sortedVrtxCtnPerFace, isort] = sort(sum(ismember(model3D.(fFields{iii}), Opt.selectedPoints), 1)) ;
        [UniqueVrtxCtnPerFace, i_unique, ~] = unique(sortedVrtxCtnPerFace) ;
        i_unique = i_unique(UniqueVrtxCtnPerFace >= 2) ;
        i_unique(end+1) = length(VrtxCtnPerFace) ;
        UniqueVrtxCtnPerFace = UniqueVrtxCtnPerFace(UniqueVrtxCtnPerFace >= 2) ;
        
        for ii = 1:length(i_unique)-1
            % get faces to convert 
            cnvrtFaces = model3D.(fFields{iii})(:, isort(i_unique(ii):(i_unique(ii+1)-1))) ;
            % convert them by excluding undesired points
            cnvrtFaces = reshape(cnvrtFaces(ismember(cnvrtFaces, Opt.selectedPoints)), [], size(cnvrtFaces, 2)) ;
                
            if UniqueVrtxCtnPerFace(ii) == 2
                fName = 'l' ;
            else
                fName = ['f' num2str(UniqueVrtxCtnPerFace(ii))] ; 
            end
            
            % create a field in model3D (if necessary)
            if ~isfield(model3D, fName) ; fFields{end+1} = fName ; end
            
            % store converted faces
            model3D.(fName)(:, end+(1:(i_unique(ii+1) - i_unique(ii)))) = cnvrtFaces ;
        end
        
        LogicalPos = ~(VrtxCtnPerFace == VertCount) ;
        
        model3D.(fFields{iii})(:, LogicalPos) = [] ;
        model3D.(fFields{iii}) = reshape(reindexer(model3D.(fFields{iii})), VertCount, []) ;
    end
    
    if isfield(model3D, 'l') ; fFields = fFields(1:end-1) ; end
    
    % regrouping all faces under one field
    nbVertFaces = cellfun(@(x) size(model3D.(x)).', fFields, 'UniformOutput', false) ;
    nbVertFaces = [nbVertFaces{:}] ;
    maxnbVert = max(prod(nbVertFaces)./nbVertFaces(2,:)) ;
    model3D.ftotal = [];
    if ~isnan(maxnbVert)
        for iii = 1:length(fFields)
            model3D.ftotal = [model3D.ftotal [model3D.(fFields{iii});...
                                                        NaN(maxnbVert-nbVertFaces(1,iii), nbVertFaces(2,iii))]] ;
        end
    end
    
    % remove vertices (and their normals) corresponding to unused meaasurement
    model3D.v = model3D.v(:, Opt.selectedPoints) ;
    model3D.vn = model3D.vn(:, Opt.selectedPoints, :) ;
    model3D.vn(:, ~Opt.frfModelLink) = NaN ;

    if any(~Opt.frfModelLink, 'all')
        model3D.vn(:, ~Opt.frfModelLink) = [] ;
    else
        model3D.vn = reshape(model3D.vn, 3, []) ;
    end
    
    % define normals associated to dof indices
    [model3D.frfModelLink(:,1), model3D.frfModelLink(:,2)] = find(Opt.frfModelLink);
end