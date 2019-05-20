function [mouse, T] = getMouse( mainDir, setMouse, varargin ) % mouse, comp, 
% getComp sets up mouse structure arrays and gathers/organizes comparison scoring files.
Nmouse = numel( setMouse );
IP = inputParser;
addRequired( IP, 'mainDir', @ischar )
addRequired( IP, 'setMouse', @iscell )
addParameter( IP, 'setj', flip(1:Nmouse), @isnumeric )
addParameter( IP, 'meta', false, @islogical )
addParameter( IP, 'astro', false, @islogical )
parse( IP, mainDir, setMouse, varargin{:} ); % mouse, comp, 
setj = IP.Results.setj;
getMeta = IP.Results.meta;
astroToggle = IP.Results.astro;
if ispc, slm = '\'; else, slm = '/'; end
fprintf('\ngetMouse... ');
mouse = struct; T = cell(1,Nmouse);
% Default imaging parameters
dist20x = [585.6, 585.6]; % dimensions of FOV for 20x objective with 1x digital zoom
res = [512,512]; digitalZoom = 8; stepDist = 0.5; 
pixDist = dist20x./res/digitalZoom; % um/pixel
% Find and read the index, containing info about each mouse
indexPath = [mainDir,'MouseIndex.xlsx'];
[~,~,mouseIndex] = xlsread(indexPath); fprintf('Read %s\n', indexPath);
varsRow = 1;
mouseVars = mouseIndex(varsRow,:);
expCol = find(strcmpi( mouseVars, 'Experiment')); idCol = find(strcmpi( mouseVars, 'ID' )); %cageCol = find(strcmpi( mouseVars, 'Cage' )); 
mouseIndex = mouseIndex( find( cellfun( @ischar, mouseIndex(:,expCol), 'UniformOutput', true ) ), :);  % disp(mouseIndex) 
% Find and read the sheet containing morphological info about each neuron
[~,~,morphIndex]= xlsread([mainDir,'MouseIndex.xlsx'], 'Morphology'); 
morphIndex(:,1) = cellfun( @num2str, morphIndex(:,1), 'UniformOutput',false);
morphVars = morphIndex(1,:); 
cellCol = find(strcmpi( morphVars, 'Cell' )); %dobCol = find(strcmpi( morphVars, 'DOB' ));
mbCol = find(strcmpi( morphVars, 'MB (um)' )); somaCol = find(strcmpi( morphVars, 'Soma (um)' ));  layerCol = find(strcmpi( morphVars, 'Layer' ));
% Set up mouse structure array and gather comparisons
IndCol = 5; NastroCol = 14; NneurCol = 15; NcellCol = 16; IsoCol = 17; BurnCol = 18; StimCol = 19; DrugCol = 20; %NoteCol =  % xCol = 6; dCol = 7; tCol = 8;
if isempty( fieldnames( mouse ) )
    mouse = repmat( struct('ID','', 'sex','', 'genotype',[], 'DOB',[], 'DOW',[], 'ages',[], 'dir','', 'figDir','', 'tiles','', 'Ntiles',NaN, 'tilename',[], ...
        'Nneur',NaN, 'Nastro',NaN,'Ncell',0,'totNeur',NaN,'totAstro',NaN,...
        'days',[], 'Ndays',NaN, 'daysDir',[], 'dates',[], 'datesStr', [], 'maxt',NaN, 'stacks',[], 'Nstacks',[], 'tif',[], 'xml',[], 'res',[], 'Nframe',[], 'dur',[], ... 
        'umPix',[], 'umStep',[], 'conv',[], ...
        'iso',NaN,'stim',NaN,'drug',[],...
        'neurTiles', [], 'NneurTiles',[], 'neurTilesInd', [],'neurSoma',[],'neurMB',[],'neurType',[],'daysComp',[], 'tracing',[], 'Ntracing',[],...   
        'filt',[],'astro',[],'Ca',[]...
        ), Nmouse,1);
end
tic;
for j = setj % Nmouse:-1:1
    fprintf('\n j = %d \n', j);
    % Get biographical info about each mouse from MouseIndex.xlsx
    indexLine(j) = find( strcmpi( cellfun(@num2str, mouseIndex(:,idCol), 'UniformOutput', 0), setMouse{j} ) );  %#ok<AGROW>
    mouse(j).ID = setMouse{j};
    mouse(j).sex = mouseIndex{indexLine(j),5};
    mouse(j).genotype = mouseIndex( indexLine(j), 7:11 ); % { strain, Cre+/-, Rosa, GFP+/-, Smo }
    if ispc
        mouse(j).DOB = datenum( mouseIndex{indexLine(j), 3}, 'mm/dd/yyyy' ); %x2mdate( mouseIndex{indexLine(j), 3} ); %  %datestr( mouse(j).DOB )  
    else
        mouse(j).DOB = x2mdate( mouseIndex{indexLine(j), 3} ); %  datenum( mouseIndex{indexLine(j), 3}, 'mm/dd/yyyy' ); %datestr( mouse(j).DOB )  
    end
    mouse(j).trim = mouseIndex{ indexLine(j), 13 };
    mouse(j).notes = mouseIndex{ indexLine(j), 24};
    % Which days and tiles have been imaged?
    mouse(j).dir = [mainDir,setMouse{j},slm]; 
    mouse(j).figDir = [mouse(j).dir,'Figures',slm]; mkdir( mouse(j).figDir );
    mouseCont = dir( mouse(j).dir );  mouseCont = mouseCont( ~strncmpi('.', {mouseCont.name}, 1) & [mouseCont.isdir] ); % mouseCont.name
    daysCont = mouseCont( strncmpi('d', {mouseCont.name}, 1 ) ); mouse(j).Ndays = numel( daysCont ); % daysCont.name
    daysDate = zeros( mouse(j).Ndays, 2 ); dateString = cell(mouse(j).Ndays,1); 
    tileCont = cell(mouse(j).Ndays,1); tileList = cell(1,50); %daysTile = cell(mouse(j).Ndays,1); 
    maxPre = 1;
    for d = 1:mouse(j).Ndays
        spacePos = strfind( daysCont(d).name, ' ' ); 
        dateString{d} = daysCont(d).name(spacePos+1:end);
        daysDate(d,:) = [ str2double( daysCont(d).name(2:spacePos-1) ), datenum( dateString{d} ) ]; 
        mouse(j).daysDir{d,1} = [mouse(j).dir,daysCont(d).name,slm];
        tileCont{d} = dir( [mouse(j).dir,daysCont(d).name] ); % tileCont{d}.name
        tileCont{d} = tileCont{d}( strncmpi([setMouse{j},'-d'],{tileCont{d}.name},numel([setMouse{j},'-d'])) & cellfun(@isempty, strfind( {tileCont{d}.name}, 'Map' )) & cellfun(@isempty, strfind( {tileCont{d}.name}, 'Cell' )) & cellfun(@isempty, strfind( {tileCont{d}.name}, 'Glia' )) ); % tileCont{d}.name
        for c = numel( tileCont{d} ):-1:1
            undPos = strfind( tileCont{d}(c).name, '_' ); dotPos = strfind( tileCont{d}(c).name, '.' ); %  
            tileName = tileCont{d}(c).name(undPos+1:dotPos-1);
            [letStart, letEnd] = regexp( tileName, '\D*'); % find the letters within the tileName
            % numbers BEFORE the tile name (prefix) indicate specific timepoints for repeated imaging of that tile on the same day
            if letStart(1) == 1
                t = 1;
            else
                t = str2double( tileName(1:letStart(1)-1) );
            end
            maxPre = max( t, maxPre );
            % numbers AFTER the tile name (suffix) indicate specific tiles from different neurons from the same mouse
            if letEnd(end) == numel(tileName)
                postNum = 1;
            else
                postNum = str2double( tileName(letEnd(end)+1:end) );
            end
            tileList{postNum} = [tileList{postNum},{tileName(letStart(1):end)}];
        end
    end
    tileList = cellfun( @unique, tileList, 'UniformOutput',false); 
    tileList( cellfun(@isempty, tileList) ) = [];
    mouse(j).tiles = [tileList{:}]; % final master list of all the tiles from this mouse
    mouse(j).Ntiles = numel( mouse(j).tiles ); % # of distinct tiles 
    % Sort tiles by which neuron they're associated with and get mophological info from MouseIndex (second sheet)
    v = 0;
    neurInd = find( strcmpi(mouse(j).ID,morphIndex(:,idCol)) )'; cellInd = cell2mat( morphIndex(neurInd,cellCol) )'; neurInd = neurInd(cellInd(~isnan(cellInd))); % ensure ordering of the cells is correct when reading morphology
    mouse(j).Nneur = numel(neurInd); % numel( tileList ); % # of neurs
    for n = 1:mouse(j).Nneur
        try
            mouse(j).neurTiles{n} = tileList{n};
            mouse(j).NneurTiles(n) = numel( mouse(j).neurTiles{n} );
            mouse(j).neurTilesInd{n} = v+1:v+mouse(j).NneurTiles(n);
            v = v + mouse(j).NneurTiles(n);
        catch
            warning('\nj = %d: %s-%d, found no matching tiles for this cell', j, mouse(j).ID, n );  
        end
        mouse(j).neurMB(n) = morphIndex{neurInd(n),mbCol};
        mouse(j).neurSoma(n) = morphIndex{neurInd(n),somaCol};
        if isnan(morphIndex{neurInd(n),layerCol}), mouse(j).neurType{n} = ''; else, mouse(j).neurType{n} = morphIndex{neurInd(n),layerCol}; end
    end
    % Figure out pertinent dates/ages
    [~,sortInd] = sort( daysDate(:,1), 'ascend' ); daysDate = daysDate(sortInd,:);
    mouse(j).datesStr = dateString( sortInd )'; %daysTile = daysTile(sortInd); % ignore tiles not used for spine scoring
    mouse(j).daysDir = mouse(j).daysDir(sortInd,:);
    mouse(j).DOW = daysDate(1,2)-daysDate(1,1); % datestr( mouse(j).DOW ) 
    mouse(j).days = daysDate(:,1)'; 
    mouse(j).dates = daysDate(:,2)';
    mouse(j).ages = mouse(j).DOW-mouse(j).DOB + mouse(j).days;
    % Gather and organize metadata about z-stacks (filepaths, resolution, # of frames)
    mouse(j).stacks = cell( mouse(j).Ntiles, mouse(j).Ndays, maxPre ); mouse(j).xml = cell( mouse(j).Ntiles, mouse(j).Ndays, maxPre ); mouse(j).filt = cell( mouse(j).Ntiles, mouse(j).Ndays, maxPre ); mouse(j).tilename = cell( mouse(j).Ntiles, mouse(j).Ndays, maxPre );
    mouse(j).absTime = nan( mouse(j).Ntiles, mouse(j).Ndays, maxPre);
    T{j} = cell( mouse(j).Ntiles, mouse(j).Ndays, maxPre );
    for d = 1:mouse(j).Ndays
        for c = numel( tileCont{d} ):-1:1 % tileCont{d}.name
            undPos = strfind( tileCont{d}(c).name, '_'); dotPos = strfind( tileCont{d}(c).name, '.'); %  tileCont{d}(c).name
            fullName = tileCont{d}(c).name(1:dotPos-1);
            tileName = tileCont{d}(c).name(undPos+1:dotPos-1);
            [letStart, ~] = regexp( tileName, '\D*'); % find the letters within the tileName
            % numbers before the tile name indicate specific timepoints for repeated imaging of that tile on the same day
            if letStart(1) == 1,  t = 1;  else,  t = str2double( tileName(1:letStart(1)-1) ); end
            x = find( strcmpi( tileName(letStart(1):end), mouse(j).tiles ) );
            %fprintf( '\n%s-d%s_%s.tif', mouse(j).ID, num2str(mouse(j).days(d),'%02.0f'), tileName )
            mouse(j).tilename{x,d,t} = tileName;
            mouse(j).stacks{x,d,t} = sprintf('%s-d%s_%s', mouse(j).ID, num2str(mouse(j).days(d),'%02.0f'), tileName );
            mouse(j).tif{x,d,t} = sprintf('%s.tif', mouse(j).stacks{x,d,t} );
            tifInfo = imfinfo( [mouse(j).daysDir{d},mouse(j).stacks{x,d,t},'.tif'] );
            %mouse(j).res(x,d,t) = tifInfo(1).Width;
            mouse(j).Nframe(x,d,t) = numel(tifInfo);
            metaFile = sprintf( '%sMetadata%s%s.xml', mouse(j).daysDir{d}, slm, fullName );
            mouse(j).umStep(x,d,t) = 0.5;
            if getMeta && exist( metaFile, 'file' ) == 2
                %fprintf('\nReading %s', metaFile);
                mouse(j).xml{x,d,t} = metaFile;
                metadata = xml2struct( metaFile );
                metaTime = metadata.PVScan.Attributes.date;
                spacePos = strfind(metaTime, ' ');
                mouse(j).time{x,d,t} = metaTime(spacePos(1)+1:end);
                mouse(j).absTime(x,d,t) = datenum( metaTime );
                %mouse(j).bit(x,d,t) = str2double( metadata.PVScan.PVStateShard.PVStateValue{2}.Attributes.value);
                mouse(j).dwell(x,d,t) = str2double( metadata.PVScan.PVStateShard.PVStateValue{4}.Attributes.value);
                mouse(j).framePeriod(x,d,t) = str2double( metadata.PVScan.PVStateShard.PVStateValue{5}.Attributes.value);
                mouse(j).pockels(x,d,t) = str2double( metadata.PVScan.PVStateShard.PVStateValue{11}.IndexedValue.Attributes.value);
                mouse(j).wavelength(x,d,t) = str2double( metadata.PVScan.PVStateShard.PVStateValue{14}.IndexedValue.Attributes.value);
                mouse(j).res(x,d,t,2) = str2double( metadata.PVScan.PVStateShard.PVStateValue{23}.Attributes.value);
                mouse(j).res(x,d,t,1) = str2double( metadata.PVScan.PVStateShard.PVStateValue{15}.Attributes.value);
                mouse(j).digitalZoom(x,d,t) = str2double( metadata.PVScan.PVStateShard.PVStateValue{22}.Attributes.value);
                mouse(j).umPix(x,d,t) = str2double( metadata.PVScan.PVStateShard.PVStateValue{17}.IndexedValue{1}.Attributes.value );
                mouse(j).conv{x,d,t} = [mouse(j).umPix(x,d,t), mouse(j).umPix(x,d,t), mouse(j).umStep(x,d,t)];
                mouse(j).PMT(x,d,t) = str2double( metadata.PVScan.PVStateShard.PVStateValue{24}.IndexedValue{2}.Attributes.value);
                mouse(j).pos(x,d,t,3) = str2double( metadata.PVScan.PVStateShard.PVStateValue{25}.SubindexedValues{3}.SubindexedValue.Attributes.value ); 
                mouse(j).pos(x,d,t,2) = str2double( metadata.PVScan.PVStateShard.PVStateValue{25}.SubindexedValues{1}.SubindexedValue.Attributes.value); 
                mouse(j).pos(x,d,t,1) = str2double( metadata.PVScan.PVStateShard.PVStateValue{25}.SubindexedValues{2}.SubindexedValue.Attributes.value);
                %preampGain = metadata.PVScan.PVStateShard.PVStateValue{27}.SubindexedValues.SubindexedValue{2}.Attributes; % double check this one
                %preampOffset = metadata.PVScan.PVStateShard.PVStateValue{28}.SubindexedValues.SubindexedValue{2}.Attributes; % double check this one
                mouse(j).laserPower(x,d,t) = str2double( metadata.PVScan.PVStateShard.PVStateValue{32}.IndexedValue.Attributes.value);
                %mouse(j).Nframe(x,d,t) =  numel(metadata.PVScan.Sequence.Frame);
                try
                    for z = mouse(j).Nframe(x,d,t):-1:1
                        T{j}{x,d,t}(z) = str2double( metadata.PVScan.Sequence.Frame{z}.Attributes.relativeTime );
                    end
                    mouse(j).dur(x,d,t) = T{j}{x,d,t}(end);
                    mouse(j).frameRate(x,d,t) = mouse(j).Nframe(x,d,t)/mouse(j).dur(x,d,t);
                catch
                    warning('Failed to get time vectors and framerate for %s', metaFile );
                end
            elseif getMeta %exist( metaFile, 'file' ) ~= 2
                mouse(j).res(x,d,t,1:2) = res;
                mouse(j).digitalZoom(x,d,t) = digitalZoom;
                mouse(j).umPix(x,d,t) = pixDist(1);
                mouse(j).umStep(x,d,t) = stepDist;
                mouse(j).conv{x,d,t} = [pixDist, stepDist];
                warning('Metadata file %s does not exist! ', metaFile );
            end
        end
    end
    mouse(j).Nstacks = sum( ~cellfun( @isempty, mouse(j).stacks ), 3 );  mouse(j).maxt = max( mouse(j).Nstacks(:) );
    % Read data from MovieTables, if available - NOTE: This will re-define Nneur
    if astroToggle
        movieIndexPath = sprintf('%s%s-MovieTable.xls', mouse(j).dir, mouse(j).ID );
        if exist(movieIndexPath,'file')
            [~,~,stackIndex] = xlsread( movieIndexPath ); fprintf('Read %s\n', movieIndexPath); % disp(stackIndex)
            stackIndex = stackIndex(2:end,:);
            mouse(j).Nastro = nan( mouse(j).Ntiles, mouse(j).Ndays, mouse(j).maxt ); mouse(j).Nneur = nan( mouse(j).Ntiles, mouse(j).Ndays, mouse(j).maxt ); mouse(j).Ncell = nan( mouse(j).Ntiles, mouse(j).Ndays, mouse(j).maxt );
            mouse(j).iso = nan( mouse(j).Ntiles, mouse(j).Ndays, mouse(j).maxt );
            mouse(j).stim = nan( mouse(j).Ntiles, mouse(j).Ndays, mouse(j).maxt );  mouse(j).drug = cell( mouse(j).Ntiles, mouse(j).Ndays, mouse(j).maxt );
            for q = 1:size(stackIndex,1)
                %stackIndex(q,:)
                %x = stackIndex{q,xCol}; d = stackIndex{q,dCol}; t = stackIndex{q,tCol};   
                Ind = stackIndex{q,IndCol};
                mouse(j).Nastro(Ind) = stackIndex{q,NastroCol}; 
                mouse(j).Nneur(Ind) = stackIndex{q,NneurCol};
                mouse(j).Ncell(Ind) = stackIndex{q,NcellCol};
                mouse(j).iso(Ind) = stackIndex{q,IsoCol};
                mouse(j).stim(Ind) = stackIndex{q,StimCol};
                if ischar(stackIndex{q,DrugCol})
                    mouse(j).drug{Ind} = stackIndex(q,DrugCol);
                else  
                    mouse(j).drug{Ind} = '';
                end
            end
            mouse(j).totAstro = nansum( max(max( mouse(j).Nastro, [], 3),[],2), 1);
            mouse(j).totNeur = nansum( max(max( mouse(j).Nneur, [], 3),[],2), 1);
        end
    end
    %}
    disp( mouse(j) )
end
fprintf('\n');
toc 
end