function [ mouse, comp, Ncomp ] = getComp( mouse, varargin ) % mouse, comp, 
% getComp sets up mouse structure arrays and gathers/organizes comparison scoring files.
Nmouse = numel( mouse );
IP = inputParser;
addRequired( IP, 'mouse', @isstruct )
addParameter( IP, 'setj', flip(1:Nmouse), @isnumeric )
parse( IP, mouse, varargin{:} ); % mouse, comp, 
setj = IP.Results.setj;
if ispc, slm = '\'; else, slm = '/'; end
fprintf('\ngetComp... ');
comp = cell(1,Nmouse); 
compStruct = struct('mouse','', 'ages',[NaN,NaN], 'days',[NaN,NaN], 'd',[NaN,NaN], 'date',[],'interval',NaN, 'timepoints', [NaN,NaN], 'scoring',[], 'Nscoring',[0,0], 'tiles',[], 'Ntiles',0, 'x',[], 'totFeat',zeros(1,7), ...
    'Nfeat',[], 'feat',[], 'um',[], 'match',[], 'Nprot',0, 'Nspine',0, 'spineFrac', NaN, 'filoFrac',NaN, 'ambigFrac', NaN, 'SF',NaN, 'TO',NaN, 'elim',NaN, 'form',NaN, ...
    'neurNfeat',[], 'neurTotFeat', [], 'neurNprot',[], 'neurNspine', [], 'neurSpineFrac',[], 'neurFiloFrac',[], 'neurAmbigFrac',[], 'neurSF',[], 'neurTO',[], 'neurElim',[], ...
    'neurForm',[], 'neurFelim',[], 'neurFform',[] ); 
tic;
for j = setj % Nmouse:-1:1 % 
    %fprintf('\n j = %d ', j);
    % Find and load scoring (.xml) files and get basic details about each comparison.
    scoreDir = [mouse(j).dir,'Scoring',slm]; %fprintf('\n%s\n', scoreDir )
    scoreCont = dir( scoreDir ); scoreCont = scoreCont( ~strncmpi('.', {scoreCont.name}, 1) & ~[scoreCont.isdir] ); % scoreCont.name
    compStruct.scoring = cell(mouse(j).Ntiles, 2); %compStruct.stack = cell(mouse(j).Ntiles, 2);
    comp{j} = repmat( compStruct, mouse(j).Ndays, mouse(j).Ndays); %repmat( compStruct, daysDate(end,1), daysDate(end,1) ); 
    linInd = [];
    for k = 1:numel( scoreCont ) % {scoreCont.name}
        T = [1, 1]; 
        charPos = getCharPos( scoreCont(k).name, {'-','_','.'} ); % scoreCont(k).name
        compName = scoreCont(k).name(charPos{1}(1)+1:charPos{3}(1)-1);
        charPos = getCharPos( compName, {'d','_'} );
        % Which tile was compared?
        tileStr = compName( charPos{2}(1)+1:end );
        dashPos = strfind(tileStr,'-');
        if isempty( dashPos )
            tileName = tileStr;
        else % If D(1) == D(2), which timepoints were compared?
            Tstr{1} = tileStr(1:dashPos-1); Tstr{2} = tileStr(dashPos+1:end);
            [letStart, ~] = regexp( Tstr{1}, '\D*'); %letEnd
            tileName = Tstr{1}(letStart:end); %
            T(1) = str2double( Tstr{1}(1:letStart-1) );
            [letStart, ~] = regexp( tileStr(dashPos+1:end), '\D*');
            T(2) = str2double( Tstr{2}(1:letStart-1) );
        end
        X = find( strcmpi( tileName, mouse(j).tiles ) );
        if isempty(X), error('Invalid tile name: %s', tileName); end
        % Which days were compared?
        dayStr = compName( charPos{1}(1)+1:charPos{2}(1)-1 );
        dashPos = strfind(dayStr,'-');
        if isempty( dashPos )
            D = str2double( dayStr )*[1,1];
            [~,~,d] = intersect( D, mouse(j).days ); %d = d'; % indices corresponding to those days
            d = d*[1,1];
        else
            D = [str2double( dayStr(1:dashPos-1)), str2double( dayStr(dashPos+1:end))];
            [~,~,d] = intersect( D, mouse(j).days ); %d = d'; % indices corresponding to those days
        end
        % Populate the comp structure arrays
        if diff(D) < 0
            d = sort( d ); % keep forward and reverse scoring files together - MAYBE UNNECESSARY
            comp{j}(d(1), d(2)).Nscoring(2) = comp{j}(d(1), d(2)).Nscoring(2) + 1;
            comp{j}(d(1), d(2)).scoring{X,2} = [scoreDir,scoreCont(k).name];
        elseif diff(D) > 0
            comp{j}(d(1), d(2)).Nscoring(1) = comp{j}(d(1), d(2)).Nscoring(1) + 1;
            comp{j}(d(1), d(2)).scoring{X,1} = [scoreDir,scoreCont(k).name];
            comp{j}(d(1), d(2)).Ntiles = comp{j}(d(1), d(2)).Ntiles + 1;
            comp{j}(d(1), d(2)).tiles{comp{j}(d(1), d(2)).Ntiles, 1} = tileName; % comp{j}(d(1), d(2)).tiles{X,1} = tileName;
        elseif diff(D) == 0 && diff(T) > 0
            comp{j}(d(1), d(2)).Nscoring(1) = comp{j}(d(1), d(2)).Nscoring(1) + 1;
            comp{j}(d(1), d(2)).scoring{X,1} = [scoreDir,scoreCont(k).name];
            comp{j}(d(1), d(2)).Ntiles = comp{j}(d(1), d(2)).Ntiles + 1;
            comp{j}(d(1), d(2)).tiles{comp{j}(d(1), d(2)).Ntiles, 1} = tileName; % comp{j}(d(1), d(2)).tiles{X,1} = tileName;
        elseif diff(D) == 0 && diff(T) < 0
            comp{j}(d(1), d(2)).Nscoring(2) = comp{j}(d(1), d(2)).Nscoring(2) + 1;
            comp{j}(d(1), d(2)).scoring{X,2} = [scoreDir,scoreCont(k).name];
        elseif diff(D) == 0 && diff(T) == 0
            comp{j}(d(1), d(1)).Nscoring(1) = comp{j}(d(1), d(2)).Nscoring(1) + 1;
            comp{j}(d(1), d(1)).scoring{X,1} = [scoreDir,scoreCont(k).name];
        end
        %{
        if sum( comp{j}(d(1), d(2)).Nscoring ) == 1 && diff(D) >= 0
            comp{j}(d(1), d(2)).mouse = mouse(j).ID; % 
            comp{j}(d(1), d(2)).days = D;
            comp{j}(d(1), d(2)).interval = diff(D);
            comp{j}(d(1), d(2)).timepoints = T;
            comp{j}(d(1), d(2)).ages = mouse(j).DOW - mouse(j).DOB + D; % ages at the 2 time points compared
            comp{j}(d(1), d(2)).date{2} = mouse(j).datesStr{d(2)};
            comp{j}(d(1), d(2)).date{1} = mouse(j).datesStr{d(1)};
        end
        %}
        linInd = [linInd, sub2ind( mouse(j).Ndays*[1,1], d(1), d(2) )]; %#ok<AGROW> % linear index corresponding to (day1,day2)
    end
    linInd = unique(linInd);
    comp{j} = comp{j}( linInd ); % convert from 2D to 1D array
    % Once all of the scoring files have been processed...
    mouse(j).daysComp = nan( mouse(j).Ndays, mouse(j).Ndays, mouse(j).Ntiles ); % index which comparisons cover which days and which tiles?
    for k = 1:numel( comp{j} )
        [d(1), d(2)] = ind2sub( [mouse(j).Ndays,mouse(j).Ndays], linInd(k) );
        comp{j}(k).mouse = mouse(j).ID; % 
        comp{j}(k).d = d';
        comp{j}(k).days = mouse(j).days(d);
        comp{j}(k).interval = diff(comp{j}(k).days);
        comp{j}(k).timepoints = T;
        comp{j}(k).ages = mouse(j).DOW - mouse(j).DOB + comp{j}(k).days; % ages at the 2 time points compared
        comp{j}(k).date{2} = mouse(j).datesStr{d(2)};
        comp{j}(k).date{1} = mouse(j).datesStr{d(1)};
        comp{j}(k).x = find(~cellfun( @isempty, comp{j}(k).scoring(:,1) )');
        for x = 1:comp{j}(k).Ntiles
            mouse(j).daysComp( find(comp{j}(k).days(1) == mouse(j).days), find(comp{j}(k).days(2) == mouse(j).days), find(strcmpi(comp{j}(k).tiles(x), mouse(j).tiles)) ) = k; %#ok<*FNDSB>
            mouse(j).daysComp( find(comp{j}(k).days(2) == mouse(j).days), find(comp{j}(k).days(1) == mouse(j).days), find(strcmpi(comp{j}(k).tiles(x), mouse(j).tiles)) ) = k;
        end
    end
end
Ncomp = cellfun( @numel, comp );
fprintf('\n');
toc 
end