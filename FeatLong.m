function [long, longMat, longMeta, longFeat ] = FeatLong( mouse, featParam, varargin ) % , intLong, volLong, sepLong
% FeatImage generates cropped summary images of longitudinally scored
% features
% -----Inputs------
% mouse 
% featParam
% minSNR = minimum SNR of a high-quality feature image
% setj
% show
% -----Outputs------ (same format repeats for axon outputs)
% long = tracking (x,y,z,click) of each feature at each timepoint
% longMat = matrix summaring the presence/absence of each feature at each timepoint
% longMeta = metadata about each longitudinal scoring file
% longFeat = cell array of structures containing summary image and various properties about each longitudinally scored feature 

Nmouse = numel( mouse );
tic;
IP = inputParser;
addRequired( IP, 'mouse', @isstruct )
addRequired( IP, 'featParam', @isstruct )
%addOptional( IP, 'dend', cell(1,Nmouse), @iscell );
%addParameter( IP, 'minSNR', 4, @isnumeric )
addParameter( IP, 'setj', flip(1:Nmouse), @isnumeric )
%addParameter( IP, 'show', false, @islogical )
parse( IP, mouse, featParam, varargin{:} );
%minSNR = IP.Results.minSNR;
setj = IP.Results.setj;
%show = IP.Results.show;
% show = false; 
fprintf('FeatLong...  ');
if ispc, slm = '\'; else, slm = '/'; end
longFcn = @(x)( ~contains( x, 'TEMP' ) &  ~contains( x, 'axon' ) ); %@(x)( ~contains(x,'TEMP') & ~contains(x,'axon') ); % prevent loading of temporary or axon-based files  %  
clearvars longMeta

long = cell(1,Nmouse); longMat = cell(1,Nmouse); longFeat = cell(1,Nmouse);
longMeta(Nmouse) = struct('x',[], 'Nstack',[], 'sTrack',[], 'Ntrack',[], 'k',[], 'd',[], 't',[], 'days',[], 'Ndays',[], 'daysMat',[], 'dUse',[], 'Nuse',[], 'useDays',[], 'useDur',[],...
    'seedID',[], 'Nlong',[], 'nanMat',[], 'Nnan',[], 'Nfull',[], 'qFull',[], 'qInc',[], 'qEmpty',[], 'Ninc',[], 'Nempty',[], 'qAll',[], ...
    'qType',[], 'Ntype',[], 'totType',[], 'qStable',[], 'Nstable',[], 'totStable',[], 'qDyn',[], 'Ndyn',[], 'totDyn',[], 'qRec',[], 'Nrec',[], 'totRec',[], 'qTrans',[], 'Ntrans',[], 'totTrans',[], ...
    'perm',[], 'imperm',[], 'tau',[], 'R2',[] ); %#ok<*AGROW>
if nargout > 3
    baseFeat = struct('source',[], 'ID',NaN, 'ind',nan(1,5), 'type',NaN, 'click',NaN, 'eye',nan(1,3), 'um',nan(1,3), 'nearDist',NaN, 'xy',nan(1,2), 'cropVec',nan(1,4), 'eyeCrop',nan(1,2), 'xyCrop',nan(1,2), 'cent',nan(1,2), 'raw',[], ...
        'filt',[], 'reg',[], 'meanShift',NaN, 'eyeInt',NaN, 'SNR',NaN, 'binary',[], 'cart',[], 'pol',[], 'edge',[], 'perim',NaN, 'Nrub',NaN, 'rubRatio',NaN, 'area',NaN, 'int',NaN, 'sat',nan(1,2), 'radius',NaN, 'vol',NaN, ...
        'dendInt',NaN, 'norm',NaN, 'seg',NaN, 'track',NaN, 'sep',NaN, 'long',NaN,...
        'shaftInt',NaN,'shaftCent',nan(1,2), 'shaftAx',nan(3,3), 'shaftArea',NaN, 'match',NaN, 'ambig',NaN, 'matchSep',NaN, 'matchErr',NaN, 'matchXYZ',nan(1,3), 'good',NaN, 'cut', '', 'd2s',NaN, 'angle',NaN);
end
tic;
for j = setj
    fprintf('\nj = %d / %d  ', j, Nmouse);
    longDir = [mouse(j).dir,'Longitudinal Analysis',slm]; 
    longFiles = FileFind( longDir, 'mat', false, longFcn );
    long{j} = cell(1,mouse(j).Ntiles); longMat{j} = cell(1,mouse(j).Ntiles); longFeat{j} = cell(1,mouse(j).Ntiles); 
    longMeta(j).k = nan(1,mouse(j).Ntiles); longMeta(j).Nlong = nan(1,mouse(j).Ntiles);
    longMeta(j).allType = cell( mouse(j).Ntiles,3); longMeta(j).Nall = nan(mouse(j).Ntiles,3);  %longMeta(j).totType = nan(1, 3);
    longMeta(j).qType = cell( mouse(j).Ntiles,3); longMeta(j).Ntype = nan(mouse(j).Ntiles,3);  longMeta(j).totType = nan(1, 3);
    longMeta(j).qStable = cell( mouse(j).Ntiles,3); longMeta(j).Nstable = nan(mouse(j).Ntiles, 3);  longMeta(j).totStable = nan(1, 3);
    longMeta(j).qDyn = cell( mouse(j).Ntiles,3); longMeta(j).Ndyn = nan(mouse(j).Ntiles,3);  longMeta(j).totDyn = nan(1,3);
    longMeta(j).qRec = cell( mouse(j).Ntiles,3); longMeta(j).Nrec = nan(mouse(j).Ntiles,3);  longMeta(j).totRec = nan(1,3);
    longMeta(j).qTrans = cell( mouse(j).Ntiles,3);  longMeta(j).Ntrans = nan(mouse(j).Ntiles,3);  longMeta(j).totTrans = nan(1, 3);
    longMeta(j).pop = nan(1,3);  longMeta(j).tau = nan(2,3); % longMeta(j).imperm = nan(1,3);
    if size( longFiles, 1 ) > 0
        for k = 1:size(longFiles,1)
            %fprintf('\n  k = %d', k); %clearvars nanMat 
            % Gather metadata about each longitudinal scoring
            L = load( longFiles{k,2} ); fprintf('\nLoaded %s\n', longFiles{k,2} ); % disp(L) % , loadVars{:}
            qBouton = reshape( L.qBouton, 1, numel(L.qBouton) ); qSpine = reshape( L.qSpine, 1, numel(L.qSpine) );  qFilo = reshape( L.qFilo, 1, numel(L.qFilo) );
            x = find( strcmpi(L.setTile, mouse(j).tiles) ); %X = S.X; % 
            longMeta(j).x = [longMeta(j).x, x];
            longMeta(j).k(x) = L.K;
            longMeta(j).Nlong(x) = size(L.long,1);
            longMeta(j).Nstack(x) = L.Nstack;
            longMeta(j).sTrack{x} = L.sTrack; %fprintf('sTrack = %d ', longMeta(j).sTrack{x} );
            longMeta(j).Ntrack(x) = numel( longMeta(j).sTrack{x} );
            longMeta(j).d{x} = L.DT(longMeta(j).sTrack{x},1);
            longMeta(j).days{x} = mouse(j).days(longMeta(j).d{x});
            longMeta(j).Ndays(x) = numel( longMeta(j).days{x} );
            longMeta(j).seedID{x} = L.seedID; % feature IDs corresponding to q indices
            longMeta(j).t{x} = L.DT(longMeta(j).sTrack{x},2)';
            %long{j}{X} = nan( longMeta(j).Nlong(X), 4, mouse(j).Ndays ); %
            long{j}{x} = L.long;  % (:,:,longMeta(j).sTrack{X})
            longMat{j}{x} = reshape( long{j}{x}(:,4,longMeta(j).sTrack{x}), longMeta(j).Nlong(x),  longMeta(j).Ntrack(x) ); % reshape( long{j}{X}(:,4,:), longMeta(j).Nlong(X), mouse(j).Ntrack(X) );
            nanMat = longMat{j}{x}; nanMat( ~isnan(nanMat) ) = 0; nanMat( isnan(nanMat) ) = 1;
            longMeta(j).nanMat{x} = nanMat;
            longMeta(j).Nnan{x} = sum(nanMat,2); % sum(nanMat(:,longMeta(j).d{x}),2);
            longMeta(j).qAll{x,3} = qFilo;  longMeta(j).qAll{x,2} = qSpine;  longMeta(j).qAll{x,1} = qBouton; % includes incompletely scored features
            longMeta(j).qFull{x} = find(longMeta(j).Nnan{x} == 0)'; % longMeta(j).Ndays(X) - longMeta(j).Ntrack(X) 
            longMeta(j).Nfull(x) = numel(longMeta(j).qFull{x});
            longMeta(j).qInc{x} = find( longMeta(j).Nnan{x} > 0 & longMeta(j).Nnan{x} < longMeta(j).Ntrack(x) )';
            longMeta(j).Ninc(x) = numel(longMeta(j).qInc{x});  
            longMeta(j).qEmpty{x} = find( longMeta(j).Nnan{x} >= longMeta(j).Ntrack(x) )';
            longMeta(j).Nempty(x) = numel(longMeta(j).qEmpty{x});
            % Get the indices of FULLY-SCORED boutons, spines and filopodia, based on their first appearance in the LONGITUDINAL scoring
            qFullProt = intersect( [qSpine, qFilo], longMeta(j).qFull{x} ); % fully-scored protrusions
            for q = qFullProt
                s1 = find( longMat{j}{x}(q,:) < 3, 1, 'first' ); % stack index when the feature was first observed
                if longMat{j}{x}(q,s1) == 1
                    longMeta(j).qType{x,2} = [longMeta(j).qType{x,2}, q];
                elseif longMat{j}{x}(q,s1) == 2
                    longMeta(j).qType{x,3} = [longMeta(j).qType{x,3}, q];
                end
            end
            longMeta(j).qType{x,1} = intersect( qBouton, longMeta(j).qFull{x} );
            longMeta(j).Ntype(x,3) = numel( longMeta(j).qType{x,3} );
            longMeta(j).Ntype(x,2) = numel( longMeta(j).qType{x,2} );
            longMeta(j).Ntype(x,1) = numel( longMeta(j).qType{x,1} );
            longMeta(j).qStable{x,3} = intersect( longMeta(j).qType{x,3},  find( nansum( longMat{j}{x}, 2 ) == longMeta(j).Ntrack(x) )' ); longMeta(j).Nstable(x,3) = numel(longMeta(j).qStable{x,3}); % stable filopodia
            longMeta(j).qStable{x,2} = intersect( longMeta(j).qType{x,2},  find( nansum( longMat{j}{x}, 2 ) == longMeta(j).Ntrack(x) )' ); longMeta(j).Nstable(x,2) = numel(longMeta(j).qStable{x,2});%Nss = numel(qSS); % stable spines
            longMeta(j).qStable{x,1} = intersect( longMeta(j).qType{x,1},  find( nansum( longMat{j}{x}, 2 ) == longMeta(j).Ntrack(x) )' ); longMeta(j).Nstable(x,1) = numel(longMeta(j).qStable{x,1});%
            longMeta(j).qDyn{x,3} = intersect( longMeta(j).qType{x,3},  find( nansum( longMat{j}{x}, 2 ) ~= longMeta(j).Ntrack(x) )' ); longMeta(j).Ndyn(x,3) = numel(longMeta(j).qDyn{x,3}); % dynamic filopodia
            longMeta(j).qDyn{x,2} = intersect( longMeta(j).qType{x,2},  find( nansum( longMat{j}{x}, 2 ) ~= longMeta(j).Ntrack(x) )' ); longMeta(j).Ndyn(x,2) = numel(longMeta(j).qDyn{x,2}); % dynamic spines
            longMeta(j).qDyn{x,1} = intersect( longMeta(j).qType{x,1},  find( nansum( longMat{j}{x}, 2 ) ~= longMeta(j).Ntrack(x) )' ); longMeta(j).Ndyn(x,1) = numel(longMeta(j).qDyn{x,1}); % dynamic boutons
            % Extract image properties for each feature/timepoint
            if nargout > 3
                longFeat{j}{x} = repmat( baseFeat, [longMeta(j).Nlong(x), longMeta(j).Ntrack(x)] );
                seed = L.seed;
                for s = 1:longMeta(j).Ntrack(x)
                    d = longMeta(j).d{x}(s); t = longMeta(j).t{x}(s);
                    stackPath = sprintf('%s%s.tif', mouse(j).daysDir{d}, mouse(j).stacks{x,d,t} );
                    rawStack = loadtiff( stackPath );  fprintf('\nLoaded %s', stackPath )
                    Nframe = size(rawStack,3); 
                    for q = 1:longMeta(j).Nlong(x)
                        %fprintf('\n[q,s] = [%d, %d]', q, s);
                        % Basic properties of the current feature on the current day 
                        try 
                            longFeat{j}{x}(q,s).source = seed(q).source;
                            longFeat{j}{x}(q,s).type = seed(q).ind(4);
                        catch 
                            longFeat{j}{x}(q,s).source = sprintf('%s-d%d-%d-%s', mouse(j).ID, mouse(j).days(1), mouse(j).days(end), mouse(j).tiles{x} ); 
                            longFeat{j}{x}(q,s).type = NaN;  %seed(q).ind(4);
                        end
                        longFeat{j}{x}(q,s).ind = [j,x,d,t,q];
                        longFeat{j}{x}(q,s).click = longMat{j}{x}(q,s);
                        if ~any( isnan(long{j}{x}(q,1:3,longMeta(j).sTrack{x}(s))) )
                            longFeat{j}{x}(q,s).eye = long{j}{x}(q,1:3,longMeta(j).sTrack{x}(s)); 
                            longFeat{j}{x}(q,s).um = featParam.conv.*longFeat{j}{x}(q,s).eye;
                            % Determine which frames to average over
                            ztemp = longFeat{j}{x}(q,s).eye(3)-featParam.Nneigh:longFeat{j}{x}(q,s).eye(3)+featParam.Nneigh;
                            ztemp(ztemp < 1) = []; ztemp(ztemp > Nframe) = []; Nztemp = numel(ztemp);
                            zref = find(ztemp == longFeat{j}{x}(q,s).eye(3));
                            tempStack = rawStack(:,:,ztemp);
                            clearvars rawCrop; % rawCrop may have different sizes, depending on whether cropped area falls completely on the image or not
                            if any( featParam.crop ~= 0 )
                                cropParam = [round(longFeat{j}{x}(q,s).eye(1:2) - featParam.crop/2), featParam.crop];  cropParam(cropParam < 1) = 1;
                                longFeat{j}{x}(q,s).cropVec = cropParam;
                                longFeat{j}{x}(q,s).eyeCrop = [longFeat{j}{x}(q,s).eye(1)-cropParam(1)+1, longFeat{j}{x}(q,s).eye(2)-cropParam(2)+1];
                                %longFeat{j}{X}(q,s).raw = imcrop( rawStack(:,:,longFeat{j}{X}(q,s).eye(3)), cropParam);
                                for z = Nztemp:-1:1
                                    rawCrop(:,:,z) = imcrop( tempStack(:,:,z), cropParam ); 
                                end
                            else
                                rawCrop = tempStack;
                                longFeat{j}{x}(q,s).eyeCrop = longFeat{j}{x}(q,s).eye(1:2);
                                %longFeat{j}{X}(q,s).raw = rawStack(:,:,longFeat{j}{X}(q,s).eye(3));
                            end
                            [~, longFeat{j}{x}(q,s).reg, longFeat{j}{x}(q,s).filt] = StackReg( rawCrop, zref ); % Average neighboring frames using registration of cropped images to generate summary images
                            longFeat{j}{x}(q,s).meanShift = nanmean( longFeat{j}{x}(q,s).reg(sum(longFeat{j}{x}(q,s).reg,2)>0, 3) );
                            % Feat = longFeat{j}{X}(q,s) % imshow( Feat.filt, [] )
                            %fprintf('\n[j,x,s,q] = [%d, %d, %d, %d]', j, X, s, q);
                            longFeat{j}{x}(q,s) = getFeatProp( longFeat{j}{x}(q,s), featParam ); % Extract feature properties    , show
                        else
                            longFeat{j}{x}(q,s).filt = zeros( featParam.crop + [1,1], 'uint16' );
                            longFeat{j}{x}(q,s).edge = cell(1,2);
                        end
                    end
                end
            end
        end
        longMeta(j).Nall = cellfun( @numel, longMeta(j).qAll );
        longMeta(j).totType = nansum( longMeta(j).Ntype );
        longMeta(j).totStable = nansum( longMeta(j).Nstable );
        longMeta(j).totDyn = nansum( longMeta(j).Ndyn );
        % which days should be used for the OVERALL survival calculations
        longMeta(j).daysMat = cell2padmat( longMeta(j).d( longMeta(j).x) )';
        [minDays,minInd] = min( sum( ~isnan(longMeta(j).daysMat ), 2 ) ); % set of days scored across most tiles
        longMeta(j).dUse = longMeta(j).daysMat(minInd,1:minDays); 
        longMeta(j).Nuse = numel(longMeta(j).dUse);
        longMeta(j).useDays = mouse(j).days( longMeta(j).dUse );
        longMeta(j).useDur = diff( longMeta(j).useDays([1,end]) );
    end
    toc
end
fprintf('\nFeatLong complete.\n');
end
