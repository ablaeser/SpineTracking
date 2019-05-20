function [ comp, feat, daysFeat ] = FeatRead( mouse, comp, featParam, varargin ) % feat, daysFeat, 
%FeatRead reads all the scoring files listed in comp structure, collects the coordinates of each feature in comp structures, and creates feat structs for each feature
% -----Inputs-------
% mouse = metadata about experiments perfoxrmed on each mouse
% comp = metadata about each comparison analyzed from each mouse
% featParam
%   .conv = conversion factors: [pixels, pixels, z-steps] -> micrometers
% -----Outputs------
% comp
% feat = structure containing information about each feature
%  .source = string indicating the stack from which the feature was identified
%  .ind = the indices that specify this feature [j,k,x,w,p,q]
%  .eye = coordinates (x,y,x) of the marker placed by eye
%  .ID = unique ID number from coordinates (3 digits each) and type: xyzw
%  .um = conversion of eye to microns
%  .nearDist = distance (um) to nearest neighboring feature of the same type
% daysFeat
Nmouse = numel( mouse ); 
IP = inputParser;
addRequired( IP, 'mouse', @isstruct )
addRequired( IP, 'comp', @iscell )
addRequired( IP, 'featParam', @isstruct )
addParameter( IP, 'setj', flip(1:Nmouse), @isnumeric )
parse( IP, mouse, comp, featParam, varargin{:} ); % feat, daysFeat, 
setj = IP.Results.setj;
fprintf('\nFeatRead...  ');
feat = cell(1,Nmouse); daysFeat = cell(1,Nmouse);
% Functions for measuring spine dynamics
SF = @(S,E)( S/(S+E) ); % stable fraction
TO = @(S,E,F)( (E+F)/(2*S+E+F) ); % turnover
% Initialize variables
baseFeat = struct('source',[], 'ind',[], 'type',[], 'eye',nan(1,3), 'ID',NaN, 'um',nan(1,3), 'nearDist',NaN, ... % fields filled in FeatRead
    'cropVec',nan(1,4), 'eyeCrop',nan(1,2), 'raw',[], 'reg',[], 'meanShift',NaN, 'filt',[],...  % fields filled in FeatImage 
    'xyCrop',nan(1,2), 'xy',nan(1,2), 'eyeInt',[], 'cent',[], 'radius',NaN, 'perim',NaN, 'area',NaN, 'vol',NaN, 'edge',[], 'int',NaN, 'sat',nan(1,2), 'norm',NaN,'d2s',NaN, 'angle',NaN, ... % fields filled in FeatProps 
    'binary',[], 'cart',[], 'pol',[], 'SNR',NaN, 'Nrub',NaN, 'rubRatio',NaN, 'shaftCent',nan(1,2), 'shaftAx',nan(3,3), 'shaftArea',NaN, 'shaftInt',NaN, ... % fields filled in FeatProps 
    'match',NaN, 'matchXYZ',nan(1,3), 'ambig',NaN, 'matchSep',NaN, 'matchErr',NaN, ... % fields filled in FeatMatch    
    'sep',NaN, 'seg',NaN, 'segPt',[], ...  % fields filled in DendRead
    'good',NaN, 'cut', '', ... % fields filled in FeatCuts
    'long',[] ... % fields filled in FeatLong
    ); 
Ncomp = cellfun(@numel, comp);

tic;
for j = setj
    fprintf('\n j = %d \n', j);
    if any( cellfun( @diff, {comp{j}.Nscoring} ) ~= 0 )
        error('Unmatched scoring file(s)'); % make sure all comparisons have forward and backward scoring files
    end
    for k = 1:Ncomp(j)
        %fprintf('\nj = %d, k = %d', j, k);
        featMat = cell(mouse(j).Ntiles,featParam.Ntype,2); feat{j}{k} = cell(mouse(j).Ntiles,featParam.Ntype,2); Nfeat = nan(mouse(j).Ntiles,featParam.Ntype,2);  %#ok<*AGROW>
        %xtemp = flip( find( ~cellfun( @isempty, comp{j}(k).scoring(:,1) ) )' );
        for x = comp{j}(k).x % xtemp
            for p = [2,1]
                [featMat(x,:,p), Nfeat(x,:,p)] = getFeatCoord( comp{j}(k).scoring{x,p}, featParam.Ntype ); % XMLpath = comp{j}(k).scoring{x,p};
                for w = find( Nfeat(x,:,p) )
                    feat{j}{k}{x,w,p} = repmat( baseFeat, Nfeat(x,w,p), 1 );
                    umMat = repmat(mouse(j).conv{x,comp{j}(k).d(1),1},Nfeat(x,w,p),1).*(featMat{x,w,p} - repmat([0,0,1],Nfeat(x,w,p),1)); % correct for fact that first frame is z = 1, not z = 0 and convert to um
                    %umMat = repmat(featParam.conv,Nfeat(x,w,p),1).*(featMat{x,w,p} - repmat([0,0,1],Nfeat(x,w,p),1)); % correct for fact that first frame is z = 1, not z = 0 and convert to um
                    distMat = squareform( pdist( umMat ) ); distMat( distMat == 0 ) = NaN;
                    for q = Nfeat(x,w,p):-1:1
                        feat{j}{k}{x,w,p}(q).source = sprintf('%s-d%d-%d-%s', mouse(j).ID, comp{j}(k).days(1), comp{j}(k).days(2), mouse(j).tiles{x});
                        feat{j}{k}{x,w,p}(q).ind = [j,k,x,w,p,q];
                        feat{j}{k}{x,w,p}(q).type = w;
                        feat{j}{k}{x,w,p}(q).eye = featMat{x,w,p}(q,:);
                        feat{j}{k}{x,w,p}(q).ID = str2double( sprintf('%03d%03d%03d%01d',feat{j}{k}{x,w,p}(q).eye(1),feat{j}{k}{x,w,p}(q).eye(2),feat{j}{k}{x,w,p}(q).eye(3),w) ); % generate unique ID number from coordinates (3 digits each) and type: xyzw
                        feat{j}{k}{x,w,p}(q).um = mouse(j).conv{x,comp{j}(k).d(1),1}.*(feat{j}{k}{x,w,p}(q).eye - [0,0,1]); 
                        if Nfeat(x,w,p) > 1, feat{j}{k}{x,w,p}(q).nearDist = min( distMat(q,:) ); end
                    end
                end
            end
        end
        % Save featMat to comp structures, compute various quantities and perform feature matching
        comp{j}(k).Nfeat = Nfeat;
        comp{j}(k).feat = featMat; 
        comp{j}(k).um = feat2um( featMat, mouse(j).conv{x,comp{j}(k).d(1),1} ); % convert to um
        comp{j}(k).totFeat = nansum( Nfeat(:,:,1), 1 );
        comp{j}(k).totFeat(4) = comp{j}(k).totFeat(4) + nansum( Nfeat(:,4,2), 1 ); % filopodia are marked in one or other scoring file, but not both
        comp{j}(k).Nprot = nansum( comp{j}(k).totFeat(1:5) ); % [1,2,4,5]
        comp{j}(k).Nspine = nansum( comp{j}(k).totFeat(1:3) );
        % Relative proportions of protrusion types on FIRST time point 
        comp{j}(k).spineFrac = nansum(comp{j}(k).totFeat([1,2,5]))/nansum(comp{j}(k).totFeat([1,2,4,5])); % comp{j}(k).Nspine/comp{j}(k).Nprot;
        comp{j}(k).filoFrac = comp{j}(k).totFeat(4)/nansum(comp{j}(k).totFeat([1,2,4,5])); % comp{j}(k).totFeat(4)/comp{j}(k).Nprot;
        comp{j}(k).ambigFrac = comp{j}(k).totFeat(5)/nansum(comp{j}(k).totFeat([1,2,4,5]));
        % Spine dyamics
        comp{j}(k).SF = SF( comp{j}(k).totFeat(1), comp{j}(k).totFeat(2) );
        comp{j}(k).TO = TO( comp{j}(k).totFeat(1), comp{j}(k).totFeat(2), comp{j}(k).totFeat(3) );
        comp{j}(k).elim = 1 - comp{j}(k).SF;
        comp{j}(k).form = comp{j}(k).totFeat(3)/(comp{j}(k).totFeat(1)+comp{j}(k).totFeat(2));
        % Filopodia dynamics
        comp{j}(k).Felim = NaN; % 
        comp{j}(k).Fform = nansum(comp{j}(k).Nfeat(:,4,2))/nansum(comp{j}(k).Nfeat(:,4,1));
        % Calculate these quantities for each individual neuron
        for n = 1:numel(mouse(j).neurTiles) %mouse(j).Nneur
            comp{j}(k).neurNfeat{n} = Nfeat(mouse(j).neurTilesInd{n},:,:);
            comp{j}(k).neurTotFeat(n,:) = nansum( comp{j}(k).neurNfeat{n}(:,:,1), 1 );
            comp{j}(k).neurTotFeat(n,4) = comp{j}(k).neurTotFeat(n,4) + nansum(comp{j}(k).neurNfeat{n}(:,4,2), 1);
            comp{j}(k).neurNprot(n) = nansum( comp{j}(k).neurTotFeat(n,1:5) );
            comp{j}(k).neurNspine(n) = nansum( comp{j}(k).neurTotFeat(n,1:3) );
            comp{j}(k).neurSpineFrac(n) = nansum(comp{j}(k).neurTotFeat(n,[1,2,5]))/nansum(comp{j}(k).neurTotFeat(n,[1,2,4,5])); % comp{j}(k).neurNspine(n)/comp{j}(k).neurNprot(n);
            comp{j}(k).neurFiloFrac(n) = comp{j}(k).neurTotFeat(n,4)/nansum(comp{j}(k).neurTotFeat(n,[1,2,4,5])); % comp{j}(k).neurTotFeat(n,4)/comp{j}(k).neurNprot(n);
            comp{j}(k).neurAmbigFrac(n) = comp{j}(k).neurTotFeat(n,5)/nansum(comp{j}(k).neurTotFeat(n,[1,2,4,5])); % comp{j}(k).neurTotFeat(n,5)/comp{j}(k).neurNprot(n);
            comp{j}(k).neurSF(n) = SF( comp{j}(k).neurTotFeat(n,1), comp{j}(k).neurTotFeat(n,2) );
            comp{j}(k).neurTO(n) = TO( comp{j}(k).neurTotFeat(n,1), comp{j}(k).neurTotFeat(n,2), comp{j}(k).neurTotFeat(n,3) );
            comp{j}(k).neurElim(n) = 1 - comp{j}(k).neurSF(n);
            comp{j}(k).neurForm(n) = comp{j}(k).neurTotFeat(n,3)/( comp{j}(k).neurTotFeat(n,1)+ comp{j}(k).neurTotFeat(n,2) );
            comp{j}(k).neurFelim(n) = NaN; %
            comp{j}(k).neurFform(n) = nansum(comp{j}(k).Nfeat(mouse(j).neurTilesInd{n},4,2))/nansum(comp{j}(k).Nfeat(mouse(j).neurTilesInd{n},4,1)); 
        end
    end
    % get the number of each feature type for each tile at each day-to-day combination
    daysFeat{j} = nan( mouse(j).Ndays, mouse(j).Ndays, mouse(j).Ntiles, featParam.Ntype ); 
    for k = 1:Ncomp(j)
        xtemp = flip( find( ~cellfun( @isempty, comp{j}(k).scoring(:,1) ) )' );
        [~,~, dInd(2)] = intersect( comp{j}(k).days(2), mouse(j).days );
        [~,~, dInd(1)] = intersect( comp{j}(k).days(1), mouse(j).days );
        for x = xtemp
            for w = 1:featParam.Ntype
                daysFeat{j}( dInd(1), dInd(2), x, w ) = comp{j}(k).Nfeat(x,w,1);
                daysFeat{j}( dInd(2), dInd(1), x, w ) = comp{j}(k).Nfeat(x,w,2);
            end
        end
    end
end
fprintf('\n');
toc;
end

%{
% How would dynamics change if filopodia were considered?
filoElim = 0.8; % assume this fraction of filopodia are eliminated 
comp{j}(k).E2 = (comp{j}(k).totFeat(2)+filoElim*nansum(comp{j}(k).Nfeat(:,4,1)))/(comp{j}(k).totFeat(1)+comp{j}(k).totFeat(2)+nansum(comp{j}(k).Nfeat(:,4,1)));
comp{j}(k).F2 = (comp{j}(k).totFeat(3)+nansum(comp{j}(k).Nfeat(:,4,2)))/(comp{j}(k).totFeat(1)+comp{j}(k).totFeat(2)+nansum(comp{j}(k).Nfeat(:,4,1)));
comp{j}(k).SF2 = 1 - comp{j}(k).E2;
comp{j}(k).TO2 = (comp{j}(k).totFeat(2)+filoElim*nansum(comp{j}(k).Nfeat(:,4,1))+comp{j}(k).totFeat(3)+nansum(comp{j}(k).Nfeat(:,4,2)))/(comp{j}(k).totFeat(2)+filoElim*nansum(comp{j}(k).Nfeat(:,4,1))+comp{j}(k).totFeat(3)+nansum(comp{j}(k).Nfeat(:,4,2))+2*comp{j}(k).totFeat(1));
%}