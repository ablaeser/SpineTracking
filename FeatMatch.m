function [comp, feat] = FeatMatch( mouse, comp, feat, varargin )
%FeatRead reads all the scoring files listed in comp structure, collects the coordinates of each feature in comp structures, and creates feat structs for each feature
% -----Inputs-------
% mouse
% comp
% feat
% -----Outputs------
% comp
%  .match = list of matched feature indices from the other side of the comparison
% feat
%  .match = [index of best-matched feature from the other side of the comparison, distance between the matched points, ambiguity (ratio of distance between nearest and next-nearest matches), translation invariant error ]
IP = inputParser;
addRequired( IP, 'mouse', @isstruct )
addRequired( IP, 'comp', @iscell )
addRequired( IP, 'feat', @iscell )
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'setj', 1:numel(mouse), @isnumeric )
addParameter( IP, 'setk', [], @isnumeric )
addParameter( IP, 'setx', [], @isnumeric )
addParameter( IP, 'setw', [], @isnumeric )
parse( IP, mouse, comp, feat, varargin{:} );
show = IP.Results.show;
setj = IP.Results.setj;
setk = IP.Results.setk; kToggle = isempty( setk ); % determine whether to set k manually or not
setx = IP.Results.setx; xToggle = isempty( setx );
setw = IP.Results.setw; wToggle = isempty( setw );
% If show is on, set up the figure
if show
    LW = 2;
    figure('units','normalized','position',[0,0,1,1],'color','w'); 
    sp(2) = subplot(1,2,2); sp(1) = subplot(1,2,1); linkaxes( sp, 'xy');
end
% Match features 
wtry = [1,5,6,8]; % check whether these feature types are matched
for j = setj 
    if kToggle
        setk = flip( 1:numel( feat{j} ) );
    end
    for k = setk 
        comp{j}(k).match = cell( size(comp{j}(k).scoring,1), size(comp{j}(k).Nfeat,2) );
        dNfeat = comp{j}(k).Nfeat(:,:,1) - comp{j}(k).Nfeat(:,:,2);
        if xToggle
            setx = find( ~isnan(dNfeat(:,1)') ); %
        end
        for x = setx   
            wmatch = find( dNfeat(x,:) == 0 & comp{j}(k).Nfeat(x,:,1) > 0 ); % which types are matched across days
            wland = intersect( wmatch, wtry ); % which matched types should be used for calculating offset (landmarks)
            Nland = sum( comp{j}(k).Nfeat(x,wland,1), 2 ); % how many total landmark feature are there?
            if Nland > 3
                % Match each feature to its nearest correspondent in the other set, without replacement.
                xyz2 = vertcat( comp{j}(k).feat{x,wland,2} ); xyz1 = vertcat( comp{j}(k).feat{x,wland,1} );
                matchOff = mean( xyz2(:,1:2) ) - mean( xyz1(:,1:2) ); % estimate (x,y) offset using matched features
                %matchOff = mean( vertcat( comp{j}(k).um{x,wland,2} ) ) - mean( vertcat( comp{j}(k).um{x,wland,1} ) ); % estimate (x,y,z) offset using matched features
                if wToggle
                    setw = wmatch; % intersect( find( comp{j}(k).Nfeat(x,:,1) ), [1,2,3,5,6,7] ); % [1,2,3,5,6,7];
                end
                for w = setw %wland %intersect( find( comp{j}(k).Nfeat(x,:,1) ), wtry )
                    %fprintf('\n[j,k,x,w] = [ %d, %d, %d, %d ]', j, k, x, w )
                    distMat = pdist2( comp{j}(k).feat{x,w,1}(:,1:2)+repmat(matchOff, comp{j}(k).Nfeat(x,w,1), 1), comp{j}(k).feat{x,w,2}(:,1:2) );
                    %distMat = pdist2( comp{j}(k).um{x,w,1}+repmat(matchOff, comp{j}(k).Nfeat(x,w,1), 1), comp{j}(k).um{x,w,2} ); 
                    E = nan(comp{j}(k).Nfeat(x,w,1),2); closeDist = nan(comp{j}(k).Nfeat(x,w,1),2);
                    if comp{j}(k).Nfeat(x,w,1) > 1
                        tempDistMat = distMat;
                        comp{j}(k).match{x,w} = zeros(comp{j}(k).Nfeat(x,w,1),4);
                        ambig = nan(size(distMat,1),1);
                        for q = 1:size(distMat,1)
                            distSort = sort( distMat(q,:), 'ascend' );
                            ambig(q) = distSort(1)/distSort(2);
                        end
                        [~,qsort] = sort(ambig, 'ascend'); % sort( minDist, 'ascend' );
                        for q = qsort' % starting from the least ambiguous pair, use distance and image error to match each feature
                            %fprintf('\n[j,k,x,w,q] = [%d, %d, %d, %d, %d]',j,k,x,w,q)
                            [tempDist, tempQ] = sort( tempDistMat(q,:) );
                            closeDist(q,:) = tempDist(1:2);
                            im{2} = feat{j}{k}{x,w,2}(tempQ(1)).filt; xy(2,:) = flip( feat{j}{k}{x,w,2}(tempQ(1)).xyCrop );
                            im{1} = feat{j}{k}{x,w,1}(q).filt; xy(1,:) = flip( feat{j}{k}{x,w,1}(q).xyCrop );
                            E(q,1) = ImageError( im, xy, false ); % 
                            im{2} = feat{j}{k}{x,w,2}(tempQ(2)).filt; xy(2,:) = flip( feat{j}{k}{x,w,2}(tempQ(2)).xyCrop );
                            E(q,2) = ImageError( im, xy, false ); %
                            [~,minInd] = min( E(q,:) );
                            %{
                            if ~any( isnan(E(q,:)) )
                                [~,minInd] = min( closeDist(q,:) ); % use distance when error cannot be calculated (unequally sized images)
                            else
                                [~,minInd] = min( E(q,:) );
                            end
                            %}
                            tempDistMat(:,tempQ(minInd)) = NaN; % prevent multiple spines from matching to the same spine
                            comp{j}(k).match{x,w}(q,:) = [ tempQ(minInd), tempDist(minInd), closeDist(q,1)/closeDist(q,2), E(q,minInd) ];
                            feat{j}{k}{x,w,1}(q).match = comp{j}(k).match{x,w}(q,:);
                            feat{j}{k}{x,w,1}(q).matchXYZ = comp{j}(k).feat{x,w,2}(tempQ(minInd),:);
                            feat{j}{k}{x,w,2}(tempQ(minInd)).match = [q, comp{j}(k).match{x,w}(q,2:4)];
                            feat{j}{k}{x,w,2}(tempQ(minInd)).matchXYZ = comp{j}(k).feat{x,w,1}(q,:);
                        end
                    else
                        %fprintf('\n[j,k,x,w,q] = [%d, %d, %d, %d, 1]',j,k,x,w)
                        closeDist = [distMat, NaN];
                        im{2} = feat{j}{k}{x,w,2}(1).filt; xy(2,:) = flip( feat{j}{k}{x,w,2}(1).xyCrop );
                        im{1} = feat{j}{k}{x,w,1}(1).filt; xy(1,:) = flip( feat{j}{k}{x,w,1}(1).xyCrop );
                        E(1) = ImageError( im, xy, false );
                        comp{j}(k).match{x,w}(1,:) = [ 1, distMat, 0, E(1) ];
                        feat{j}{k}{x,w,1}(1).match = comp{j}(k).match{x,w}(1,:);
                        feat{j}{k}{x,w,1}(1).matchXYZ = comp{j}(k).feat{x,w,2}(1,:);
                        feat{j}{k}{x,w,2}(1).match = [1, comp{j}(k).match{x,w}(1,2:4)];
                        feat{j}{k}{x,w,2}(1).matchXYZ = comp{j}(k).feat{x,w,1}(1,:);
                    end    
                    % Show the matched features side by side (optional)
                    if show
                        for q = 1:comp{j}(k).Nfeat(x,w,1) % 
                            fprintf('\n[j,k,x,w] = [ %d, %d, %d, %d ], q = %d', j, k, x, w, q )
                            qmatch = feat{j}{k}{x,w,1}(q).match(1);
                            subplot(sp(1));
                            imshow( imadjust( feat{j}{k}{x,w,1}(q).filt ), [] ); hold on; axis off;
                            try 
                                plot( feat{j}{k}{x,w,1}(q).eyeCrop(1), feat{j}{k}{x,w,1}(q).eyeCrop(2), 'gx' );
                                plot( feat{j}{k}{x,w,1}(q).xyCrop(1), feat{j}{k}{x,w,1}(q).xyCrop(2), 'bd' );  % 'gx'
                                for z = 1:numel( feat{j}{k}{x,w,1}(q).edge{1} )
                                    plot( feat{j}{k}{x,w,1}(q).edge{1}{z}(:,1), feat{j}{k}{x,w,1}(q).edge{1}{z}(:,2), 'b', 'LineWidth', LW);
                                end
                                for z = 1:numel( feat{j}{k}{x,w,1}(q).edge{2} )
                                    plot( feat{j}{k}{x,w,1}(q).edge{2}{z}(:,1), feat{j}{k}{x,w,1}(q).edge{2}{z}(:,2), 'r', 'LineWidth', LW);
                                end
                            catch
                                fprintf('\nFailed to display xyCrop or edge for p=1 feature')
                            end
                            hold off;
                            title( sprintf('%s [j,k,x] = [%d, %d, %d], w = %d, q = %d / %d', feat{j}{k}{x,w,1}(q).source, ....
                                feat{j}{k}{x,w,1}(q).ind(1), feat{j}{k}{x,w,1}(q).ind(2), feat{j}{k}{x,w,1}(q).ind(3), feat{j}{k}{x,w,1}(q).ind(4), feat{j}{k}{x,w,1}(q).ind(6), comp{j}(k).Nfeat(x,w,1) ))
                            subplot(sp(2)); 
                            imshow( imadjust( feat{j}{k}{x,w,2}(qmatch).filt ), [] ); hold on; axis off;
                            try
                                plot( feat{j}{k}{x,w,2}(qmatch).eyeCrop(1), feat{j}{k}{x,w,2}(qmatch).eyeCrop(2), 'gx' );
                                plot( feat{j}{k}{x,w,2}(qmatch).xyCrop(1), feat{j}{k}{x,w,2}(qmatch).xyCrop(2), 'bd' );  % 'gx'
                                for z = 1:numel( feat{j}{k}{x,w,2}(qmatch).edge{1} )
                                    plot( feat{j}{k}{x,w,2}(qmatch).edge{1}{z}(:,1), feat{j}{k}{x,w,2}(qmatch).edge{1}{z}(:,2), 'b', 'LineWidth', LW);
                                end
                                for z = 1:numel( feat{j}{k}{x,w,2}(qmatch).edge{2} )
                                    plot( feat{j}{k}{x,w,2}(qmatch).edge{2}{z}(:,1), feat{j}{k}{x,w,2}(qmatch).edge{2}{z}(:,2), 'r', 'LineWidth', LW);
                                end
                            catch
                                fprintf('\nFailed to display xyCrop or edge for p=2 feature')
                            end
                            title( sprintf('qmatch = %d, sep = [%2.2f, %2.2f], ambig = %2.2f, E = [%2.2f, %2.2f]', feat{j}{k}{x,w,1}(q).match(1), closeDist(q,1), closeDist(q,2), feat{j}{k}{x,w,1}(q).match(3), E(q,1), E(q,2) ) )
                            hold off;
                            impixelinfo;
                            pause;
                        end
                    end
                end
            else
                fprintf('\nSkipped %s-d%d-%d_%s - only %d landmark features.', comp{j}(k).mouse, comp{j}(k).days(1), comp{j}(k).days(2), mouse(j).tiles{x}, Nland )
            end
        end
    end
end
fprintf('\n');
end