function offset = getFeatMatch( feat1, feat2, pixDist, stepDist)
%Match features from one time point to those of a second time point
% -----Inputs------
% feat1, feat2 = coordinates of matched features at 2 time points
% pixDist = microns per pixel
% stepDist = microns per Z-step
% -----Outputs------
%

% feat1 = feat(3,:,1); feat2 = feat(3,:,2);
Nfeat1 = cellfun(@numel, feat1)/3; Nfeat2 = cellfun(@numel, feat2)/3;
% Determine which types of feats to use for alignment
matchedTypes = []; %[2,3,5]; % these types should always be matched
for w = [6,5,3,2] %[1,6] % these types might be matched
    if Nfeat1(w) > 0 && Nfeat2(w) == Nfeat1(w)
        matchedTypes = [matchedTypes, w]; %#ok<*AGROW>
    end
end
if isempty( matchedTypes )
    error('No matched features')
end
matchedTypes = sort( matchedTypes, 'descend' ); % start with boutons or ambiguous, since these are probably the most accurate registrations 
% Estimate z-offset (# of frames)
allFeat1 = vertcat( feat1{matchedTypes} ); Z1 = allFeat1(:,3); XY1 = allFeat1(:,1:2); 
allFeat2 = vertcat( feat2{matchedTypes} ); Z2 = allFeat2(:,3); XY2 = allFeat2(:,1:2); 
Zoff = median(Z2) - median(Z1); 
% Perform independent alignments for each feature type 
XYoff = round( mean( XY2 ) - mean( XY1 ) ); % initialize a running estimate of XY offset 
for w = matchedTypes
    offsetAll{w} = zeros(0,5);
    if Nfeat(w,1) == Nfeat(w,2)
        distMat = zeros( Nfeat(w,1) ); minDist = zeros( Nfeat(w,1), 3 );
        for n = 1:Nfeat(w,1)
            for m = 1:Nfeat(w,2)
                XYsep = feat{w,2}(m,1:2) - feat{w,1}(n,1:2) - XYoff;
                Zsep = feat{w,2}(m,3) - feat{w,1}(n,3) - Zoff;
                distMat(n,m) = norm( [pixDist, stepDist].*[XYsep, Zsep]);
            end
            [tempDist,ind] = min( distMat(n,:) ); % index of closest match
            minDist(n,1) = n; minDist(n,2) = ind; minDist(n,3) = tempDist; % minDist = [ day 1 feature #, matched day 2 feature #, separation ]
        end

        % Match the closest available features, correcting for offsets in XY and Z
        matched = zeros( Nfeat(w,1), 1 ); offset = zeros( Nfeat(w,1), 5);
        [~,sortInd] = sort( minDist(:,3), 'ascend' ); sortInd = sortInd';
        for n = sortInd
            [matchDist,matched(n)] = min( distMat(n,:) ); % match dist is the _offset-corrected_ distance between matched features
            offset(n,:) = [feat{w,1}(n,3), feat{w,2}(matched(n),:) - feat{w,1}(n,:), matchDist]; % [z1, (x2-x1), (y2-y1), (z2-z1), matchDist]
            distMat(:,matched(n)) = NaN; % prevent multiple features from matching to one feature
        end
        if any( offset(:,5) > 10 )
            warning('%s-d%d-%d_%s has %d %s nearest match(es) very far apart', mouse, day(1), day(2), tile{x}, numel( find(offset(:,5) > 10) ), types{w});
        end
        %disp( matched ); %disp( distMat );  disp( minDist );       
        offsetAll{w} = [offsetAll{w}; offset]; % [XYoffset, Zoffset];
        XYoff = round( median( offsetAll{w}(:,2:3) ) );


    else
        error('%d-%d_%s has unequal number of %s features', day(1), day(2), tile{x}, types{w} );
    end
end
tempOff = vertcat( offsetAll{:} );
offset{x} = round( median( tempOff(:,2:4), 1 ) );


%%
for w = matchedTypes
    comp.offsetAll{x,w} = zeros(0,5);
    if comp.Nfeat(x,w,1) == comp.Nfeat(x,w,2)
        distMat = zeros( comp.Nfeat(x,w,1) ); minDist = zeros( comp.Nfeat(x,w,1), 3 );
        for n = 1:comp.Nfeat(x,w,1)
            for m = 1:comp.Nfeat(x,w,2)
                XYsep = comp.feat{x,w,2}(m,1:2) - comp.feat{x,w,1}(n,1:2) - XYoff;
                Zsep = comp.feat{x,w,2}(m,3) - comp.feat{x,w,1}(n,3) - Zoff;
                distMat(n,m) = norm( [pixDist, stepDist].*[XYsep, Zsep]);
            end
            [tempDist,ind] = min( distMat(n,:) ); % index of closest match
            minDist(n,1) = n; minDist(n,2) = ind; minDist(n,3) = tempDist; % minDist = [ day 1 feature #, matched day 2 feature #, separation ]
        end

        % Match the closest available features, correcting for offsets in XY and Z
        matched = zeros( comp.Nfeat(x,w,1), 1 ); offset = zeros( comp.Nfeat(x,w,1), 5);
        [~,sortInd] = sort( minDist(:,3), 'ascend' ); sortInd = sortInd';
        for n = sortInd
            [matchDist,matched(n)] = min( distMat(n,:) ); % match dist is the _offset-corrected_ distance between matched features
            offset(n,:) = [comp.feat{x,w,1}(n,3), comp.feat{x,w,2}(matched(n),:) - comp.feat{x,w,1}(n,:), matchDist]; % [z1, (x2-x1), (y2-y1), (z2-z1), matchDist]
            distMat(:,matched(n)) = NaN; % prevent multiple features from matching to one feature
        end
        if any( offset(:,5) > 10 )
            warning('%s-d%d-%d_%s has %d %s nearest match(es) very far apart', comp.mouse, comp.day(1), comp.day(2), comp.tile{x}, numel( find(offset(:,5) > 10) ), types{w});
        end
        %disp( matched ); %disp( distMat );  disp( minDist );       
        comp.offsetAll{x,w} = [comp.offsetAll{x,w}; offset]; % [XYoffset, Zoffset];
        XYoff = round( median( comp.offsetAll{x,w}(:,2:3) ) );


    else
        error('%d-%d_%s has unequal number of %s features', comp.day(1), comp.day(2), comp.tile{x}, types{w} );
    end
end
tempOff = vertcat( comp.offsetAll{x,:} );
comp.offset{x} = round( median( tempOff(:,2:4), 1 ) );
% Predict locations of nonmatched features
nonmatchedTypes = 1:size(comp.Nfeat,2); 
nonmatchedTypes( [4,matchedTypes] ) = []; % filopodia are never matched
for w = nonmatchedTypes
    if comp.Nfeat(x,w,1) > 0
        comp.feat{x,w,2} = comp.feat{x,w,1} + repmat( comp.offset{x}, [comp.Nfeat(x,w,1),1] );
        % Make sure matched comp.feat{2} lives within the actual stack
        comp.feat{x,w,2}( comp.feat{x,w,2}(:,1) < 1, 1) = 1; 
        comp.feat{x,w,2}( comp.feat{x,w,2}(:,1) > comp.resolution(1), 1) = comp.resolution(1);
        comp.feat{x,w,2}( comp.feat{x,w,2}(:,2) < 1, 2) = 1;
        comp.feat{x,w,2}( comp.feat{x,w,2}(:,2) > comp.resolution(2), 2) = comp.resolution(2);
        comp.feat{x,w,2}( comp.feat{x,w,2}(:,3) < 1, 3) = 1;
        comp.feat{x,w,2}( comp.feat{x,w,2}(:,3) > comp.Nframe(x,2), 3) = comp.Nframe(x,2);
    end
end

end