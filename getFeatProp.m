function Feat = getFeatProp( Feat, featParam ) % , minEcc
%getFeatProp estimates various properties of the feature from a 2D summary image
% -----Inputs-------
% Feat = feat structure 
% featRad = 1) radius to search for local maxima 2) radius around that local maxima calculate intensity from 
% minArea = minimum area of shaft region to be counted
% Rmax = maximum radius (in pixels) around centroid of the feature to consider as shaft, for measuring distance from feature to shaft.
% show = whether to display the analysis process or not
% -----Outputs------Xfeat = NaN; Yfeat = NaN; Tfeat = NaN; Rfeat = NaN; Xshaft = NaN; Yshaft = NaN; Tshaft = NaN; Rshaft = NaN; 
% Feat   
%  .eyeInt = intensity in initial search region 
%  .SNR = signal-to-noise ratio of cropped image, using multilevel thresholding to distinguish S and N. 
%  .binary
%  .cart
%  .pol
%  .int = estimated intensity of the feature
%  .shaftInt = estimated intensity of the shaft
%  .norm = Feat.int/Feat.shaftInt
%  .area = # of suprathreshold pixels within [ feature, shaft ]
%  .radius = effective radius of the spine, ie the radius, in um, of a cirle with  area equal to that measured directly
%  .vol = approximate volume (um^3) of the spine, ie (4/3)pi*(Feat.radius)^3
%  .edge = pixels on edge of feature and neighboring shaft
%  .cent = centroid of the second-pass feature mask
%  .shaftCent = centroid of the shaft mask < Rmax
%  .shaftAx = principal axes of the local shaft
%  .d2s = distance, in um, along the second principal axis, from centroid of the feature to the shaft
% Threshold the cropped image to identify neurite pixels (feature and shaft)
thresh = multithresh( Feat.filt, featParam.threshN); % using multilevel thresholding is more sensitive to dimmish signal than graythresh
cropBin = logical(imquantize( Feat.filt, thresh ) - 1); % ignore the subthreshold pixels  
%Feat.binary = cropBin;
% Define first-pass mask 
cropSize = size( Feat.filt );
[rr, cc] = meshgrid( 1:cropSize(2), 1:cropSize(1) );
if Feat.type == 6, featRad = featParam.axonRad; else, featRad = featParam.spineRad; end
eyeMask = sqrt((rr-Feat.eyeCrop(1)).^2+(cc-Feat.eyeCrop(2)).^2) <= featRad(1); 
Feat.eyeInt = mean( Feat.filt(eyeMask) );
eyeBin = eyeMask & cropBin; 
%figure; subplot(1,3,3); imshow( eyeBin ); subplot(1,3,2); imshow( eyeMask ); subplot(1,3,1); imshow( cropBin ); hold on; plot( Feat.eyeCrop(1), Feat.eyeCrop(2), 'w*' ); impixelinfo;
% Use first-pass centroid to analyze the feature
if any( eyeBin(:) )
    [eyeBinY,eyeBinX] = ind2sub( cropSize, find(eyeBin) ); %eyeBinY = cropSize(1) - eyeBinY;
    Feat.xyCrop = round( wmean( [eyeBinX,eyeBinY], double(repmat(Feat.filt(eyeBin),[1,2])), 1 ) ); 
    offset = Feat.xyCrop - Feat.eyeCrop;
    Feat.xy = Feat.eye(1:2) + offset;
    % Define second-pass feature mask around Feat.xyCrop
    featCirc = sqrt( (rr-Feat.xyCrop(1)).^2+(cc-Feat.xyCrop(2)).^2 ) <= featRad(2); % imshow( featCirc )
    featBin = featCirc & cropBin;  % imshow( featBin )
    Feat.area = sum( featBin(:) ); % & cropBin(:)
    [Yfeat,Xfeat] = ind2sub( size(Feat.filt), find(featBin) ); 
else
    featBin = false(cropSize);
    Feat.edge{1} = [];
    Yfeat = NaN; Xfeat = NaN;
end
% If a usable feature is found, calculate various properties
if Feat.area > 0 
    Feat.int = mean( Feat.filt( featBin ) );
    Feat.sat(1) = numel( find(Feat.filt(featBin) >= 4094) )/Feat.area;
    Feat.edge{1} = bwboundaries( featBin, 4, 'noholes' );
    Feat.edge{1} = cellfun( @flip, Feat.edge{1}, repmat({2},numel(Feat.edge{1}),1), 'UniformOutput', false); % bwboundaries output has [y,x] columns
    Feat.cent = round( wmean([Xfeat,Yfeat],double(repmat( Feat.filt(featBin), 1,2 )),1 ) );
    Feat.perim = size( vertcat( Feat.edge{1}{:} ),1 );
    Feat.radius = featParam.conv(1)*sqrt(Feat.area/pi);
    Feat.vol = (4*pi/3)*(Feat.radius)^3;
    %[Tfeat,Rfeat] = cart2pol( Xfeat-Feat.cent(1), cropSize(1) - Yfeat - (cropSize(1)-Feat.cent(2)) );
else
    Feat.cent = nan(1,2);
    %Tfeat = NaN; Rfeat = NaN;
end
%{
% Identify putative shaft region(s)
shaftBin = cropBin & ~featBin; % imshow( shaftBin ); % exclude feature pixels from consideration in upcoming shaft detection
CC = bwconncomp(shaftBin,4); 
RP = regionprops( CC, 'Area' ); % ,'Orientation' , 'Eccentricity'
CC = CC.PixelIdxList; 
CC = CC( [RP.Area] > featParam.minArea(1) ); % only keep sufficiently large CCs 
shaftBin(:) = false; shaftBin( vertcat(CC{:})' ) = true; % imshow( cropBin ) % remake shaftBin after removing all small (< minArea) blobs 
Feat.edge{2} = bwboundaries( shaftBin, 4, 'noholes' ); Nshaft = numel( Feat.edge{2} );

% Get cartesian and polar coordinates of each signal pixel (for Sholl-ish analysis)   
Feat.binary = featBin | shaftBin;
[Ybin,Xbin] = ind2sub( cropSize, find( Feat.binary ) ); %Y = cropSize(1) - Y;
[Tbin,Rbin] = cart2pol( Xbin-Feat.cent(1), cropSize(1) - Ybin - (cropSize(1)-Feat.cent(2)) );
Feat.cart = single([Xbin,Ybin]);
Feat.pol = single([Tbin,Rbin]);
% Get properties of shaft region
if Nshaft > 0
    Feat.edge{2} = cellfun( @flip, Feat.edge{2}, repmat({2},Nshaft,1), 'UniformOutput', false); % bwboundaries output has [y,x] columns
    Feat.shaftArea = sum( shaftBin(:) ); 
    % Analyze the shaft pixels within Rmax of Feat.cent
    [Yshaft,Xshaft] = ind2sub( size(Feat.filt), find(shaftBin) );
    [Tshaft,Rshaft] = cart2pol( Xshaft-Feat.cent(1), cropSize(1) - Yshaft-(cropSize(1)-Feat.cent(2)) );
    badShaft = Rshaft > featParam.Rmax; %round(featParam.maxD2S/featParam.conv(1));
    Xshaft( badShaft ) = []; Yshaft( badShaft ) = []; Tshaft( badShaft ) = []; Rshaft( badShaft ) = [];
    goodShaftInd = sub2ind( size(Feat.filt), Yshaft, Xshaft );
    if ~isempty( Xshaft )
        % Get the distances from each feature pixel to each good shaft pixel, and find how many shaft pixels are 'rubbing' against the feature region
        if Feat.area > 0 
            allDist = pdist2( [Xfeat, Yfeat], [Xshaft,Yshaft] ); % imshow( featBin ); hold on; pause; imshow( shaftBin );
            minDist = min(allDist); %min( allDist, [], 2 );
            Feat.Nrub = numel( find( minDist < 1.1 ) );
            Feat.rubRatio = Feat.Nrub/Feat.perim;
        end
        %shaftIm = zeros( cropSize ); shaftIm( Yshaft, Xshaft ) = Feat.filt( Yshaft, Xshaft ); % imshow( shaftIm, [] ); imshow( Feat.filt, [] )
        shaftIm = zeros( cropSize ); shaftIm(goodShaftInd) = Feat.filt(goodShaftInd);  % imshow( shaftIm, [] ); 
        [Feat.shaftCent, ~, imMom, Feat.shaftAx ] = ImMoments( shaftIm );
        Feat.shaftInt = mean( shaftIm(shaftIm > 0 )); % mean( Feat.filt( shaftBin(:) ) );
        Feat.norm = double(Feat.int)/double(Feat.shaftInt); %NOTE: THIS VALUE IS NOW OVERWRITTEN BY getDendFeat
        Feat.sat(2) = numel( find(shaftIm >= 4094) )/Feat.shaftArea;
        f2s = Feat.shaftCent-Feat.cent; % the vector that goes between the feature and the shaft centroid
        Feat.angle = acosd( sum(f2s.*Feat.shaftAx(1:2,1)')/norm(f2s) ); % angle (degrees) between shaft's principal axis and f2s 
        Feat.d2s = abs(featParam.conv(1)*sum(f2s.*Feat.shaftAx(1:2,2)')); % how far (um) does the feature protrude out from the shaft (dot product of f2s and second principal axis)
        signal = double( mean(Feat.filt( Feat.binary )) - mean(Feat.filt( ~Feat.binary )) ); % signal = mean foreground value - mean background value
        noise = std( double( Feat.filt(~Feat.binary) ) ); % noise = std dev of background pixel values
        Feat.SNR = signal/noise;
    end
else
    Feat.edge{2} = [];
    Feat.angle = NaN; Feat.d2s = NaN; Feat.sat(2) = NaN;
    Tshaft = NaN; Rshaft = NaN;
end
if show
    opt = {[0.08,0.001], [0.08,0.04], 0.01};  % {gapVector (between subplots) [vert, horz], margin height [bottom, top], margin width [left, right] } 
    LW = 1;
    arrowLength = 10;
    close all;
    figure('units','normalized','Position',[0,0,1,1],'Color','w'); %[0.1 0.1 0.8 0.8]
    % Filtered image with spine and shaft outlined
    subtightplot(2,3,[2,5],opt{:}); 
    imshow( Feat.filt, [], 'InitialMagnification', 'fit' ); 
    hold on;
    axis off; impixelinfo; 
    plot( Feat.eyeCrop(1), Feat.eyeCrop(2), 'mx' );  % 'mx'
    %eyeEdge = bwboundaries( eyeMask, 4, 'noholes' );
    %plot( eyeEdge{1}(:,2), eyeEdge{1}(:,1), 'm', 'LineWidth', LW  ); % bwboundaries returns cell of edge pixels [row, column]
    if ~any(isnan(Feat.xyCrop))
        for z = 1:numel( Feat.edge{1} )
            plot( Feat.edge{1}{z}(:,1), Feat.edge{1}{z}(:,2), 'b', 'LineWidth', LW);
        end
        plot( Feat.xyCrop(1), Feat.xyCrop(2), 'bd' )
        plot( Feat.cent(1), Feat.cent(2), 'yx' )
    end
    for z = 1:Nshaft
        plot( Feat.edge{2}{z}(:,1), Feat.edge{2}{z}(:,2), 'r', 'LineWidth', LW);
    end
    if ~isempty( Feat.shaftCent )
        plot( Feat.shaftCent(1), Feat.shaftCent(2), 'rx' );
        line([Feat.shaftCent(1),Feat.cent(1)],[Feat.shaftCent(2),Feat.cent(2)],'Color','k','LineStyle','--', 'LineWidth',1.5*LW)
        line([Feat.shaftCent(1), Feat.shaftCent(1)+arrowLength*Feat.shaftAx(1,1)],[Feat.shaftCent(2), Feat.shaftCent(2)+arrowLength*Feat.shaftAx(2,1)],'Color','r','LineWidth',1.5*LW)
        line([Feat.shaftCent(1), Feat.shaftCent(1)+arrowLength*Feat.shaftAx(1,2)],[Feat.shaftCent(2), Feat.shaftCent(2)+arrowLength*Feat.shaftAx(2,2)],'Color','g','LineWidth',1.5*LW)
    end
    title( sprintf('area = %d, angle = %2.1f, d2s = %2.1f um, rub = %d, sat = [%1.2f, %1.2f]', Feat.area, Feat.angle, Feat.d2s, Feat.Nrub, Feat.sat(1), Feat.sat(2) ) ); % mom = [%2.0f, %2.0f] , imMom(1,1), imMom(2,2) , vol = %2.1f um^3
    % Intensity histograms of full image, feature and shaft
    subtightplot(2,3,4,opt{:});%subplot(1,3,3);
    binWidth = 93; histEdge = 0:binWidth:4095;
    cropHist = histc( Feat.filt(:), histEdge );
    featHist = histc( Feat.filt(featBin(:)), histEdge );
    shaftHist = histc( Feat.filt(shaftBin(:)), histEdge );
    plot( histEdge, cropHist, 'k' ); hold on;
    plot( histEdge, featHist, 'b' ); 
    plot( histEdge, shaftHist, 'r' );
    xlabel('Pixel Value'); ylabel('Frequency');
    title( sprintf('SNR: %2.1f, intensity: %2.0f / %2.0f = %2.2f',Feat.SNR, Feat.int, Feat.shaftInt, Feat.norm ) );
    set(gca,'yscale','log'); axis square;
    legend('All pix','Feature','Shaft');
    % Binarized image in polar coords, showing feature, shaft and extraneous pixels
    subtightplot(2,3,1,opt{:});
    polar(Feat.pol(:,1),Feat.pol(:,2),'k.'); hold on;
    polar(Tfeat,Rfeat,'b.'); 
    polar(Tshaft,Rshaft,'r.');
    if numel( Feat.ind ) == 6
        title(sprintf('%s [j,k,x,w,p,q] = [%d, %d, %d, %d, %d, %d]', Feat.source, Feat.ind(1), Feat.ind(2), Feat.ind(3), Feat.ind(4), Feat.ind(5), Feat.ind(6) ));
    elseif numel( Feat.ind ) == 5
        title(sprintf('%s [j,x,d,t,q] = [%d, %d, %d, %d, %d]', Feat.source, Feat.ind(1), Feat.ind(2), Feat.ind(3), Feat.ind(4), Feat.ind(5) ));
    end
    axis square;
    subtightplot(2,3,3,opt{:});
    polTemp = Feat.pol( Feat.pol(:,2) < featParam.Rmax,:);
    polMed = median( polTemp );
    polMean = mean( polTemp );
    rose( polTemp(:,1) ); hold on;
    polar( polMed(1), polMed(2), 'k.' )
    polar( polMean(1), polMean(2), 'b.' )
    subtightplot(2,3,6,opt{:});
    Nr = hist( polTemp(:,2), 0:1:featParam.Rmax+1 ); xlim([0,featParam.Rmax+1]);
    plot( 0:1:featParam.Rmax+1, Nr, 'k' ); hold on;
    plot( 1:featParam.Rmax+1, diff(Nr), 'b'); 
    plot([featRad(2),featRad(2)],[0,max(Nr)],'r--');
    set( gca, 'Xtick', 0:2:featParam.Rmax+1 );
    xlim([0,20]);
    xlabel('Pixels from center'); ylabel('Suprathreshold Pixels'); axis square;
    %spaceplots;
    pause; 
end
%}
end

% Show the process (optional)
%{
figure; 
subplot(1,6,1); imshow(Feat.filt,[]); title('Feat.filt');
subplot(1,6,2); imshow( cropBin, [] ); title('cropBin');
subplot(1,6,3); imshow( eyeBin ); title('eyeBin');
subplot(1,6,4); imshow( featBin); title('featBin');
subplot(1,6,5); imshow( shaftBin ); title('shaftBin');
subplot(1,6,6); imshow( Feat.binary ); title('Feat.binary');
impixelinfo;
%}