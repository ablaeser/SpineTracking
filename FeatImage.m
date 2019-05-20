function feat = FeatImage( mouse, comp, feat, featParam, varargin )
% FeatImage generates cropped summary images of feature 
% -----Inputs------
% mouse = metadata about experiments performed on each mouse
% comp = metadata about each comparison analyzed from each mouse
% feat = structure containing information about each feature
% featParam
%   .crop = dimensions to crop around feature
%   .Nneigh = how many frames before/after the feature frame to average over
% -----Outputs------
% feat = cell array of structures containing summary image and various properties about each feature 
%   .cropVec
%   .eyeCrop
%   .reg
%   .meanShift
%   .filt

Nmouse = numel( comp );
tic;
IP = inputParser;
addRequired( IP, 'mouse', @isstruct )
addRequired( IP, 'comp', @iscell )
addRequired( IP, 'feat', @iscell )
addRequired( IP, 'featParam', @isstruct )
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'setj', flip(1:Nmouse), @isnumeric )
addParameter( IP, 'setk', [], @isnumeric )
addParameter( IP, 'setw', [], @isnumeric )
addParameter( IP, 'setx', [], @isnumeric )
addParameter( IP, 'setp', [], @isnumeric )
addParameter( IP, 'setq', [], @isnumeric )
parse( IP, mouse, comp, feat, featParam, varargin{:} );
show = IP.Results.show;
% determine whether to set indices manually or not
setj = IP.Results.setj;
setk = IP.Results.setk; 
setw = IP.Results.setw; 
setx = IP.Results.setx; 
setp = IP.Results.setp; 
setq = IP.Results.setq; 
% setk = []; setw = []; setx = []; setp = []; setq = [];
kToggle = isempty( setk ); % determine whether to set k manually or not
xToggle = isempty( setx );
wToggle = isempty( setw );
pToggle = isempty( setp );
qToggle = isempty( setq );

fprintf('Generating feature summary images...  ');
if show
    close all;
    rawFig = figure('units','normalized', 'OuterPosition',[-0.0036 0.43 0.5 0.55], 'Color','w', 'name', 'raw');
    summFig = figure('units','normalized', 'OuterPosition',[0.5 0.05 0.52 0.93], 'Color','w', 'name','summary');
end
for j = setj
    fprintf('\nj = %d / %d (%s)  ',j, Nmouse, mouse(j).ID);
    try
        if kToggle
            setk = flip( 1:numel( feat{j} ) );
        end
        for k = setk 
            %fprintf('\n    k = %d / %d ',k, numel(feat{j}));
            if xToggle
                setx = find( nansum(comp{j}(k).Nfeat(:,:,1), 2) )'; %flip( find(~cellfun(@isempty, feat{j}{k}(:,1)))' )
            end
            for x = setx  
                %fprintf('\n        x = %d',x) 
                if pToggle
                    setp = [2,1];
                end
                for p = setp %[2,1]
                    %fprintf('\n[j,k,x,p] = [%d, %d, %d, %d]',j,k,x,p)
                    d = find( comp{j}(k).days(p) == mouse(j).days );
                    t = comp{j}(k).timepoints(p);
                    stackPath = sprintf('%s%s.tif',mouse(j).daysDir{d}, mouse(j).stacks{x,d,t});
                    rawStack = loadtiff( stackPath ); % fprintf('\nLoaded %s', stackPath )
                    [~,~, Nframe] = size(rawStack); 
                    if wToggle
                        setw = find( cellfun( @numel, feat{j}{k}(x,:,p) ) ); 
                    end
                    for w = setw
                        if qToggle
                            setq = numel( feat{j}{k}{x,w,p} ):-1:1;
                        end
                        for q = setq 
                            %fprintf('\nj = %d, k = %d, x = %d, p = %d, w = %d, q = %d\n',j,k,x,p,w,q)
                            % Determine which frames to use and make a substack (tempStack)
                            ztemp = feat{j}{k}{x,w,p}(q).eye(3)-featParam.Nneigh:feat{j}{k}{x,w,p}(q).eye(3)+featParam.Nneigh;
                            ztemp(ztemp < 1) = []; ztemp(ztemp > Nframe) = []; Nztemp = numel(ztemp);
                            zref = find(ztemp == feat{j}{k}{x,w,p}(q).eye(3));
                            tempStack = rawStack(:,:,ztemp);
                            % Crop the raw images (tempStack) (optional)
                            clearvars rawCrop; % rawCrop may have different sizes, depending on whether cropped area falls completely on the image or not
                            if all(featParam.crop ~= 0)
                                cropParam = [round(feat{j}{k}{x,w,p}(q).eye(1:2) - featParam.crop/2), featParam.crop];  cropParam(cropParam < 1) = 1;
                                feat{j}{k}{x,w,p}(q).cropVec = cropParam;
                                feat{j}{k}{x,w,p}(q).eyeCrop = [feat{j}{k}{x,w,p}(q).eye(1)-cropParam(1)+1, feat{j}{k}{x,w,p}(q).eye(2)-cropParam(2)+1];
                                %feat{j}{k}{x,w,p}(q).raw = imcrop( rawStack(:,:,feat{j}{k}{x,w,p}(q).eye(3)), cropParam);
                                for z = Nztemp:-1:1
                                    rawCrop(:,:,z) = imcrop( tempStack(:,:,z), cropParam ); 
                                end
                            else
                                rawCrop = tempStack;
                                feat{j}{k}{x,w,p}(q).eyeCrop = feat{j}{k}{x,w,p}(q).eye(1:2);
                                %feat{j}{k}{x,w,p}(q).raw = rawStack(:,:,feat{j}{k}{x,w,p}(q).eye(3));
                            end
                            rawMean = mean( rawCrop, 3); rawMax = max( rawCrop, [], 3);
                            % stackIn = rawCrop;
                            [regCrop, shift, regMean, regMax] = StackReg( rawCrop, zref, 'show', false ); % Average neighboring frames using registration of cropped images
                            feat{j}{k}{x,w,p}(q).reg = shift;
                            feat{j}{k}{x,w,p}(q).meanShift = nanmean( shift(sum(shift,2)>0, 3) );
                            feat{j}{k}{x,w,p}(q).filt  = regMax;
                            % Show the results
                            if show
                                opt = { [0.04,0.02], [0.05,0.03], [0.01,0.01]}; % {gapVector (between subplots) [vert, horz], margin height [bottom, top], margin width [left, right] } 
                                % SUMMARY IMAGES
                                figure( summFig );
                                projPlot(4) = subtightplot(2,2,4, opt{:});
                                imshow( regMax, [] ); 
                                title('Crop-Registered Max Projection'); 
                                projPlot(3) = subtightplot(2,2,2, opt{:});
                                imshow( rawMax, [] ); 
                                title('Naive Max Projection');                           
                                projPlot(2) = subtightplot(2,2,3, opt{:});
                                imshow( regMean, [] ); 
                                title('Crop-Registered Averaging'); 
                                projPlot(1) = subtightplot(2,2,1, opt{:});
                                imshow( rawMean, []); title('Naive Average');
                                hold on; plot( feat{j}{k}{x,w,p}(q).eyeCrop(1), feat{j}{k}{x,w,p}(q).eyeCrop(2), 'r.', 'MarkerSize', 10  ); hold off;
                                linkaxes( projPlot, 'xy' ); 
                                % RAW AND REGISTERED SUBIMAGES
                                figure( rawFig );
                                for z = Nztemp:-1:1
                                    tempInd = sub2ind( [Nztemp, 2], z, 2 );
                                    sp(tempInd) = subtightplot(2,Nztemp,tempInd, opt{:});
                                    imshow( regCrop(:,:,z), [] ); hold on; %plot( cropCoord{w}(q,1), cropCoord{w}(q,2), 'r+', 'MarkerSize', MS  );
                                    title(sprintf('[%2.1f, %2.1f, %2.1f, %2.1f]', shift(z,1), shift(z,2), shift(z,3), shift(z,4)));
                                    tempInd = sub2ind( [Nztemp, 2], z, 1);
                                    sp(tempInd) = subtightplot(2,Nztemp,tempInd, opt{:});
                                    imshow( rawCrop(:,:,z), [] ); hold on; %plot( cropCoord{w}(q,1), cropCoord{w}(q,2), 'r+', 'MarkerSize', MS  );
                                end
                                linkaxes( sp, 'xy' ); 
                                impixelinfo;
                                pause;
                            end
                        end
                    end
                end
            end
        end
        toc
    catch
        fprintf(' - Failed');
    end
end
fprintf('\n');
end