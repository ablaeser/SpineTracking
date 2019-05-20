clear; clc;
% Determine which mouse, tile and days to analyze, and generate paths to save files and figures to
setMouse = '2505'; 
setTile = 'C';
mainDir = 'D:\Andy\Drexel\'; 
slm = '\';
% Feat Parameters
res = [512,512];  blankIm = zeros( res, 'uint16' ); 
dist20x = [585.6, 585.6]; % dimensions of FOV for 20x objective with 1x digital zoom
featZoom = 8;
featDist = dist20x./res/featZoom; % um/pixel
featStep = 0.5; 
pixDist = dist20x./res/featZoom; % um/pixel
spineRad = [3,3];
axonRad = [3,8];
featParam = struct('type',[], 'Ntype',NaN, 'conv',[pixDist, featStep], 'crop',127*[1,1], 'Nneigh',0, 'threshN',3, 'minArea',[200,100], 'spineRad',[3,3], 'Rmax',round(5/pixDist(1)), 'axonRad',[3,8], 'maxD2S',5, 'minLength',30, 'minSpine', 70, 'edge',4);  
featParam.type = {'Stable','Elim','Formed','Filo','Ambig','Bouton','Tip'}; featParam.Ntype = numel(featParam.type);
screen = get(0,'ScreenSize');  
figDir = 'D:\Andy\Drexel\Figures'; %  [figMain,datestr(today,'yyyy-mm-dd'),slm]; mkdir(figDir);
[mouse, T] = getMouse( mainDir, {setMouse}, 'meta',true ); 
Nmouse = numel(setMouse); J = find( strcmpi( setMouse, {mouse.ID} ) ); 
X = find( strcmpi( setTile, [mouse.tiles] ) );
setDays = mouse.days; %
[~,D,~] = intersect( mouse(J).days, setDays ); D = D'; % indices of the days included/excluded from analysis
%TP = repmat({1},1,mouse(J).Ndays); 
saveDir = [mouse(J).dir,'Longitudinal Analysis',slm]; mkdir( saveDir );
saveRoot = [saveDir, setMouse,'_',setTile,'-'];
tempMatPath = [saveRoot,'longTEMP.mat']; % [saveRoot,'axonTEMP.mat'];
matPath = [saveRoot,'long.mat'];
pdfRoot = [mouse(J).figDir, setMouse,'_',setTile,'-'];
% Load all relevant tif stacks and stackInfo
tic;
Nstack = sum(mouse(J).Nstacks(X,D)); %sum(mouse(J).Nstacks(X,D).*(setDays./setDays));
s = 0; 
Istack = cell(Nstack,1); Nframe = zeros(Nstack,1); S = nan(1,Nstack); DT = nan(Nstack,2);
for d = find( ~cellfun(@isempty, mouse.stacks(X,:,1)) )  % D   
    for t = 1:mouse(J).Nstacks(X,d)% TP{d}
        s = s + 1;
        %S(s) = sub2ind( size(mouse(J).stacks), X, d ); fprintf('s = %d, S = %d',s,S(s));  % linear indices of the stacks used for analysis
        S(s) = sub2ind( size(mouse(J).stacks), X, d, t );
        DT(s,:) = [d,t];
        tifPath = sprintf('%s%s.tif',mouse(J).daysDir{d}, mouse(J).stacks{S(s)});
        Istack{s} = loadtiff( tifPath ); fprintf('\ns = %d. Loaded %s.    ', s, [mouse(J).daysDir{d},mouse(J).stacks{S(s)}] );
        Nframe(s) = size(Istack{s},3);
        toc
    end
end
sTrack = 1:Nstack;
Ntrack = numel( Nframe(Nframe > 0 )); Strack = S(sTrack);
neighOpt = round( (screen(3)*Ntrack/screen(4) - 1)/2 ); % optimal # of neighboring frames to fill the entire screen while keeping the images as large as possible
%stackInfo = StackFacts( mouse(J), 'load',true); % load metadata about the stacks, including alignment 
%StackAlign( mouse, stackInfo, 'show', true, 'setx', X); %
%%
% Use the earliest comparison of setTile, with the most spines, for seed spines.
[mouse, comp, Ncomp] = getComp( mouse ); 
[comp, feat, daysFeat] = FeatRead( mouse, comp, featParam );
totSpine = nansum(daysFeat{J}(:,:,X,1:3),4); 
dSeed(1) = find( setDays(1) == mouse(J).days ); 
[~,dSeed(2)] = max( totSpine(dSeed(1),:) ); % second day with most spines identified in a dSeed(1)-based comparison
K = mouse(J).daysComp(dSeed(1),dSeed(2),X); 
P = find( comp{J}(K).days == mouse(J).days(dSeed(1)) ); % comparison for dSeed(1)->dSeed(2)? comp{J}(K).tiles{X,:,P};
pTemp = [1,2];
tSeed = [1,1]; %
feat = FeatImage( mouse, comp, feat, featParam, 'show',false, 'setj', J, 'setk', K ); % , 'setj', 10:12 , 'setj', 13
feat = FeatProps( feat, featParam, 'setj', J, 'setk', K ); % , 'show',true, 'setw', 6
[comp, feat] = FeatMatch( mouse, comp, feat, 'show',false );
% Select the marking with the greatest number of features to use as the seed day/timepoint
Sseed(2) = sub2ind( size(mouse(J).stacks), X, dSeed(2), tSeed(2) );
Sseed(1) = sub2ind( size(mouse(J).stacks), X, dSeed(1), tSeed(1) );
[~,~,sSeed] = intersect( Sseed, S ); sSeed = sSeed';

seed = []; 
seed = [seed; vertcat( feat{J}{K}{X,6,1} )];
seed = [seed; vertcat( feat{J}{K}{X,1,1} )]; 
seed = [seed; vertcat( feat{J}{K}{X,2,1} )];
seed = [seed; vertcat( feat{J}{K}{X,3,2} )];
seed = [seed; vertcat( feat{J}{K}{X,4,:} )]; 
Nseed = size(seed,1); 

seedID = [seed.ID];
seedXYZ = vertcat( seed.eye );
seedUm = repmat( featParam.conv, Nseed, 1 ).*seedXYZ;
seedSep = squareform( pdist( seedUm ) ); seedSep( seedSep == 0 ) = NaN;
seedType = reshape( [seed.ind], 6, [] )'; seedType = seedType(:,4);

qBouton = find( seedType == 6 )'; Nbouton = numel( qBouton );
qSpine = find( seedType > 0 & seedType <4 )'; Nspine = numel( qSpine );
qFilo = find( seedType == 4 )'; Nfilo = numel( qFilo );

% Initialize variables for use in scoring (only create these variables the first time scoring a given tile - otherwise use the loaded version)
zAlign = nan(Nseed, Nstack); 
long = nan( Nseed, 4, Nstack ); Nlong = size(long,1); 
% Fill pre-scored features
qStable = find(seedType == 1)'; qElim = find(seedType == 2)'; qForm = find(seedType == 3)';
% Boutons
long(qBouton, 1:3, sSeed(1) ) = seedXYZ(qBouton,:);  
long( qBouton, 4, sSeed(1) ) = 1;
% Stable spines
%tempStable = vertcat( feat{J}{K}{X,1,2} );
%long( qStable, 1:3, sSeed(2) ) = vertcat( tempStable.eye );  %long( qStable, 1:4, sSeed(2) ) = NaN;
long( qStable, 1:3, sSeed(1) ) = seedXYZ(qStable,:);  long( qStable, 4, sSeed(1) ) = 1;
long( qStable, :, sSeed(2) ) = NaN;
% Eliminated spines
tempElim = vertcat( feat{J}{K}{X,2,2} );
if ~isempty( tempElim )
    long( qElim, 1:3, sSeed(2) ) = vertcat( tempElim.eye );  long( qElim, 4, sSeed(2) ) = 3;
    long( qElim, 1:3, sSeed(1) ) = seedXYZ(qElim,:);  long( qElim, 4, sSeed(1) ) = 1;
end
% Formed spines
tempForm = vertcat( feat{J}{K}{X,3,1} );
if ~isempty( tempForm )
    long( qForm, 1:3, sSeed(1) ) = vertcat( tempForm.eye );  long( qForm, 4, sSeed(1) ) = 3;
    long( qForm, 1:3, sSeed(2) ) = seedXYZ(qForm,:);  long( qForm, 4, sSeed(2) ) = 1;
end

longMat = nan( Nlong, Nstack );
Nnan = NaN; qNan = NaN; qFull = NaN; Nfull = NaN; qInc = NaN; Ninc = NaN; qStable = NaN; qUnstable = NaN; Nstable = NaN; Nunstable = NaN;
%long(qBouton,1:3,sSeed(1))
q = 1; % which feature to start scoring at
saveVars = {'long','Nlong','longMat','q','setTile','X','K','setDays','pTemp','D','DT','S','totSpine','dSeed','tSeed','sSeed','pTemp','sTrack','Strack','Nstack','Nscore','Ntrack','Nframe','zAlign',...
    'seed','seedID','seedXYZ','seedUm','seedSep','seedType','Nseed','qNan','Nnan','qFull','Nfull','qInc','Ninc','qStable','Nstable','qUnstable','qSpine','qFilo','qBouton','Nunstable'}; % 

%% Load previous analysis
load( matPath ); fprintf('\nLoaded %s\n', matPath); % load final version
%load( tempMatPath ); fprintf('\nLoaded %s\n', tempMatPath); % load temporary version
saveVars = {'long','Nlong','longMat','q','setTile','X','K','setDays','pTemp','D','DT','S','totSpine','dSeed','tSeed','sSeed','pTemp','sTrack','Strack','Nstack','Nscore','Ntrack','Nframe','zAlign',...
    'seed','seedID','seedXYZ','seedUm','seedSep','seedType','Nseed','qNan','Nnan','qFull','Nfull','qInc','Ninc','qStable','Nstable','qUnstable','qSpine','qFilo','qBouton','Nunstable'}; % 

%% Add new stack(s)

sNew = 6:9; fprintf('\nAdding stack %d\n', sNew ); 
Nstack = numel(Istack);
sTrack = [sTrack, sNew]; Ntrack = numel(sTrack); 
long(:,:,sNew) = NaN; 
for s = sNew
    Nframe(s) = size(Istack{s},3); 
    DT(s,:) = [s,1];
end
setDays(sNew) = mouse.days(sNew);


%% Find and reset bad entries
[qBad, sBad] = ind2sub( size(squeeze( long(:,4,:))), find( squeeze( long(:,4,:) ) > 3 ) );
Nbad = numel(qBad);
fprintf('\nFound %d bad entries\n   qBad   sBad\n', Nbad); disp( [qBad,sBad] )
for n = 1:Nbad
    long( qBad(n), 1:4, sBad(n) ) = NaN;
end

%% Reset specific entries
qReset = [15];% [203,206,236,239,240];% find( long(:,3,3) > 170 );
sReset = [10]; %[1:Ntrack]; % ; %2; % 
long( qReset,:,sReset) = NaN; fprintf('\nReset %d\n', qReset);


%% Display features across all time points and mark
% Input parameters of scoring
sScore = [1,3,4,6,7]; %1:Nstack; % [1,4,10,11]; % [1,5:9];% 1:5;% [1,5,6]; %  [1,4:8];  %1:4; % %   which stacks to actively score here
Nscore = numel(sScore); 
neighOpt = 5; % round( (screen(3)*Nscore/screen(4) - 1)/2 ); % 5 % 3; % optimal # of neighboring frames to fill the entire screen while keeping the images as large as possible
neighFrame = 2; %neighOpt; %neighOpt;  %5;%
jumpFrame = 1; % show images this many frames apart
maxSep = 20;
Nsurr = 9; 
Nsave = 2; % interval between saves
%
longMat = reshape( long(:,4,:), Nlong, Nstack );
longMat = longMat(:,sTrack);
nanMat = longMat(:,1:Ntrack); nanMat( ~isnan(nanMat) ) = 0; nanMat( isnan(nanMat) ) = 1;
Nnan = sum(nanMat,2);
qNan = find( Nnan > 1 & Nnan < Ntrack )';   %& Nnan <= 5 find( ~isnan( AstroMat(:,end) ) )'; % find( isnan( sum( AstroMat(:,dSet), 2 ) ) )';  & Nnan < Nstack
qFull = find(Nnan==0)'; Nfull = numel(qFull);  % fully analyzed  % =Nstack-Ntrack
qEmpty = find( Nnan == Ntrack )'; Nempty = numel(qEmpty); % not analyzed at all 
qInc = setdiff( 1:Nseed, [qFull,qEmpty] ); Ninc = numel(qInc); % partially analyzed
qStable = find( sum(longMat,2) == Ntrack )'; Nstable = numel(qStable);
qUnstable = setdiff( qFull, qStable ); Nunstable = numel(qUnstable);
%{
sHas = 1:5;% sScore; %[1,7];  %
sHasnot = 6:9; %[2:6]; %
qHas = find( sum(nanMat(:,sHas),2) == 0 )'; 
qHasnot = find( sum(nanMat(:,sHasnot),2) > 0 )';
qHHN = intersect( qHas, qHasnot );
%}
%{
umMat = squareform( pdist( repmat(featParam.conv,Nlong,1).*long(:,1:3,1) ) ); 
umMat( umMat==0 ) = NaN; 
[minDist, ~] = min( umMat, [], 2 );
[sortDist,qDist] = sort( minDist, 'descend');
qDist( isnan(sortDist) ) = [];
%}
%qStable = find(seedType == 1)'; qElim = find(seedType == 2)'; qForm = find(seedType == 3)';
% Commence scoring
close all; 
longFig = figure('units','normalized','OuterPosition',[0,0,1.1,1],'color',0.6*[1,1,1],'MenuBar','none');  % [1-0.1*Nstack,0,0.1*Nstack,1]'k'
mymenu = uimenu('Parent',longFig,'Label','Hot Keys');
uimenu('Parent',mymenu,'Label','Zoom','Accelerator','z','Callback',@(src,evt)zoom(longFig,'on'));
uimenu('Parent',mymenu,'Label','Unzoom','Accelerator','x','Callback',@(src,evt)zoom(longFig,'out'));
uimenu('Parent',mymenu,'Label','Pointer','Accelerator','c','Callback',@(src,evt)zoom(longFig,'off'));
totFrame = 2*neighFrame+1; % show totFrame images for each day
opt = {[0.001,0.001], [0.001,0.02], 0};  % {gapVector (between subplots) [vert, horz], margin height [bottom, top], margin width [left, right] } 
clickColor = {'b','m','g'};
colorMat = distinguishable_colors(Nbouton+Nspine+Nfilo);  % Nbouton+Nspine+Nfilo  Nsurr
spSize = [totFrame, Nscore]; Nsp = totFrame*Nscore; spPos = nan( Nsp, 4 ); 
sp(Nsp) = subtightplot( Nscore, totFrame, Nsp, opt{:} ); 
qTemp = intersect(qUnstable, qSpine);
 
for q = 67 %qUnstable %intersect(qInc,qSpine) % qFilo %qInc %[qInc(qInc>65)] %qBouton%35:Nlong % 1:Nlong %qInc%qFilo %:Nlong%qSpine %  qSpine %24:Nlong% qSpine % qSpine%  % 
    fprintf('\nq = %d / %d', q, Nseed ); 
    [tempSep, qNear] = sort( seedSep(q,:), 'ascend' ); % sort features from closest to furthest
    qNear = qNear( tempSep < maxSep );
    p = seed(q).ind(5); fprintf(', sSeed = %d', sSeed(p) );
    for s = sScore
        qAlign = qNear( long(qNear,4,s) == 1 ); % only keep features that have been observed on the day in question days
        try
            if long(q,3,s) > 0 % feature has already been marked in MATLAB
                zAlign(q,s) = long(q,3,s);
            elseif s == sSeed(p) %s == sSeed % feature has already been marked in ImageJ
                zAlign(q,s) = seedXYZ(q,3);
            elseif numel( qAlign ) > 6 % enough NEARBY features have already been marked in MATLAB
                zNear = long(qAlign,3,s);
                zAlign(q,s) = round( seedXYZ(q,3) - nanmedian(long(qAlign,3,sSeed(p))-zNear) );        
            elseif isnan(long(q,3,s)) && numel( find(~isnan(long(:,3,sSeed(p))-long(:,3,s) )) ) > 4 % enough features have already been marked in MATLAB
                zAlign(q,s) = round( seedXYZ(q,3) - nanmedian(long(:,3,sSeed(p))-long(:,3,s)) ); 
            else 
                zAlign(q,s) = round( Nframe(s)*seedXYZ(q,3)/Nframe(sSeed(p)) );
            end
            fprintf('\n    s = %d: zAlign = %d', s, zAlign(q,s) ); 
        catch
            warning('Failed to align: s = %d', s)
        end
        for f = -neighFrame:neighFrame
            zShift = zAlign(q,s) + jumpFrame*f;
            if zShift <= Nframe(s) && zShift > 0
                % CONTRAST ENHANCEMENT OPTIONS
                tempIm = Istack{s}(:,:,zShift); % always on
                %tempIm = medfilt2( tempIm, 2*[1,1], 'symmetric' ); % median filter to reduce salt and pepper noise
                %tempIm = mean( Istack{s}(:,:,zshift), 2  ); % median filter to reduce salt and pepper noise
                tempIm = imadjust( tempIm, [],[], 0.7 ); % adjust image for better contrast
                %tempIm = imadjust( tempIm, stretchlim(tempIm, 0.5) ); % adjust image for better contrast
                %tempIm = adapthisteq( tempIm, 'Range','original', 'Distribution','rayleigh', 'NumTiles',[16,16], 'ClipLimit', 0.001 ); % adaptively adjust image for better contrast
            else 
                fprintf('\nUsing blankIm: q=%d, s= %d f=%d, %d',q, s, f );
                tempIm = blankIm; 
            end
            spInd = sub2ind( spSize, neighFrame+f+1, find(s==sScore) );
            sp(spInd) = subtightplot( Nscore, totFrame, spInd, opt{:} );
            cla;
            spPos(spInd,:) = get( sp(spInd), 'Position' ); % use position vector of each axis to identify which frames are picked
            imshow( tempIm, [] ); hold on; axis off;
            %title( sprintf('z = %d', zShift ), 'FontSize',8, 'HorizontalAlignment','left' );
            if zShift == long(q,3,s)
                plot( long(q,1,s), long(q,2,s), '+', 'Color', 'b', 'LineWidth',2 ); % clickColor{long(q,4,s)}  , 'MarkerSize', 7,
                %plot( long(q,1,s), long(q,2,s), 'o', 'Color', 'b', 'LineWidth',2 );  %  , 'MarkerSize', 7
            end
            %pause;
            % LABEL SURROUNDING CELLS
            c = 0;
            for Q = qNear(1:min(Nsurr,numel(qNear)))
                %plot( long(Q,1,s), long(Q,2,s), 'x', 'Color', 'y', 'MarkerSize', 3 );
                if ismember(long(Q,4,s), [1,2])
                    c = c + 1;
                    %text( long(Q,1,s), long(Q,2,s), sprintf('%d',Q), 'Color', 0*[1,1,1], 'HorizontalAlignment','center', 'FontSize',6,'FontWeight','bold'); % clickColor{long(Q,4,s)}
                    %text( long(Q,1,s), long(Q,2,s), sprintf('%d',Q), 'Color', 1*[1,1,1], 'HorizontalAlignment','center', 'FontSize',8.5,'FontWeight','bold');
                    %plot( long(Q,1,s), long(Q,2,s), '*', 'color', colorMat(Q,:), 'LineWidth',1.5 ); % , 'MarkerSize', colorMat(c,:) , 'MarkerSize', 7
                end
            end
        end
    end
    linkaxes( sp(:), 'xy' ); 
    %set(gca,'Xlim',[302,423], 'Ylim',[137,258])
    %ipi = impixelinfo; set(ipi, 'units','normalized','Position',[0.5, 0.982, 0.1, 0.1] )
    zoom(longFig,'out');
    if find(sSeed(p)==sScore)
        subtightplot( Nscore, totFrame, sub2ind(spSize, (totFrame+1)/2, find(sSeed(p)==sScore) ), opt{:} ); 
        %plot( seedXYZ(q,1), seedXYZ(q,2), 'rx',  'MarkerSize', 5); %plot( seedXYZ(q,1), seedXYZ(q,2), 'ro','MarkerSize', 34);
    end
    subplot( sp(1) );  title( sprintf('q = %d / %d, click 0 / %d', q, Nseed, Nscore), 'FontSize',8, 'HorizontalAlignment','left' );
    toc
    pause; % pause after each click to allow for zooming 
    for c = 1:Nscore
        try
            [xClick, yClick, clickBut, ~, clickPos] = ginputax( 1 ); 
            fprintf('\n  c = %d: [X,Y,Button,Pos] = [ %05.1f, %05.1f, %1d, %02.4f]', c,xClick,yClick,clickBut, sum( clickPos, 2 )); %
            if xClick > 0 && xClick < res(1) && yClick > 0 && yClick < res(2)
                [fClick, sClick] = ind2sub( spSize, find( sum( clickPos, 2 ) == sum( spPos, 2 ) ) ); %determine which subplot was clicked
                Zclick = zAlign(q,sScore(sClick)) + jumpFrame*(fClick-1-neighFrame);
                long(q,:,sScore(sClick)) = [ round(xClick), round(yClick), Zclick, clickBut ]; 
                subtightplot( Nscore, totFrame, sub2ind(spSize, fClick, sClick ), opt{:} );
                plot( long(q,1,sScore(sClick)), long(q,2,sScore(sClick)), '.', 'Color', clickColor{clickBut}, 'MarkerSize', 10 ); 
            end
        catch
            fprintf('\n  c = %d: Missing click',c); 
        end
        if c < Nscore
            subplot( sp(1) )
            title( sprintf('q = %d / %d, click c = %d / %d', q, Nseed, c, Nscore ), 'FontSize',8, 'HorizontalAlignment','left' );
            pause;
        end
    end
    if rem(q,Nsave) == 0,  save( tempMatPath, saveVars{:} ); fprintf('\nSaved %s\n', tempMatPath);  end
    if any( long(q,4,:) > 3 ), error('\nq = %d, s = %d: scoring value greater than 3!\n', q, find( ~ismember(squeeze(long(q,4,long(q,1,:) > 0)),1:3))'  ); end
end
close all;
tic;
Nlong = size(long,1); % disp(long);
longMat = reshape( long(:,4,:), Nlong, Nstack );
longMat = longMat(:,sTrack);
nanMat = longMat; nanMat( ~isnan(nanMat) ) = 0; nanMat( isnan(nanMat) ) = 1;
Nnan = sum(nanMat,2);
qNan = find( Nnan > 1 & Nnan < 3 )';  % find( ~isnan( AstroMat(:,end) ) )'; % find( isnan( sum( AstroMat(:,dSet), 2 ) ) )';  & Nnan < Nstack
qFull = find(Nnan==0)'; Nfull = numel(qFull);  % fully analyzed  % =Nstack-Ntrack
qEmpty = find( Nnan == Ntrack )'; Nempty = numel(qEmpty); % not analyzed at all 
qInc = setdiff( 1:Nseed, [qFull,qEmpty] ); Ninc = numel(qInc); % partially analyzed
qStable = find( sum(longMat,2) == Ntrack )'; Nstable = numel(qStable);
qUnstable = setdiff( qFull, qStable ); Nunstable = numel(qUnstable);
save( matPath, saveVars{:} ); fprintf('\nSaved %s\n', matPath); 

%print(gcf,'-dtiff', 'D:\Andy\Drexel\2505\Figures\TrackingExampleNoNeigh2.tif')

%% Check the results
q = 1;
stretchTol = 0.001;
clip = 0.01;
MS = 40;
opt = {[0.003,0.001], [0.03,0.02], 0}; % {gapVector (between subplots) [vert, horz], margin height [bottom, top], margin width [left, right] } 
clickColor = {'r','m','g'};
%close all;
M = 3; %round(screen(3)/screen(4)); %1; %
N = 4; %ceil(Ntrack/M); % 3; %
figure('units','normalized','position',[0,0,1,1],'color','w'); %  [0.18,0.04,0.64,0.87] [0.25,0,0.5,1]
for q = 9 % qUnstable % intersect( find( sum(longMat,2) ~= Nstack  )', qBouton ) %find( sum(longMat,2) ~= Nstack & ~isnan(sum(longMat,2)) )'%qUnstable%qStable %1:Nlong %[qBouton,qSpine,qFilo] % [qSpine,qFilo] %
    for s = find( squeeze( nansum( long(:,3,:) ) )' ) %flip(sTrack)
        p = seed(q).ind(5);
        sp(s) = subtightplot( M, N, find(s==sTrack), opt{:} );
        if ~isnan(long(q,3,s)) && long(q,3,s) > 0
            tempIm = Istack{s}(:,:,long(q,3,s));
            %tempIm = imadjust( tempIm, stretchlim(tempIm, stretchTol) ); % adjust image for better contrast
            %tempIm = adapthisteq(tempIm, 'Range','original', 'Distribution','rayleigh', 'NumTiles',[8,8], 'ClipLimit', 0.02);
            imshow( tempIm, [] ); hold on; 
            plot( long(q,1,s), long(q,2,s), 'o','MarkerSize',MS, 'Color', clickColor{long(q,4,s)} ); 
        else
            imshow( blankIm ); hold on;
        end
        axis off;
        if s == sSeed(p)
            plot( seedXYZ(q,1), seedXYZ(q,2), 'ro','MarkerSize',MS);
            %title( sprintf('z = %d', long(q,s,3) );
        end
        if s == 1
            title( sprintf('q = %d / %d, ', q, size(long,1)) );
        end
    end
    linkaxes( sp(:),'xy'); impixelinfo; 
    pause; % pause to allow for zooming 
    clf;
end
close all;

%% Show alignments
opt = {[0.09,0.07], [0.07,0.05], [0.05,0.01]};  % {gapVector (between subplots) [vert, horz], margin height [bottom, top], margin width [left, right] } 
figure('units','normalized','position',[0,0,1,1],'color','w'); % 
for s = sTrack
    % Use the matched points to do a linear regression alignment
    qGood = find( ~isnan(long(:,4,sSeed(1)) ) & ~isnan(long(:,4,s) ) )';
    qBad = setdiff( 1:Nseed, qGood );
    xMatch = long(qGood,3,sSeed(1)); yMatch = long(qGood,3,s);
    p = polyfit( xMatch, yMatch, 1);
    xTemp = 1:Nframe(sSeed(1)); %1:length(stackInfo.Zalign{DT(sSeed(1),1),DT(s,1),X,DT(sSeed(1),2),DT(s,2)});
    Zmatch{s} = polyval(p, xTemp);
    % Make the plot
    subtightplot(2,ceil(Ntrack/2),find(s==sTrack), opt{:});
    plot( long(qGood,3,sSeed(1)), long(qGood,3,s), 'k.', 'MarkerSize', 9 ); hold on;
    plot( long(qBad,3,sSeed(1)), long(qBad,3,s), 'ko', 'MarkerSize', 4 );
    plot( xTemp, Zmatch{s}, 'k' );
    %plot( stackInfo.Zalign{DT(sSeed(1),1),DT(s,1),X,DT(sSeed(1),2),DT(s,2)}, 'r' ); 
    xlim([0,Nframe(sSeed(1))+1]); ylim([0,max(Nframe)+1]); axis square;
    xlabel( sprintf('z (%s)', mouse(J).stacks{X,DT(sSeed(1),1), DT(sSeed(1),2)}), 'Interpreter','none' ); ylabel(sprintf('Z (%s)', mouse(J).stacks{X,DT(s,1), DT(s,2)}), 'Interpreter','none');% title( sprintf( 
    title( sprintf('s = %d: N = %d features', s, numel(qGood) ) );
    %legend('Matched','Unmatched','Regression','Location','NorthWest');
end
%% Summarize the results (max proj)
clickColor = {'m','y','w'};
close all;
longFullSummary = figure('units','normalized','position',[0,0,1,1],'color','w','PaperOrientation','landscape','PaperPositionMode','manual','PaperUnits','normalized','PaperPosition',[0,0,1,1]); %  [0.18,0.04,0.64,0.87] [0.25,0,0.5,1]
M = 3; %round(screen(3)/screen(4)); %1; %
N = 4; %ceil(Ntrack/M); % 3; %
opt = {[0.02,0.001], [0.02,0.02], [0.001,0.001]}; % {gapVector (between subplots) [vert, horz], margin height [bottom, top], margin width [left, right] } 
for s = flip(sTrack) % 1:4 % 
    try
        d = DT(s,1); t = DT(s,2); 
        zMin = min(long(:,3,s)); zMax = max(long(:,3,s));
        sp(s) = subtightplot( M, N, find(s==sTrack), opt{:} );
        zTemp = long(:,3,s); zTemp(isnan(zTemp)) = []; zTemp = unique(zTemp);
        tempIm = max( Istack{s}(:,:,zTemp), [], 3 );  % max( Istack{s}(:,:,zMin:zMax), [], 3 );  %  round(Zmatch{s})
        %tempIm = adapthisteq(tempIm);  %  adjust image for better contrast , 'Range','original', 'Distribution','rayleigh', 'NumTiles',[8,8], 'ClipLimit', 0.02 
        imshow( tempIm, [] ); hold on; axis off;
        hold on; 
        title( sprintf('%s-d%d %s', setMouse, mouse(J).days(d), setTile), 'FontSize', 10, 'VerticalAlignment','bottom' );
        for q = Nlong:-1:1
            if ~isnan(long(q,4,s))
                %plot( long(q,1,s), long(q,2,s), '.', 'Color', clickColor{long(q,4,s)},'MarkerSize',6); axis off;
                %text(long(q,1,s) - textOffset(1), long(q,2,s) + textOffset(2), num2str(q), 'Color', clickColor{long(q,4,s)}, 'FontSize', FS, 'HorizontalAlignment', 'right' ); %
                text( long(q,1,s), long(q,2,s), sprintf('%d',q), 'Color', clickColor{long(q,4,s)} , 'HorizontalAlignment','center', 'FontSize',16, 'FontWeight','bold'); % clickColor{long(q,4,s)}
                %text( long(q,1,s), long(q,2,s), sprintf('%d',q), 'Color', 'k' , 'HorizontalAlignment','center', 'FontSize',7,'FontWeight','bold'); % clickColor{long(q,4,s)}
            end
        end
    catch
        fprintf('\ns = %d', s);
    end
end
linkaxes( sp(:),'xy'); 
pdfPath = [pdfRoot,'longFullSummary.pdf']; % [pdfRoot,'axonFullSummary.pdf'];
export_fig( pdfPath, '-pdf', '-painters','-r600','-q101', longFullSummary ); fprintf('\nSaved %s\n', pdfPath);
impixelinfo;


%% Generate a multi-page PDF summarizing the trajectories of each feature
Nneigh = 0;
maxSep = 10;
stretchTol = 0.001;
perPage = 4; %round( (screen(4)*Ntrack/screen(3) - 1) ); % optimal # of neighboring sections to use to fill the entire screen while keeping the images as large as possible % # of features to show on each page of the PDF
neighColor = distinguishable_colors(10);
clickColor = {'g','r','w'};
opt =  {[0.015,0.004], [0.005,0.02], [0.001,0.001]}; % {gapVector (between subplots) [vert, horz], margin height [bottom, top], margin width [left, right] } 
Npage = ceil(Nlong/perPage); %ceil(Nbouton/perPage); % 
pdfPath = [pdfRoot,'featTrajectories.pdf']; % [figRoot,'axonTrajectories.pdf'];
if exist( pdfPath, 'file') == 2, delete( pdfPath ); end % erase the preexisting file if needed (otherwise it will just append to that file)
close all;
featTrajectories = figure('units','normalized','position',[0,0,1,1],'color','w','PaperOrientation','landscape','PaperPositionMode','manual','PaperUnits','normalized','PaperPosition',[0,0,1,1]);
for p = 1:Npage %1%:
    qPage = (1:perPage)+(p-1)*perPage;     qPage = qPage( qPage <= Nlong );
    for s = sTrack
        for q = 1:length(qPage)
            d = DT(s,1); t = DT(s,2);
            subtightplot( perPage, Ntrack, sub2ind([Ntrack,perPage], find(s==sTrack), q ), opt{:} );
            if ~isnan(long(qPage(q),4,s)) && long(qPage(q),4,s) > 0
                cropParam = [round(long(qPage(q),1:2,s) - featParam.crop/2), featParam.crop];  cropParam(cropParam < 1) = 1;
                xyCrop = long(qPage(q),1:2,s)-cropParam(1:2)+[1,1];%[long(qPage(q),1,s)-cropParam(1)+1, long(qPage(q),2,s)-cropParam(2)+1];
                [tempSep, qNear] = sort( seedSep(qPage(q),:), 'ascend' ); % sort features from closest to furthest
                qNear = qNear( tempSep < maxSep ); Nnear = numel(qNear);
                nearCrop = long(qNear,1:2,s) - repmat(cropParam(1:2)+[1,1], Nnear, 1);
                zProj = long(qPage(q),3,s)-Nneigh:long(qPage(q),3,s)+Nneigh; zProj(zProj<1 | zProj > size(Istack{s},3)) = [];
                [~,~,~,tempIm] = StackReg( Istack{s}(:,:,zProj), -1 ); 
                %tempIm = Istack{s}(:,:,long(qPage(q),3,s));
                tempIm = imcrop( tempIm, cropParam ); 
                %tempIm = imadjust( tempIm, stretchlim(tempIm, stretchTol) ); % adjust image for better contrast
                imshow( tempIm, [] ); hold on;
                plot( xyCrop(1), xyCrop(2), '.', 'MarkerSize',15, 'Color', 'k' ); 
                plot( xyCrop(1), xyCrop(2), 'o', 'MarkerSize',3, 'Color', clickColor{long(qPage(q),4,s)} ); % ,
                for n = 1:min([Nnear,10])
                    plot( nearCrop(n,1), nearCrop(n,2), '+', 'MarkerSize',4, 'Color', neighColor(n,:) )
                end
            else
                imshow( blankIm ); %imcrop( blankIm, cropParam)
            end
            axis off;
            if ismember(long(qPage(q),4,s), [1,2]), featStatus = 'present';  elseif long(qPage(q),4,s) == 3, featStatus = 'absent';  else, featStatus = 'not scored'; end
            if ismember(qPage(q), qBouton), featType = 'bouton';  elseif long(qPage(q),4,s) == 1, featType = 'spine';  elseif long(qPage(q),4,s) == 2, featType = 'filopodia';  else featType = 'protrusion'; end %#ok<*SEPEX>
            title( sprintf('%s_%s-d%d: q = %d %s %s',mouse(J).ID, mouse(J).tilename{X,d,t}, setDays(d), qPage(q), featType, featStatus ), 'FontSize', 8, 'FontName','Arial', 'Color','k', 'HorizontalAlignment', 'center', 'interpreter','none' );
        end
    end
    %pause
    export_fig( pdfPath, '-pdf', '-painters','-r600','-q101', '-append', featTrajectories ); fprintf( 'Saved %s (page %d / %d) \n',pdfPath, p, Npage );  
    clf;
end
close( featTrajectories );

%% Summarize positions of each feature in z, across all timepoints
Zsummary = figure('units','normalized','position',[0,0,1,1],'color','w','PaperOrientation','landscape','PaperPositionMode','manual','PaperUnits','normalized','PaperPosition',[0,0,1,1]);
plot( -1*squeeze( long(qInc,3,:) )', 'r' ); hold on; 
plot( -1*squeeze( long(qUnstable,3,:) )', 'b' ); 
plot( -1*squeeze( long(qStable,3,:) )', 'k', 'LineWidth', 1.5 );
%set(gca, 'Xtick', 1:3, 'XtickLabel', {'d06','d15','d24'} );
ylabel('Z (um below first frame)');
title( longName{L} );

%pdfPath = [figDir,'astroZsummary.pdf']; % [pdfRoot,'axonFullSummary.pdf'];
%export_fig( pdfPath, '-pdf', '-painters','-r600','-q101', Zsummary ); fprintf('\nSaved %s\n', pdfPath);


