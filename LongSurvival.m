function [ longMeta, transMat, Ntrans, surv, totSurv ] = LongSurvival( mouse, longMeta, longMat, varargin )
% FeatSurvival estimates the survival and non-survival (transition) of spines/boutons, and filopodia
% -----Inputs------
% mouse 
% dFeat = days over which longitudinal scoring of features was conducted
% Nfeat = # of features considered
% featMat = 
% -----Outputs------
% transMat = [transition type (1 through 9), observed longevity of the feature PRIOR to transition]
% Ntrans = # of features that exhibited each transition type
% surv = # of features, of each type, that survived to each timepoint
% dUse
% totSurv = fraction of features still present at each timepoint
Nmouse = numel(mouse);
IP = inputParser;
addRequired( IP, 'mouse', @isstruct )
addRequired( IP, 'longMeta', @isstruct )
addRequired( IP, 'longMat', @iscell )
%addOptional( IP, 'dTot', {}, @iscell );
addParameter( IP, 'fit', '', @ischar ) %addParameter( IP, 'fit', true, @islogical )
addParameter( IP, 'setj', flip(1:Nmouse), @isnumeric )
addParameter( IP, 'show', false, @islogical )
parse( IP, mouse, longMeta, longMat, varargin{:} );
setj = IP.Results.setj;
show = IP.Results.show;
fitType = IP.Results.fit;
if isempty(fitType), fitToggle = false; else, fitToggle = true; end
if show
    FS = 16;
    %load( 'colorblind_colormap.mat' );
    mouseColor = distinguishable_colors(Nmouse);  %colorblind;  % 
    TL = [0.005,0];
    textPad = 0.08;
    opt = {[0.05,0.02], [0.05,0.05], [0.06,0.04]};  % { gap [vert, horz], margin [bottom, top], margin [left, right] } 
end
% Determine which type of model to use to fit the survival data, if any
if fitToggle
    % modelStart = {[0.9,-20], [0.5,-1,-0.1], [1,-1]}; % initial points for model parameters (a,b,c) to search around
    switch fitType
        case '1'
            modelSurv = @(a,b,x)(a*exp(b*x) + (1-a)); % 1 time-constant kinetic model: a = impermanent fraction, b = time constant (1/days)
            modelStart = {[0.1,-0.05], [0.5,-1], [1,-1]}; % initial points for model parameters (a,b) to start search around
        case '2'
            modelSurv = @(a,b,c,x)( a*exp(b*x) + (1-a)*exp(c*x) ); % 2 time-constant kinetic model: a = short-lived fraction, b = short time constant (1/days), c = long time constant (1/days)  
            modelStart = {[0.1,-0.01,-0.005], [0.5,-1,-0.1], [1,-1,-0.005]}; % initial points for model parameters (a,b,c) to start search around
    end
end


fprintf('\nLongSurvival\n');
transMat = cell(1,Nmouse); Ntrans = cell(1,Nmouse);  surv = cell(1,Nmouse); totSurv = cell(1,Nmouse);
for j = setj %1:Nmouse %find( cellfun( @nansum, Nin ) > 0 ) %
    daysTemp = longMeta(j).useDays-longMeta(j).useDays(1);
    totSurv{j}{3} = zeros(1, longMeta(j).Nuse); totSurv{j}{2} = zeros(1, longMeta(j).Nuse); totSurv{j}{1} = zeros(1, longMeta(j).Nuse);
    for x = longMeta(j).x
        transMat{j}{x} = nan( longMeta(j).Nlong(x), 2 );
        for w = 1:3 % 1 = bouton, 2 = spine, 3 = filo
            surv{j}{x,w} = zeros( 1, longMeta(j).Ndays(x) );
            for q = longMeta(j).qType{x,w} %longMeta(j).qFull{x}
                %fprintf('\n q = %d:  ', q); 
                dTrans = find( longMat{j}{x}(q,:) ~= longMat{j}{x}(q,1), 1 );
                if isempty( dTrans )
                    dTrans = longMeta(j).Ntrack(x); %dFeat{j}{x}(end);
                    transMat{j}{x}(q,2) = longMeta(j).days{x}(dTrans)-longMeta(j).days{x}(1); %mouse.days(dTrans)-mouse.days(dFeat{j}{x}(1));
                else
                    transMat{j}{x}(q,2) = longMeta(j).days{x}(dTrans-1)-longMeta(j).days{x}(1);
                end
                transStr = sprintf('%d%d', longMat{j}{x}(q,[1,dTrans]) ); %fprintf('\n[j,x,w,q] = [%d, %d, %d, %d] %s :', j,x,w,q, transStr);
                switch transStr
                    case '11' % stable spine/bouton
                        transMat{j}{x}(q,1) = 1;
                        surv{j}{x,w} = surv{j}{x,w} + ones(1,longMeta(j).Ntrack(x));
                    case '12' % spine->filo
                        transMat{j}{x}(q,1) = 2;
                        surv{j}{x,w}(1:dTrans-1) = surv{j}{x,w}(1:dTrans-1) + 1;
                    case '13' % elim spine/bouton
                        transMat{j}{x}(q,1) = 3;
                        surv{j}{x,w}(1:dTrans-1) = surv{j}{x,w}(1:dTrans-1) + 1;
                    case '21' % filo->spine
                        transMat{j}{x}(q,1) = 4;
                        surv{j}{x,w}(1:dTrans-1) = surv{j}{x,w}(1:dTrans-1) + 1;
                    case '22' % stable filo
                        transMat{j}{x}(q,1) = 5;
                        surv{j}{x,w} = surv{j}{x,w} + ones(1,longMeta(j).Ntrack(x));
                    case '23' % elim filo
                        transMat{j}{x}(q,1) = 6;
                        surv{j}{x,w}(1:dTrans-1) = surv{j}{x,w}(1:dTrans-1) + 1;
                    case '31' % formed spine/bouton
                        transMat{j}{x}(q,1) = 7;
                    case '32' % formed filo
                        transMat{j}{x}(q,1) = 8;
                    case '33' % non-feature 
                        transMat{j}{x}(q,1) = 9;
                end
                %fprintf('  type = %d,  dur = %d\n', transMat{j}{x}(q,1), transMat{j}{x}(q,2) );
                % Distinguish transient and recurrent spines
                if w == 2 && ismember( transMat{j}{x}(q,1), [2,3] )
                    %fprintf('\n[j,x,w,q] = [%d, %d, %d, %d]  ', j,x,w,q);
                    if any( find( longMat{j}{x}(q,dTrans:end) == 1  ) )
                        %disp( [longMeta(j).qRec{x,w}, q] )
                        longMeta(j).qRec{x,w} = [longMeta(j).qRec{x,w}, q]; % fprintf('recurrent');
                    else
                        %disp( [longMeta(j).qTrans{x,w}, q] );
                        longMeta(j).qTrans{x,w} = [longMeta(j).qTrans{x,w}, q]; % fprintf('transient');
                    end
                end
                %disp( longMat{j}{x}(q,:) )
                %disp( transMat{j}{x}(q,:) )
            end
            longMeta(j).Nrec(x,w) = numel( longMeta(j).qRec{x,w} ); longMeta(j).Ntrans(x,w) = numel( longMeta(j).qTrans{x,w} );
            %Ntrans{j}{w}{x} = hist( transMat{j}{w}{x}(:,1), 9 );
             % Calculate overall survival
            %[~,dTemp,~] = intersect( longMeta(j).d{x}', longMeta(j).dUse ); dTemp = dTemp';
            totSurv{j}{w} = totSurv{j}{w} + surv{j}{x,w}; %totSurv{j}{w} = totSurv{j}{w} + surv{j}{x,w}(dTemp');
        end
        % Estimate the stability of spines that first appeared on the SECOND imaging session
        longMeta(j).qForm{x} = intersect( find((transMat{j}{x}(:,1) == 4 | transMat{j}{x}(:,1) == 7) & transMat{j}{x}(:,2) == 0), [longMeta(j).qType{x,[2,3]}] )'; % 2->1 and 3->1 transitions for dendritic features only
        %longMat{j}{x}(longMeta(j).qForm{x},:)
        longMeta(j).Nform(x) = numel(longMeta(j).qForm{x});
        longMeta(j).formDur{x} = nan(longMeta(j).Nform(x), 1); 
        c = 0;
        for q = longMeta(j).qForm{x}
            c = c+1;
            %longMat{j}{x}(q,1:end)
            dTrans = find( longMat{j}{x}(q,3:end) ~= 1, 1 ) + 2; % day when the spine disappeared. Add 2 to correct for exclusion of first 2 timepoints
            if isempty( dTrans )
                longMeta(j).formDur{x}(c) = longMeta(j).useDays(end) - longMeta(j).useDays(2);
            else
                longMeta(j).formDur{x}(c) = longMeta(j).useDays(dTrans-1) - longMeta(j).useDays(2);
            end
        end
        longMeta(j).qMatur{x} = find( transMat{j}{x}(:,1) == 4 )';  longMeta(j).Nmatur(x) = numel(longMeta(j).qMatur{x});
        for q = longMeta(j).qMatur{x}
            maturTemp = longMat{j}{x}(q,:);  dMat = find( maturTemp == 1 );
            matDur(q) = longMeta(j).days{x}(dMat(end)) - longMeta(j).days{x}(dMat(1));
        end
    end
    longMeta(j).totMatur = nansum( longMeta(j).Nmatur );  longMeta(j).maturFrac = longMeta(j).totMatur/nansum(longMeta(j).Ntype(:,3)); % 
    longMeta(j).totRec = nansum(longMeta(j).Nrec, 1);  longMeta(j).totTrans = nansum(longMeta(j).Ntrans, 1); % numbers of transient and recurrent spines
    
    % Fit the overall survival curve to the model defined above 
    if fitToggle
        for w = 1:3
            clearvars tempFit
            if totSurv{j}{w}(1) > 10
                [tempFit, tempGof] = fit( daysTemp', totSurv{j}{w}'/totSurv{j}{w}(1), modelSurv, ...
                    'StartPoint', modelStart{w}, 'display','final', 'MaxIter',1000, 'MaxFunEvals',1000, 'TolFun',10^-10, 'TolX',10^-8 );
                longMeta(j).pop(1,w) = tempFit.a; longMeta(j).pop(2,w) = 1-longMeta(j).pop(1,w); 
                longMeta(j).tau(1,w) = -1/tempFit.b; 
                if strcmpi( fitType, '2' ), longMeta(j).tau(2,w) = -1/tempFit.c; else, longMeta(j).tau(2,w) = NaN; end
                longMeta(j).R2(w) = tempGof.rsquare;
            else
                longMeta(j).pop(w) = NaN;  longMeta(j).imperm(w) = NaN;  longMeta(j).tau(w) = NaN;  longMeta(j).R2(w) = NaN;
            end
        end
    end
    %plot( f, x, y ); hold all;
    %{    
    if show
        % Plot the survival curves for boutons, spines and filopodia
        FS = 16;
        load( 'colorblind_colormap.mat' );
        mouseColor = distinguishable_colors(Nmouse);  %colorblind;  % 
        TL = [0.005,0];
        textPad = 0.08;
        close all; clearvars h
        SurvivalCurves = figure('units','normalized','outerposition',[0 0 1 1], 'color', 'w');
        for x = longMeta(j).x
            % Boutons
            plot( longMeta(j).days{x}-longMeta(j).days{x}(1), surv{j}{x,1}/surv{j}{x,1}(1), 'Color', mouseColor(j,:), 'LineStyle','-' ); hold on;
            plot( longMeta(j).days{x}-longMeta(j).days{x}(1), surv{j}{x,1}/surv{j}{x,1}(1), '.', 'Color', mouseColor(j,:), 'MarkerSize', 10); 
            text( longMeta(j).days{x}(end)-longMeta(j).days{x}(1)+textPad, (surv{j}{x,1}(end)/surv{j}{x,1}(1)), sprintf('%s (n = %d)', mouse(j).ID, surv{j}{x,1}(1) ), 'HorizontalAlignment','left', 'Color', mouseColor(j,:) );
            % Spines
            plot( longMeta(j).days{x}-longMeta(j).days{x}(1), surv{j}{x,2}/surv{j}{x,2}(1), 'Color', mouseColor(j,:), 'LineStyle','-.' ); %hold on;
            plot( longMeta(j).days{x}-longMeta(j).days{x}(1), surv{j}{x,2}/surv{j}{x,2}(1), '.', 'Color', mouseColor(j,:), 'MarkerSize', 10); 
            text( longMeta(j).days{x}(end)-longMeta(j).days{x}(1)+textPad, (surv{j}{x,2}(end)/surv{j}{x,2}(1)), sprintf('%s (n = %d)', mouse(j).ID, surv{j}{x,2}(1) ), 'HorizontalAlignment','left', 'Color', mouseColor(j,:) );
            % Filopodia
            plot( longMeta(j).days{x}-longMeta(j).days{x}(1), surv{j}{x,3}/surv{j}{x,3}(1), 'Color', mouseColor(j,:), 'LineStyle','--' ); %hold on;
            plot( longMeta(j).days{x}-longMeta(j).days{x}(1), surv{j}{x,3}/surv{j}{x,3}(1), '.', 'Color', mouseColor(j,:), 'MarkerSize', 10); 
            text( longMeta(j).days{x}(end)-longMeta(j).days{x}(1)+textPad, (surv{j}{x,3}(end)/surv{j}{x,3}(1)), sprintf('%s (n = %d)', mouse(j).ID, surv{j}{x,3}(1) ), 'HorizontalAlignment','left', 'Color', mouseColor(j,:) );
            ylabel('Survival', 'FontSize', FS); xlabel('Imaging Day');
            title( sprintf('%s-%s', mouse(j).ID, mouse(j).tiles{x} ), 'FontSize', FS)
            ylim([0,1.001]); xlim([-0.01,longMeta(j).days{x}(end)-longMeta(j).days{x}(1)+2*textPad]);
            %legend(h(Jlong), Slegend(Jlong), 'Location','SouthWest'); %  'NorthEastOutside'
            set(gca,'TickDir','out','box', 'off','TickLength',TL,'Xtick',longMeta(j).days{x}-longMeta(j).days{x}(1), 'FontSize', FS);
            pause; cla;
        end
        % Boutons
        h(1) = plot( longMeta(j).days{x}-longMeta(j).days{x}(1), totSurv{j}{1}/totSurv{j}{1}(1), 'Color', mouseColor(j,:), 'LineStyle','-' ); hold on;
        plot( longMeta(j).days{x}-longMeta(j).days{x}(1), totSurv{j}{1}/totSurv{j}{1}(1), '.', 'Color', mouseColor(j,:), 'MarkerSize', 10); 
        text( longMeta(j).days{x}(end)-longMeta(j).days{x}(1)+textPad, (totSurv{j}{1}(end)/totSurv{j}{1}(1)), sprintf('%s (n = %d)', mouse(j).ID, totSurv{j}{1}(1) ), 'HorizontalAlignment','left', 'Color', mouseColor(j,:) );
        % Spines
        h(2) = plot( longMeta(j).days{x}-longMeta(j).days{x}(1), totSurv{j}{2}/totSurv{j}{2}(1), 'Color', mouseColor(j,:), 'LineStyle','-.' ); %hold on;
        plot( longMeta(j).days{x}-longMeta(j).days{x}(1), totSurv{j}{2}/totSurv{j}{2}(1), '.', 'Color', mouseColor(j,:), 'MarkerSize', 10); 
        text( longMeta(j).days{x}(end)-longMeta(j).days{x}(1)+textPad, (totSurv{j}{2}(end)/totSurv{j}{2}(1)), sprintf('%s (n = %d)', mouse(j).ID, totSurv{j}{2}(1) ), 'HorizontalAlignment','left', 'Color', mouseColor(j,:) );
        % Filopodia
        h(3) = plot( longMeta(j).days{x}-longMeta(j).days{x}(1), totSurv{j}{3}/totSurv{j}{3}(1), 'Color', mouseColor(j,:), 'LineStyle','--' ); %hold on;
        plot( longMeta(j).days{x}-longMeta(j).days{x}(1), totSurv{j}{3}/totSurv{j}{3}(1), '.', 'Color', mouseColor(j,:), 'MarkerSize', 10); 
        text( longMeta(j).days{x}(end)-longMeta(j).days{x}(1)+textPad, (totSurv{j}{3}(end)/totSurv{j}{3}(1)), sprintf('%s (n = %d)', mouse(j).ID, totSurv{j}{3}(1) ), 'HorizontalAlignment','left', 'Color', mouseColor(j,:) );
        ylabel('Survival', 'FontSize', FS); title( 'Feature Survival', 'FontSize', FS)
        title( sprintf('%s: Overall', mouse(j).ID), 'FontSize', FS)
        ylim([0,1.001]); xlim([-0.01,longMeta(j).days{x}(end)-longMeta(j).days{x}(1)+2*textPad]);
        legend(h, {'Bouton','Spine','Filopodia'}, 'Location','SouthWest'); %  'NorthEastOutside'
        set(gca,'TickDir','out','box', 'off','TickLength',TL,'Xtick',longMeta(j).days{x}-longMeta(j).days{x}(1), 'FontSize', FS);
        pause;
        pdfPath = [mouse(j).figDir, mouse(j).ID, '-SurvivalCurves.pdf'];
        %export_fig( pdfPath, '-painters','-r600','-q101', SurvivalCurves ); fprintf( 'Saved %s \n',pdfPath );  
        close all;
        %close all; clearvars h
        %TransitionsPie = figure('units','normalized','outerposition',[0 0 1 1], 'color', 'w');
    end 
    %}
end
% Plot the survival curves for boutons, spines and filopodia across mice
if show
    close all; clearvars h
    SurvivalCurves = figure('units','normalized','outerposition',[0 0 1 1], 'color', 'w');
    sp(3) = subtightplot(3,1,3,opt{:}); sp(2) = subtightplot(3,1,2,opt{:}); sp(1) = subtightplot(3,1,1,opt{:});
    linkaxes(sp,'x');
    for j = setj
        daysTemp = longMeta(j).useDays-longMeta(j).useDays(1);
        % Boutons
        subplot(sp(1))
        h(1) = plot( daysTemp, totSurv{j}{1}/totSurv{j}{1}(1), 'Color', mouseColor(j,:) ); hold on;
        plot( daysTemp, totSurv{j}{1}/totSurv{j}{1}(1), '.', 'Color', mouseColor(j,:), 'MarkerSize', 10); 
        text( mouse(j).days(longMeta(j).dUse(end))-longMeta(j).useDays(1)+textPad, (totSurv{j}{1}(end)/totSurv{j}{1}(1)), sprintf('%s (n = %d)', mouse(j).ID, totSurv{j}{1}(1) ), 'HorizontalAlignment','left', 'Color', mouseColor(j,:), 'FontSize',FS/2 );
        if strcmpi( fitType, '1' ) && ~isnan(longMeta(j).tau(1,1))
            modelCurve = modelSurv( longMeta(j).pop(1,1), -1/longMeta(j).tau(1,1), daysTemp );
            plot( daysTemp, modelCurve, 'Color', mouseColor(j,:), 'LineWidth',1, 'LineStyle','--' );
        elseif strcmpi( fitType, '2' ) && ~isnan(longMeta(j).tau(1,1)+longMeta(j).tau(2,1))
            modelCurve = modelSurv( longMeta(j).pop(1,1), -1/longMeta(j).tau(1,1), -1/longMeta(j).tau(2,1), daysTemp );
            plot( daysTemp, modelCurve, 'Color', mouseColor(j,:), 'LineWidth',1, 'LineStyle','--' );
        end

        ylabel('Survival', 'FontSize', FS, 'FontWeight','bold'); title('Boutons', 'FontSize', FS);
        set(gca, 'Xtick',[],'TickDir','out','box', 'off','TickLength',TL, 'FontSize', FS);
        % Spines
        subplot(sp(2))
        h(2) = plot( daysTemp, totSurv{j}{2}/totSurv{j}{2}(1), 'Color', mouseColor(j,:) ); hold on;
        plot( daysTemp, totSurv{j}{2}/totSurv{j}{2}(1), '.', 'Color', mouseColor(j,:), 'MarkerSize', 10); 
        %text(mouse(j).days(longMeta(j).dUse(end))-longMeta(j).useDays(1)+textPad, (totSurv{j}{2}(end)/totSurv{j}{2}(1)), sprintf('%s (n = %d)', mouse(j).ID, totSurv{j}{2}(1) ), 'HorizontalAlignment','left', 'Color', mouseColor(j,:) );
        if strcmpi( fitType, '1' ) && ~isnan(longMeta(j).tau(1,2))
            modelCurve = modelSurv( longMeta(j).pop(1,2), -1/longMeta(j).tau(1,2), daysTemp );
            plot( daysTemp, modelCurve, 'Color', mouseColor(j,:), 'LineWidth',1, 'LineStyle','--' );
            tempText = sprintf('%s (P%d, n = %d, tau_1 = %2.2f days, perm = %2.2f, R^2 = %2.2f)', mouse(j).ID, mouse(j).ages(1), totSurv{j}{2}(1), longMeta(j).tau(1,2), longMeta(j).pop(2,2), longMeta(j).R2(1,2) );
            text( daysTemp(end), modelCurve(end), tempText, 'HorizontalAlignment','left', 'Color', mouseColor(j,:), 'FontSize',FS/2 );
        elseif strcmpi( fitType, '2' ) && ~isnan(longMeta(j).tau(1,2)+longMeta(j).tau(2,2))
            modelCurve = modelSurv( longMeta(j).pop(1,2), -1/longMeta(j).tau(1,2), -1/longMeta(j).tau(2,2), daysTemp );
            plot( daysTemp, modelCurve, 'Color', mouseColor(j,:), 'LineWidth',1, 'LineStyle','--' );
            tempText = sprintf('%s (P%d, n = %d, tau_1 = %2.2f days, tau_2 = %2.2f days, perm = %2.2f, R^2 = %2.2f)', mouse(j).ID, mouse(j).ages(1), totSurv{j}{2}(1), longMeta(j).tau(1,2), longMeta(j).tau(2,2),longMeta(j).pop(2,2), longMeta(j).R2(1,2) );
            text( daysTemp(end), modelCurve(end), tempText, 'HorizontalAlignment','left', 'Color', mouseColor(j,:), 'FontSize',FS/2 );
        end

        ylabel('Survival', 'FontSize', FS, 'FontWeight','bold'); title('Spines', 'FontSize', FS);
        set(gca,'Xtick',[],'TickDir','out','box', 'off','TickLength',TL, 'FontSize', FS);
        % Filopodia
        subplot(sp(3))
        h(3) = plot( daysTemp, totSurv{j}{3}/totSurv{j}{3}(1), 'Color', mouseColor(j,:) ); hold on;
        plot( daysTemp, totSurv{j}{3}/totSurv{j}{3}(1), '.', 'Color', mouseColor(j,:), 'MarkerSize', 10); 
        %text( mouse(j).days(longMeta(j).dUse(end))-longMeta(j).useDays(1)+textPad, (totSurv{j}{3}(end)/totSurv{j}{3}(1)), sprintf('%s (n = %d)', mouse(j).ID, totSurv{j}{3}(1) ), 'HorizontalAlignment','left', 'Color', mouseColor(j,:), 'FontSize',FS/2 );
        ylabel('Survival', 'FontSize', FS, 'FontWeight','bold'); title( 'Filopodia', 'FontSize', FS)
        xlabel('Days Imaged', 'FontSize', FS, 'FontWeight','bold');
        ylim([0,1.001]);  xlim([-0.01, 10*ceil(max([ longMeta.useDur ])/5)]);
        if strcmpi( fitType, '1' ) && ~isnan(longMeta(j).tau(1,3))
            modelCurve = modelSurv( longMeta(j).pop(1,3), -1/longMeta(j).tau(1,3), daysTemp );
            plot( daysTemp, modelCurve, 'Color', mouseColor(j,:), 'LineWidth',1, 'LineStyle','--' );
            tempText = sprintf('%s (P%d, n = %d, tau_1 = %2.2f days, perm = %2.2f, R^2 = %2.2f)', mouse(j).ID, mouse(j).ages(1), totSurv{j}{3}(1), longMeta(j).tau(1,3), longMeta(j).pop(2,3), longMeta(j).R2(1,3) );
            text( daysTemp(end), modelCurve(end), tempText, 'HorizontalAlignment','left', 'Color', mouseColor(j,:), 'FontSize',FS/2 );
        elseif strcmpi( fitType, '2' ) && ~isnan(longMeta(j).tau(1,3)+longMeta(j).tau(2,3))
            modelCurve = modelSurv( longMeta(j).pop(1,3), -1/longMeta(j).tau(1,3), -1/longMeta(j).tau(2,3), daysTemp );
            plot( daysTemp, modelCurve, 'Color', mouseColor(j,:), 'LineWidth',1, 'LineStyle','--' );
            tempText = sprintf('%s (P%d, n = %d, tau_1 = %2.2f days, tau_2 = %2.2f days, perm = %2.2f, R^2 = %2.2f)', mouse(j).ID, mouse(j).ages(1), totSurv{j}{3}(1), longMeta(j).tau(1,3), longMeta(j).tau(2,3),longMeta(j).pop(2,3), longMeta(j).R2(1,3) );
            text( daysTemp(end), modelCurve(end), tempText, 'HorizontalAlignment','left', 'Color', mouseColor(j,:), 'FontSize',FS/2 );
        end
        
        %set(gca, 'Xtick',[0,2,4,7,10,14:7:70], 'TickDir','out','box', 'off','TickLength',TL, 'FontSize', FS); % 'Xtick',longMeta(j).days{x}-longMeta(j).days{x}(1), daysTemp
        pause;
    end
    %figDir = [mouse(j).figDir,datestr(today,'yyyy-mm-dd'),slm]; mkdir(figDir);
    %pdfPath = [figDir,'SurvivalCurves.pdf'];
    %export_fig( pdfPath, '-painters','-r600','-q101', SurvivalCurves ); fprintf( 'Saved %s \n',pdfPath );  
    %close all;
    
    % Proportion of spines/filopodia undergoing transitions
    %{ 
    clearvars h; %close all; 
    figure('units','normalized','outerposition',[0 0 1 1], 'color', 'w');
    for j = 1:Nmouse
        %x = longMeta(j).x(1);
        % Boutons
        h(1) = plot( xTemp, totSurv{j}{1}/totSurv{j}{1}(1), 'Color', mouseColor(j,:), 'LineStyle','-' ); hold on;
        plot( xTemp, totSurv{j}{1}/totSurv{j}{1}(1), '.', 'Color', mouseColor(j,:), 'MarkerSize', 10); 
        text( mouse(j).days(longMeta(j).dUse(end))-longMeta(j).useDays(1)+textPad, (totSurv{j}{1}(end)/totSurv{j}{1}(1)), sprintf('%s (n = %d)', mouse(j).ID, totSurv{j}{1}(1) ), 'HorizontalAlignment','left', 'Color', mouseColor(j,:) );
        % Spines
        h(2) = plot( xTemp, totSurv{j}{2}/totSurv{j}{2}(1), 'Color', mouseColor(j,:), 'LineStyle','-.' ); %hold on;
        plot( xTemp, totSurv{j}{2}/totSurv{j}{2}(1), '.', 'Color', mouseColor(j,:), 'MarkerSize', 10); 
        text(mouse(j).days(longMeta(j).dUse(end))-longMeta(j).useDays(1)+textPad, (totSurv{j}{2}(end)/totSurv{j}{2}(1)), sprintf('%s (n = %d)', mouse(j).ID, totSurv{j}{2}(1) ), 'HorizontalAlignment','left', 'Color', mouseColor(j,:) );
        % Filopodia
        h(3) = plot( xTemp, totSurv{j}{3}/totSurv{j}{3}(1), 'Color', mouseColor(j,:), 'LineStyle','--' ); %hold on;
        plot( xTemp, totSurv{j}{3}/totSurv{j}{3}(1), '.', 'Color', mouseColor(j,:), 'MarkerSize', 10); 
        text( mouse(j).days(longMeta(j).dUse(end))-longMeta(j).useDays(1)+textPad, (totSurv{j}{3}(end)/totSurv{j}{3}(1)), sprintf('%s (n = %d)', mouse(j).ID, totSurv{j}{3}(1) ), 'HorizontalAlignment','left', 'Color', mouseColor(j,:) );
        ylabel('Survival', 'FontSize', FS); title( 'Filopodia Survival', 'FontSize', FS)
        ylim([0,1.001]); 
        legend(h, {'Bouton','Spine','Filopodia'}, 'Location','NorthEast'); %  'NorthEastOutside'
        set(gca,'TickDir','out','box', 'off','TickLength',TL, 'FontSize', FS); % 'Xtick',longMeta(j).days{x}-longMeta(j).days{x}(1),
        pause;
    end
    xlim([-0.01,Inf]);
    %}
end 
end