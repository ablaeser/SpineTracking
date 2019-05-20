function feat = FeatProps( feat, featParam, varargin )
% FeatProps scans through all comparisons and collects details about each feature
% -----Inputs------
% feat 
% featParam = conversions from um per pixel/z-step
% featRad = 
%
% -----Outputs------
% feat
Nmouse = numel( feat );
IP = inputParser;
addRequired( IP, 'feat', @iscell )
addRequired( IP, 'featParam', @isstruct )
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'setj', flip(1:Nmouse), @isnumeric )
addParameter( IP, 'setk', [], @isnumeric )
addParameter( IP, 'setw', [], @isnumeric )
addParameter( IP, 'setx', [], @isnumeric )
addParameter( IP, 'setp', [], @isnumeric )
addParameter( IP, 'setq', [], @isnumeric )
parse( IP, feat, featParam, varargin{:} );

show = IP.Results.show;
setj = IP.Results.setj;
setk = IP.Results.setk; kToggle = isempty( setk ); % determine whether to set k manually or not
setw = IP.Results.setw; wToggle = isempty( setw );
setx = IP.Results.setx; xToggle = isempty( setx );
setp = IP.Results.setp; pToggle = isempty( setp );
setq = IP.Results.setq; qToggle = isempty( setq );
%
tic;
fprintf('Extracting feature properties...  ');
for j = setj
    fprintf('\nj = %d / %d  ', j, Nmouse );
    if kToggle
        setk = flip( 1:numel( feat{j} ) );
    end
    for k = setk
        %fprintf('\n    k = %d / %d ',k, numel(feat{j}));
        if xToggle
            setx = find(~cellfun(@isempty, feat{j}{k}(:,1)))'; % flip( 
        end
        for x = setx
            if pToggle
                setp = [2,1];
            end
            for p = setp
                if wToggle
                    setw = find( cellfun( @numel, feat{j}{k}(x,:,p) ) );
                end
                for w = setw 
                    if qToggle
                        setq = numel( feat{j}{k}{x,w,p} ):-1:1;
                    end
                    for q = setq % setq = 1:size(feat{j}{k}{x,w,p});
                        %fprintf('\n[j,k,x,p,w,q] = [%d, %d, %d, %d, %d, %d]',j,k,x,p,w,q)
                        %imshow( Feat.filt, [], 'InitialMagnification','fit' )
                        feat{j}{k}{x,w,p}(q) = getFeatProp( feat{j}{k}{x,w,p}(q), featParam ); % , minEcc connMap, 
                    end
                end
            end
        end
    end
    toc
end
fprintf('\n');
end