function [ stackOut, shift, varargout ] = StackReg( stackIn, varargin )
% StackReg Performs dftregistration on a stack. Can be set to run either  
% -----Inputs------
% stackIn = stack to register
% zref = which frame to use as reference. If set to 0 (default), each frame is registered to the previous frame. If set to a negative value, each frame is registered to the middle frame. 
%       If set to a positive integer, each frame registered to that frame.
% upscale = 1/(desired resolution of registration). Setting upscale = 10 (default) yields accuracy to 0.1 pixels. See dftregistration for details.
% save = path to save stackOut to (default = '', in which case the stack is not saved.)
% -----Outputs------
% stackOut = cell array of registered images
% shift = (Nframe x 4)array of [ xshift, yshift, shift distance, error ] 
% stackMean = mean projection of registered stack
% stackMax = max projection of registered stack
 
stackCheck = @(x)( any([isnumeric(x(:)), ischar(x(:))] ) );
IP = inputParser;
addRequired( IP, 'stackIn', stackCheck )
addOptional( IP, 'zref', 0, @isnumeric )
addParameter( IP, 'crop', nan(1,4), @isnumeric )
addParameter( IP, 'upscale', 10, @isnumeric )
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'save', '', @ischar )
parse( IP, stackIn, varargin{:} );
zref = IP.Results.zref;
crop = IP.Results.crop;
upscale = IP.Results.upscale;
show = IP.Results.show;
savePath = IP.Results.save;
% crop = nan(1,4); upscale = 10; show = true;

% If stackIn is a filepath, load that file
if ischar( stackIn )
    stackIn = loadtiff(stackIn);
end
% crop the stack (optional)
if all( ~isnan( crop ) ) 
    stackIn = stackIn( crop(2):crop(2)+crop(4), crop(1):crop(1)+crop(3), : );
end
[resX, resY, Nframe] = size(stackIn);

% Perform registration 
stackOut = nan(size(stackIn)); shift = nan(Nframe,4); %zmid = NaN;
if zref >= 0 % Fixed-frame registration (default)
    if zref == 0
        zref = ceil( Nframe/2 );
    end
    fftref = fft2( stackIn(:,:,zref) ); 
    for z = 1:Nframe
        fftsource = fft2( stackIn(:,:,z) );
        [output, fftreg] = dftregistration( fftref, fftsource, upscale);
        stackOut(:,:,z) = abs( ifft2(fftreg) );
        shift(z,1) = output(4); shift(z,2) = output(3); shift(z,3) = norm(output(3:4)); shift(z,4) = output(1); 
    end
elseif isnan(zref) % Frame-by-frame registration starting from the first frame
    stackOut(:,:,1) = stackIn(:,:,1);  shift(1,:) = [0, 0, 0, 0];
    for z = 2:Nframe
        fftref = fft2( stackIn(:,:,z-1) ); 
        fftsource = fft2( stackIn(:,:,z) );
        [output, fftreg] = dftregistration( fftref, fftsource, upscale);
        stackOut(:,:,z) = abs( ifft2(fftreg) );
        shift(z,1) = output(4); shift(z,2) = output(3); shift(z,3) = norm(output(3:4)); shift(z,4) = output(1); 
        % Show the process (optional)
        %{
        if show
            subplot( sp(1) );  imshow( stackIn(:,:,z-1), [] )
            title( sprintf('Reference image: z = %d', z-1) );
            subplot( sp(2) );  imshow( stackIn(:,:,z), [] );
            title( sprintf('Comparison image: z = %d / %d', z,Nframe) );
            subplot( sp(3) );  imshow( stackOut(:,:,z), [] ); impixelinfo;
            title( sprintf('Registered image: shift = [%2.2f, %2.2f, %2.2f, %2.2f]',shift(z,1), shift(z,2), shift(z,3), shift(z,4)) );
            subplot( sp(4) );
            imshowpair( stackIn(:,:,z), stackOut(:,:,z), 'falsecolor' );
            title('Original vs. Registered');
            pause;
        end
        %}
    end
elseif zref < 0 % Frame-by-frame registration starting from the middle frame
    zmid = round( median( 1:Nframe ) );
    stackOut(:,:,zmid) = stackIn(:,:,zmid);  shift(zmid,:) = [0, 0, 0, 0];
    for z = zmid+1:Nframe
        fftref = fft2( stackIn(:,:,z-1) ); 
        fftsource = fft2( stackIn(:,:,z) );
        [output, fftreg] = dftregistration( fftref, fftsource, upscale);
        stackOut(:,:,z) = abs( ifft2(fftreg) );
        shift(z,1) = output(4); shift(z,2) = output(3); shift(z,3) = norm(output(3:4)); shift(z,4) = output(1); 
        % Show the process (optional)
        %{
        if show
            subplot( sp(1) );  imshow( stackIn(:,:,z-1), [] )
            title( sprintf('Reference image: z = %d', z-1) );
            subplot( sp(2) );  imshow( stackIn(:,:,z), [] );
            title( sprintf('Comparison image: z = %d / %d', z,Nframe) );
            subplot( sp(3) );  imshow( stackOut(:,:,z), [] ); impixelinfo;
            title( sprintf('Registered image: shift = [%2.2f, %2.2f, %2.2f, %2.2f]',shift(z,1), shift(z,2), shift(z,3), shift(z,4)) );
            subplot( sp(4) );
            imshowpair( stackIn(:,:,z), stackOut(:,:,z), 'falsecolor' );
            title('Original vs. Registered');
            pause;
        end
        %}
    end
    for z = zmid-1:-1:1
        fftref = fft2( stackIn(:,:,z+1) ); 
        fftsource = fft2( stackIn(:,:,z) );
        [output, fftreg] = dftregistration( fftref, fftsource, upscale);
        stackOut(:,:,z) = abs( ifft2(fftreg) );
        shift(z,1) = output(4); shift(z,2) = output(3); shift(z,3) = norm(output(3:4)); shift(z,4) = output(1); 
    end
end
stackOut = cast( stackOut, 'like', stackIn ); % Cast stackOut to same class as stackIn

% get weighted average and max projections over the registered frames
if nargout > 2
    regWeight = reshape( 1-shift(:,4), 1, [] ); 
    totWeight = sum(regWeight);
    varargout{1} = uint16( sum(repmat(reshape(regWeight,1,1,Nframe),resX,resY).*double(stackOut),3)/totWeight ); 
    varargout{2} = uint16( max(stackOut,[],3) );
end

% Save stackOut to tif file (optional)
if ~isempty( savePath )
    mkdir( fileparts( savePath ) );
    saveastiff( stackOut, savePath ); fprintf('Saved %s \n', savePath );
end
% Track the shift across frames (optional)
if show
    %{
    MS = 10;
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,[1,3])
    plot( shift(:,1), shift(:,2) ); hold on;
    plot( shift(1,1), shift(1,2), 'g.', 'MarkerSize', MS ); plot( shift(end,1), shift(end,2), 'r.', 'MarkerSize', MS );
    xlabel('X-Shift (pix)'); ylabel('Y-Shift (pix)'); axis square;
    subplot(2,2,2);
    plot( shift(:,3) );
    ylabel('Total Shift (pix)'); % xlabel('Frame'); 
    subplot(2,2,4);
    plot( shift(:,4) );
    xlabel('Frame'); ylabel('Error');
    %}
    opt = { [0.04,0.02], [0.05,0.03], [0.03,0.01]}; % {gapVector (between subplots) [vert, horz], margin height [bottom, top], margin width [left, right] } 
    figure('units','normalized','outerposition',[0 0 1 1]);
    for z = 1:Nframe
        subtightplot( 2, Nframe, z, opt{:} );
        imshow( stackIn(:,:,z), [] );
    end
    for z = 1:Nframe
        subtightplot( 2, Nframe, z+Nframe, opt{:} );
        imshow( stackOut(:,:,z), [] );
        title(sprintf('[dX, dY, dist, err] = [%2.1f, %2.1f, %2.1f, %2.1f]', shift(z,1), shift(z,2), shift(z,3), shift(z,4)));
    end
    impixelinfo;
end
%toc;
end
