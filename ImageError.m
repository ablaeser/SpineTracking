function [ error, diffphase ] = ImageError( im, varargin )
%Calculate the translation invariant normalized RMS error

% Parse inputs
if nargin > 1
    xy = varargin{1}; % coordinates of fixed points to crop around: [y1,x1; y2,x2]
else
    xy = [round( size(im{1} )/2); round( size(im{2} )/2) ];
end
if any( isnan(xy(1,:)) ) 
    xy(1,:) = round( size(im{1} )/2);
end
if any( isnan(xy(2,:)) ) 
    xy(2,:) = round( size(im{2} )/2);
end
if nargin > 2
    show = varargin{2};
else
    show = false;
end
% Show the input images (optional)
if show
    figure;
    subplot(2,2,1);
    imshow( im{1}, []); hold on;
    plot( xy(1,2), xy(1,1), 'r*' )
    subplot(2,2,2);
    imshow( im{2}, []); hold on;
    plot( xy(2,2), xy(2,1), 'r*' )
    impixelinfo;
end
% If the images are unequal size, crop them to equal size, centered around their fixed points.
ydim = 1;
height = [size(im{1},ydim); size(im{2},ydim)];
ymargin = [xy(1,ydim), xy(2,ydim); height(1)-xy(1,ydim), height(2)-xy(2,ydim)]; % get the margins around the fixed points for each image [top1, bottom1; top2, bottom2]
if any( diff( ymargin )  ) %diff(height) ~= 0
    minYmargin = min( ymargin, [],  2 );
    im{2} = im{2}(xy(2,ydim)-minYmargin(1)+1:xy(2,ydim)+minYmargin(2),:);
    im{1} = im{1}(xy(1,ydim)-minYmargin(1)+1:xy(1,ydim)+minYmargin(2),:);
end
xdim = 2;
width = [size(im{1},xdim); size(im{2},xdim)];
xmargin = [xy(1,xdim), xy(2,xdim); width(1)-xy(1,xdim), width(2)-xy(2,xdim)];
if any( diff( xmargin )  ) % diff(width) ~= 0
    minXmargin = min( xmargin, [],  2 );
    im{2} = im{2}(:,xy(2,xdim)-minXmargin(1)+1:xy(2,xdim)+minXmargin(2));
    im{1} = im{1}(:,xy(1,xdim)-minXmargin(1)+1:xy(1,xdim)+minXmargin(2));
end
% Calculate error and phase difference
if isequal( size(im{1}), size(im{2}) )
    buf1ft = fft2( im{1} ); 
    buf2ft = fft2( im{2} );
    CCmax = sum(sum(buf1ft.*conj(buf2ft))); 
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2); 
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero); 
    error = sqrt(abs(error));
    diffphase = atan2(imag(CCmax),real(CCmax)); 
else
    error = NaN; diffphase = NaN;
end

if show
    subplot(2,2,1);
    title(sprintf('Error = %2.2f',error));
    subplot(2,2,3)
    imshow(im{1}, []);
    subplot(2,2,4)
    imshow(im{2}, []);
    impixelinfo;
end
end

%{


Itest = Istack{1}(:,:,34);
Iwin = uint16( window.*double(Itest));
figure;
subplot(1,4,1);
imshow( Itest, [] )
subplot(1,4,2)
Wham = hamming(512)*hamming(512).';
imshow( uint16( Wham.*double(Itest)), [] )
subplot(1,4,3)
Whan = hanning(512)*hanning(512).';
imshow( uint16( Whan.*double(Itest)), [] );
subplot(1,4,4)
Wblack = blackman(512)*blackman(512).';
imshow( uint16( Wblack.*double(Itest)), [] );

impixelinfo;
%}