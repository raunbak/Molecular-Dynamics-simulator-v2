function ImagePlane = HistogramToImageBlur(Histogram)

OffSet = 0;
% rescaling values

Histogram = Histogram/max(max(max(Histogram)))*100;

RotatedHistogram = imrotate(reshape(Histogram(:,:,1),size(Histogram,1),size(Histogram,2)),45,'bilinear','crop');

RotatedHistogram = zeros(size(RotatedHistogram,1),size(RotatedHistogram,2),size(RotatedHistogram,3));

for Znumber = 1:size(Histogram,3)
    RotatedHistogram(:,:,Znumber) = imrotate(reshape(Histogram(:,:,Znumber),size(Histogram,1),size(Histogram,2)),45,'bilinear','crop');
end
Histogram = RotatedHistogram;

% Blurring and Projecting to 2d plane
ImagePlane = zeros(size(Histogram,1), size(Histogram,3));

BinLength = 0.89;
lambda = 397; % flourecent wave length in nm
%w0 = 2.87; % value from image calibration in microns
%w0 = 2.4;
w0=1.7;
z0 = (pi*w0^2/(lambda*10^(-3)));

disp('projecting and blurring')

for zPlane = 1:size(Histogram,2)
    Zdistance = abs(zPlane*BinLength-BinLength*size(Histogram,2)/2);
   % blur = w0*sqrt(1+Zdistance^2/z0^2)*sqrt(2)-w0+OffSet;
    blur = w0*sqrt(1+Zdistance^2/z0^2)/sqrt(2);
    H = fspecial('gaussian', 480, blur);
    BlurredPlane = imfilter(reshape(Histogram(:,zPlane,:),size(Histogram,1),size(Histogram,3)),H,'replicate');
    ImagePlane=ImagePlane+BlurredPlane;
end

disp('Done blurring...')

