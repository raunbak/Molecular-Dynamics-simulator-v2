clear all;clc;
disp('Running Simulation... not really just loading data')
data = importdata('HistogramData.txt',',',1);
%% Change variables here
sizeX = 250;
sizeY = 250;
sizeZ = 250;

%% Creating histogram 
puredata = data.data;

disp('Creating histogram from data')
Hist = zeros(sizeX,sizeY,sizeZ);
size(Hist)
n = 1;
for i = 1:size(Hist,1);
    for j = 1:size(Hist,2);
       for k = 1:size(Hist,3);
           
           
            Hist(i,j,k) = puredata(n);
            n = n + 1;
       end
    end
end

disp('Creating image')
Image = HistogramToImageBlur(Hist);
disp('Adjusting gray scale')
Image = mat2gray(Image);
%% Showing the image.
f_1 = figure;
imshow(Image)
