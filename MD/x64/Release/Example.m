clear all;clc;
disp('Running Simulation... not really just loading data')
data = importdata('HistogramData.txt',',',1);
%%

puredata = data.data;
%A = zeros(20, 10, 3);
%size(A)

disp('Creating histogram from data')
Hist = zeros(250,250,250);
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
%%
f_1 = figure;
imshow(Image)

export_fig(f_1,'1ion','-pdf','-nocrop','-transparent')

