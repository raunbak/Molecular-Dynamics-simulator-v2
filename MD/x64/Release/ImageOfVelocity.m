clear all;clc;
disp('Running Simulation... not really just loading data')
data = importdata('CountHistogramData.txt',',',1);
Vdata = importdata('VelocityHistogramData.txt',',',1);
%%

puredata = data.data;
pureV = Vdata.data;

disp('Creating histogram from data')
Hist = zeros(250,250,250);
VHist = zeros( max(pureV(:,2))+1, max(pureV(:,3))+1, max(pureV(:,4))+1);


size(Hist)
n = 1;
for i = 1:size(Hist,1);
    for j = 1:size(Hist,2);
       for k = 1:size(Hist,3);
           
            if (puredata(n) ~= 0)
                Hist(i,j,k) = pureV(n)/puredata(n);
            else
                Hist(i,j,k) = 0;    
            end
            n = n + 1;
       end
    end
end

%%
Gap = 120:130;
Slice = Hist(:,:,Gap);
SumSlice = zeros(250,250);
for i = 1:250
   for j = 1:250
       for k = 1:11
        SumSlice(i,j) = SumSlice(i,j) + Slice(i,j,k);   
       end
   end
end

SumSlice = SumSlice/11;
%sum(Slice,3)/length(Gap);
%%
%disp('Creating image')
%Image = HistogramToImageBlur(Hist);
%disp('Adjusting gray scale')
Image = mat2gray(SumSlice);
%%
f_1 = figure;
hold on 
set(gca,'FontSize',15)
imshow(Image)
xlabel('y-retning')
ylabel('x-retning')
colormap(jet(256));
colorbar('FontSize',15)