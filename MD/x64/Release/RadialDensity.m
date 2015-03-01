clear all;clc;
disp('Running Simulation... not really just loading data')
data = importdata('HistogramData.txt',',',1);
%%
puredata = data.data;
%A = zeros(20, 10, 3);
%size(A)

disp('Creating histogram from data')
Hist = zeros(200,200,520);
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

%% Picking out the slice and summing

Slice = Hist(:,:,110:130);
SumSlice = sum(Slice,3);



%%
Counts = zeros(50,1);
r = 1;
n = 1;
while ( r < 100)
     
    for i = 100-r:r+100;
        for j = 100-r:r+100;
           
                d = sqrt( (i-100)^2 + (j-100)^2);
                if( d <= r  && d >= r - 2)  
                    Counts(n) = Counts(n) +  SumSlice(i,j) ;
                end
        end
    end

    %if ( n > 1)
    %Counts(n) = Counts(n) - Counts(n-1);
    %end
    n = n + 1;
    r = r + 2
    
end

R = 1:2:size(Hist,1)/2;
r = 0:2:((size(Hist,1)/2)-1);
V = 2* pi * 2* R';

Density = Counts ./ V ;

%% Plotting 
f_1 = figure;
hold on
%plot(R,Density,'k')
%plot(R,Density,'xk')
bar(R,Density,'k')
xlabel('Radius [bins]');
ylabel('Time aveage density [arb]');

%plot([100.6742 100.6742],[0 40],'r')

%plot([94 94]',[0 40],'r')
axis([0 100 0 1800])

hold off
export_fig(f_1,'Sphere57','-pdf','-nocrop','-transparent')