clear all;clc; close all;
disp('Data load')
data = importdata('TemperatureData.txt',',',1);
%% Change values here
T = 0.009337;
dim = 2; % dont use 1 
%- 2 is Total T, 3 is Tz, 4 is Tx, 5 is Ty.

A = 0;    %Start index of plot, keep at zero if you want all of the plot
L = 4285; % Number of cycles in plot. Rule of tumb should be TimeSteps/StepsPrRFperiode

%Low is when you want to start calculating temperature,get it to fit with
%simulating using the same rule of tumb.
Low = L-2000; 
High = L;

%% Create plot
Steps = 105; % Number of steps in RF periode.
f_1 = figure;
hold on
set(gca,'FontSize',12)

puredata = data.data;
Total = [];
count = 0;

for i = A:L
    plot(i,puredata(round(Steps/2)+Steps*i,dim),'.') 
    %This is same rule as used in the program, looking at the data past
    %1000 RF periodes, this should hit the lowest point in the periodes
    %everytime.
    
    %plot([round(Steps/2)+Steps*i round(Steps/2)+Steps*i]',[0 40],'g')
    
    if(i == Low || i == High)
        plot([i i],[0 2],'k')
    end
    if (i >= Low && i<= High)
    Total = [Total  puredata(round(Steps/2)+Steps*i,dim)];
    count = count + 1;
    end
end

plot(A:L,T*ones(1,L-A+1),'r');


AvgTemp = sum(Total) / count;
AvgTempstd = std(Total);


AvgTemp
AvgTempstd
procentile_off = (T - AvgTemp) / T * 100

xlabel('Tid [Cycles]');
ylabel('T [K]');
title('Temperature of crystal')
hold off


