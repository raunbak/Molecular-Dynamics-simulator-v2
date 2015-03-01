clear all;clc; close all;
disp('Data load')
data = importdata('TemperatureData.txt',',',1);
%%
Steps = 105;
T = 0.02170;
puredata = data.data;
dim = 2;
%%
f_1 = figure;
hold on
set(gca,'FontSize',12)
%plot(puredata(:,dim),'.')

%%
Total = [];
count = 0;
A = 0%11000;
%L = 10000;
L = 4285;%22289;%8707;

Low = L-2000;%21507;%6007;
 High = L;%23420;%8007;
for i = A:L
    plot(i,puredata(round(Steps/2)+Steps*i,dim),'.')
    %plot([round(Steps/2)+Steps*i round(Steps/2)+Steps*i]',[0 40],'g')
    
    if(i == Low || i == High)
        plot([i i],[0 2],'k')
    end
    if (i >= Low && i<= High)
    Total = [Total  puredata(round(Steps/2)+Steps*i,dim)];
    count = count + 1;
    end
end
%plot(puredata(:,1),T*ones(length(puredata(:,1)),1),'r');
plot(A:L,T*ones(1,L-A+1),'r');

%LowestPartOfPeriode = floor((Steps/2) + 0.5);
%
%for t = 1:24830
%    if ( mod(t,LowestPartOfPeriode) == 0 && mod(t,Steps) ~= 0)
%    plot([t t],[0 100])
%    end
%end
%axis([0 5*10^4 0 200])
AvgTemp = sum(Total) / count;
AvgTempstd = std(Total);


AvgTemp
AvgTempstd
procent_afvigelse = (T - AvgTemp) / T * 100
%text(2.442e4,78,['Procent afvigelse ',num2str(procent_afvigelse),'%'])

%axis([19507 21420 0 1.5])
xlabel('Tid [Skridt]');
ylabel('T [K]');
title('Temperatur af krystal')
hold off

%export_fig(f_1,'Temperatur','-pdf','-nocrop','-transparent')


%%
