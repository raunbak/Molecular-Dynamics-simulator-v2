clear all;clc;
disp('Running Simulation... not really just loading data')
data = importdata('CountHistogramData.txt',',',1);
Vdata = importdata('VelocityHistogramData.txt',',',1);
%%

simT = 0.02170;

puredata = data.data;
pureV = Vdata.data;
%A = zeros(20, 10, 3);
%size(A)

disp('Creating histogram from data')
Hist = zeros(250,250,250);
VHist = zeros( max(pureV(:,2))+1, max(pureV(:,3))+1, max(pureV(:,4))+1);


size(Hist)
n = 1;
for i = 1:size(Hist,1);
    for j = 1:size(Hist,2);
       for k = 1:size(Hist,3);
           
            Hist(i,j,k) = puredata(n);
            VHist(i,j,k) = pureV(n);
            n = n + 1;
       end
    end
end
%%
Counts = zeros(size(Hist,1)/2 -1,1);
Vel = zeros(size(Hist,1)/2 -1,1);
Howmany = zeros(size(Hist,1)/2,1);
Radii = zeros(size(Hist,1)/2,1);
r = 1;
n = 1;
while ( r < 125)
    
    
    for i = 125-r:r+125;
        for j = 125-r:r+125;
            for k = 125-r:r+125;
           
                d = sqrt( (i-125)^2 + (j-125)^2 + (k-125)^2 );
                if( d <= r  && d > r - 1)
                    
                    Counts(n) = Counts(n) +  Hist(i,j,k) ;
                    Vel(n) = Vel(n) + VHist(i,j,k);
                    
                    
                    Howmany(n) = Howmany(n) + 1;
                    
                    
                end
                
            end
        end
    end
    
    if(Counts(n) == 0)
         Vel(n) =0;
    else
         Vel(n) = Vel(n) / Counts(n); %* 1/(r* 0.89e-6);
    end
    %if ( n > 1)
    %Counts(n) = Counts(n) - Counts(n-1);
    %end
    
    Radii(n) = (r*0.89e-6);
    
    n = n + 1;
    r = r + 1
    
end

%%
massInu = 40;
Tokg = massInu * 1.66053878283e-27;
kb = 1.380650424e-23;
T = ( Tokg*Vel) / (3*kb);

%%
f_1 = figure;
hold on
set(gca,'FontSize',12)
%title('\Gamma_p = 200')
%plot(Radii(1:end-1),T,'k')
plot(Radii(1:end-1),T,'xk')
xlabel('Radius [m]');
ylabel('T [K]');
plot([8.96*10^(-5) 8.96*10^(-5)]',[0 0.025],'b--')
plot([0 0.12e-3]',[simT simT],'r')
%axis([0 0.12e-3 0 0.025])
legend('Data','Radius af krystal fra MD','Simuleret T','Location','Best')

hold off
export_fig(f_1,'Vel','-pdf','-nocrop','-transparent')
