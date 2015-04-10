clear all;clc;
disp('Running Simulation... not really just loading data')
disp('This file should only be used on crystal spheres')

%% The following data files are used.
data = importdata('CountHistogramData.txt',',',1);
Vdata = importdata('VelocityHistogramData.txt',',',1);
data2 = importdata('IonData.txt',',',2); %Import of radius of krystal.
a = data2.textdata(2);
b = char(a);
b = str2num(b(1,4:15));
%% Change stuff here for different simulations

% Value of simulated temperature
simT = 0.02170;
massInu = 40;

SizeXofHistogram = 250;
SizeYofHistogram = 250;
SizeZofHistogram = 250;
PixelLength = 0.89e-6;

RadiusOfKrystalCalculatedInProgram = b;%8.96*10^(-5);

%% Create histograms in matlab
puredata = data.data;
pureV = Vdata.data;

disp('Creating histogram from data')
Hist = zeros(SizeXofHistogram,SizeYofHistogram,SizeZofHistogram);
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
%% Finding how many are at each simple sphere-shel at radius r
Counts = zeros(size(Hist,1)/2 -1,1);% Summed total number of ions at radius index
Vel = zeros(size(Hist,1)/2 -1,1);   % Summed total vel at this radius index
Howmany = zeros(size(Hist,1)/2,1);  % Array for how many ions at a given radius index
Radii = zeros(size(Hist,1)/2,1);    % Array for radius length 

r = 1;
n = 1;

% Loop for filling out arrays. 
%Could be done alot cleaner/faster with marking indexes and such. 
while ( r < 125)
    
    for i = 125-r:r+125;
        for j = 125-r:r+125;
            for k = 125-r:r+125;
           
                d = sqrt( (i-125)^2 + (j-125)^2 + (k-125)^2 ); %Distance of ion
                if( d <= r  && d > r - 1)
                    
                    Counts(n) = Counts(n) +  Hist(i,j,k) ;
                    Vel(n) = Vel(n) + VHist(i,j,k);
                    
                    
                    Howmany(n) = Howmany(n) + 1;
                    
                    
                end
                
            end
        end
    end
    
    %Since we dont check the same radius again we normalize here.
    if(Counts(n) == 0)
         Vel(n) =0;
    else
         Vel(n) = Vel(n) / Counts(n); %* 1/(r* 0.89e-6);
    end
    
    Radii(n) = (r*PixelLength);
    
    n = n + 1;
    r = r + 1
    
end

%% Convert to Kelvin.

Tokg = massInu * 1.66053878283e-27;
kb = 1.380650424e-23;
T = ( Tokg*Vel) / (3*kb);

%% Plot
f_1 = figure;
hold on
set(gca,'FontSize',12)
plot(Radii(1:end-1),T,'xk')
xlabel('Radius [m]');
ylabel('T [K]');
plot([RadiusOfKrystalCalculatedInProgram RadiusOfKrystalCalculatedInProgram]',[0 simT*2],'b--') % Radius of crystal from program.
plot([0 0.12e-3]',[simT simT],'r')

legend('Data','Radius of crystal from MD','Sim T','Location','Best')
hold off
