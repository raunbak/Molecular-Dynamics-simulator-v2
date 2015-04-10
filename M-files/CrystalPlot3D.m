% For use with two ion species.

Data = importdata('MDpos.xyz','\t',2);
points = Data.data;
hold on
scatter3(points(1:end-100,3),points(1:end-100,1),points(1:end-100,2),'.')
scatter3(points(end-100:end,3),points(end-100:end,1),points(end-100:end,2),'r.')
axis('equal')

%% For use with one ion species
figure
Data = importdata('MDpos.xyz','\t',2);
points = Data.data;
hold on
scatter3(points(1:end,3),points(1:end,1),points(1:end,2),'.')
axis('equal')