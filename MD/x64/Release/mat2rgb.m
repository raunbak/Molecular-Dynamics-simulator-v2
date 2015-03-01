function rgb = mat2rgb(X, map, numel)
if nargin < 3
    numel = 64;
end

if nargin == 1
    m = jet(numel);
    colormap(m);
    map = colormap('jet');
else
    m = jet(numel);
    colormap(m);
    map = colormap(map);
end

g = mat2gray(X);
x = gray2ind(g, size(map,1));
rgb = ind2rgb(x,map);
end