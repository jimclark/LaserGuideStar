nstars = 259;

idxs = 1:nstars;
phi = acos(1-2.*(idxs-0.5)./nstars);
theta = pi*(1+sqrt(5))*(idxs-0.5);
x = 1.0*cos(theta).*sin(phi);
y = 1.0*sin(theta).*sin(phi);
z = 1.0*cos(phi);

xyz = [x' y' z'];

resultStruct = tsp_nn('xy',xyz);

%%
resultStruct2 = tsp_ga(resultStruct);

%%
save('tsp_fib.mat','resultStruct','resultStruct2');

%%
srcs = resultStruct2.optSolution(1:nstars);
dests = resultStruct2.optSolution(2:nstars+1);

dists = sqrt((x(srcs)-x(dests)).^2 + (y(srcs)-y(dests)).^2) + (z(srcs)-z(dests)).^2;

disp(rad2deg(mean(dists)))



figureGlobe2 = figure;
axesm('globe','Grid','on','GLineWidth',2,'MeridianLabel','on','MLabelParallel','equator','ParallelLabel','on','PLabelMeridian','prime')
set(gca, 'fontsize', 14,'linewidth',2)

plot3(x,y,z,'*','LineWidth',2)

for i = 1:nstars
    plot3(1.02*[x(srcs(i)) x(dests(i))],1.02*[y(srcs(i)) y(dests(i))],1.02*[z(srcs(i)) z(dests(i))],'k-','LineWidth',2);
end

base = zeros(180,360);
baseR = georefcells([-90 90],[0 360],size(base));
copperColor = [0.62 0.38 0.24];
geoshow(base,baseR,'FaceColor',copperColor)
camlight right
material([.8 .9 .4])

%%
saveas(figureGlobe2,'skyglobe2.png')