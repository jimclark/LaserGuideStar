% Run StarkSkymap first to initialize starlats and starlons
close all;

nstars = numel(starlats);

sc_max_acc = sc_max_thrust_nom/sc_mass_opt_tot;

idxs = nchoosek(1:nstars,2);

[arclens az] = distance(starlats(idxs(:,1)),starlons(idxs(:,1)),...
    starlats(idxs(:,2)),starlons(idxs(:,2)));
arcrads = deg2rad(arclens);

costs = 2*sqrt(range_LGS.*sc_max_acc.*arcrads);

lendist = length(costs);

G = graph(idxs(:,1),idxs(:,2));

figureChart = figure;
hGraph = plot(G,'XData',starlons-360.*(starlons>180),'YData',starlats,'LineStyle','none','NodeLabel',{});

figureMap = figure;
axesm('MapProjection','robinson','Grid','on','GLineWidth',2)
p1 = scatterm(starlats,starlons,'*', 'linewidth', 2,'DisplayName','Stark 2015 targets');

%%

Aeq = spalloc(nstars,length(idxs),nstars*(nstars-1)); % Allocate a sparse matrix
for ii = 1:nstars
    whichIdxs = (idxs == ii); % Find the trips that include stop ii
    whichIdxs = sparse(sum(whichIdxs,2)); % Include trips where ii is at either end
    Aeq(ii,:) = whichIdxs'; % Include in the constraint matrix
end
beq = 2*ones(nstars,1);

intcon = 1:lendist;
lb = zeros(lendist,1);
ub = ones(lendist,1);
%%
opts = optimoptions('intlinprog');
[x_tsp,costopt,exitflag,output] = intlinprog(costs,intcon,[],[],Aeq,beq,lb,ub,opts);

x_tsp = logical(round(x_tsp));
Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2));

%%
highlight(hGraph,Gsol,'LineStyle','-')

tourIdxs = conncomp(Gsol);
numtours = max(tourIdxs); % number of subtours
fprintf('# of subtours: %d\n',numtours);

%%

A = spalloc(0,lendist,0); % Allocate a sparse linear inequality constraint matrix
b = [];
while numtours > 1 % Repeat until there is just one subtour
    % Add the subtour constraints
    b = [b;zeros(numtours,1)]; % allocate b
    A = [A;spalloc(numtours,lendist,nstars)]; % A guess at how many nonzeros to allocate
    for ii = 1:numtours
        rowIdx = size(A,1) + 1; % Counter for indexing
        subTourIdx = find(tourIdxs == ii); % Extract the current subtour
%         The next lines find all of the variables associated with the
%         particular subtour, then add an inequality constraint to prohibit
%         that subtour and all subtours that use those stops.
        variations = nchoosek(1:length(subTourIdx),2);
        for jj = 1:length(variations)
            whichVar = (sum(idxs==subTourIdx(variations(jj,1)),2)) & ...
                       (sum(idxs==subTourIdx(variations(jj,2)),2));
            A(rowIdx,whichVar) = 1;
        end
        b(rowIdx) = length(subTourIdx) - 1; % One less trip than subtour stops
    end

    % Try to optimize again
    [x_tsp,costopt,exitflag,output] = intlinprog(costs,intcon,A,b,Aeq,beq,lb,ub,opts);
    x_tsp = logical(round(x_tsp));
    Gsol = graph(idxs(x_tsp,1),idxs(x_tsp,2));
    
    % Visualize result
    hGraph.LineStyle = 'none'; % Remove the previous highlighted path
    highlight(hGraph,Gsol,'LineStyle','-')
    drawnow
    
    % How many subtours this time?
    tourIdxs = conncomp(Gsol);
    numtours = max(tourIdxs); % number of subtours
    fprintf('# of subtours: %d\n',numtours)
end

%%

save('stark_skymap_tsp.mat','Gsol')

%%

dists = 1:numel(Gsol.Edges);

for i = 1:numel(Gsol.Edges)
    dists(i) = distance(starlats(Gsol.Edges{i,1}(1)),starlons(Gsol.Edges{i,1}(1)),...
        starlats(Gsol.Edges{i,1}(2)),starlons(Gsol.Edges{i,1}(2)));
end

dists = deg2rad(dists);

times = 2*sqrt(range_LGS.*dists./sc_max_acc);

days = times/86400;

dvs = times.*sc_max_acc;

%%


figureGlobe = figure;
axesm('globe','Grid','on','GLineWidth',2,'MeridianLabel','on','MLabelParallel','equator','ParallelLabel','on','PLabelMeridian','prime')
p1 = scatterm(starlats,starlons,'*', 'linewidth', 2,'DisplayName','Stark 2015 targets');
p2 = scatterm(deeplats,deeplons,'rv', 'linewidth', 2);
p3 = scatterm(brightlats,brightlons,'g+', 'linewidth', 2);
% legend([p1 p3 p2],{'Stark 2015 targets','Magnitude 2 stars','Hubble/Chandra deep fields'})
set(gca, 'fontsize', 14,'linewidth',2)

for i = 1:numel(Gsol.Edges)
    plot3m([starlats(Gsol.Edges{i,1}(1)) starlats(Gsol.Edges{i,1}(2))],...
        [starlons(Gsol.Edges{i,1}(1)) starlons(Gsol.Edges{i,1}(2))],...
        [0.02 0.02],'k-','LineWidth',2);
end

base = zeros(180,360);
baseR = georefcells([-90 90],[0 360],size(base));
copperColor = [0.62 0.38 0.24];
geoshow(base,baseR,'FaceColor',copperColor)
camlight right
material([.8 .9 .4])
saveas(figureGlobe,'skyglobe.png')

%%
figureBar = figure;
bar(days)

