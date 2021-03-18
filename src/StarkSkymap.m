close all

fileID = fopen('simbad-trim.csv','r');
C = textscan(fileID,['HIP' '%d' ';' '%f' '%f' '%f' '%f' '%f' '%f']);

starids = zeros(size(C{1}));
starlats = zeros(size(C{1}));
starlons = zeros(size(C{1}));

for i = 1:size(C{1},1)
    
    starids(i) = C{1}(i);
    
    rahr = C{2}(i);
    ramn = C{3}(i);
    rasc = C{4}(i);
    
    rasc = 15*rasc;
    
    carry = floor(rasc/60);
    rasc = rasc - 60*carry;
    
    ramn = 15*ramn + carry;
    
    carry = floor(ramn/60);
    ramn = ramn - 60*carry;
    
    radg = 15*rahr + carry;
    
    starlons(i) = dms2degrees([radg ramn rasc]);
    
    dcdg = C{5}(i);
    dcmn = C{6}(i);
    dcsc = C{7}(i);
    starlats(i) = dms2degrees([dcdg dcmn dcsc]);
end

fileID = fopen('bright_stars_simbad_trim.csv','r');
C = textscan(fileID,['%f' '%f' '%f' '%f' '%f' '%f']);

brightlats = zeros(size(C{1}));
brightlons = zeros(size(C{1}));

for i = 1:size(C{1},1)
    
    rahr = C{1}(i);
    ramn = C{2}(i);
    rasc = C{3}(i);
    
    rasc = 15*rasc;
    
    carry = floor(rasc/60);
    rasc = rasc - 60*carry;
    
    ramn = 15*ramn + carry;
    
    carry = floor(ramn/60);
    ramn = ramn - 60*carry;
    
    radg = 15*rahr + carry;
    
    brightlons(i) = dms2degrees([radg ramn rasc]);
    
    dcdg = C{4}(i);
    dcmn = C{5}(i);
    dcsc = C{6}(i);
    brightlats(i) = dms2degrees([dcdg dcmn dcsc]);
end

% Hubble Deep Field (north), HDF South, HU(X)DF/Chandra South
deeplons = [189.2058,338.2343,53.1625];
deeplats = [62.2161,-60.5507,-27.7914];

[brightlatmat,starlatmat] = ndgrid(brightlats,starlats);
[brightlonmat,starlonmat] = ndgrid(brightlons,starlons);
[arclens,azs] = distance(brightlatmat,brightlonmat,starlatmat,starlonmat);
[closest_to_each_bright,idx_to_each_bright] = min(arclens,[],2);

fileID = fopen('LUVOIR-Architecture_A-NOMINAL_OCCRATES-observations-trim.csv','r');
% HIP,Visit #,Visit dt (years),Exp Time (days)
C = textscan(fileID,'%d %d %f %f','Delimiter',',');

obs_ids = C{1}; % which star is being visited
obs_cts = C{2}; % which visit number to this star
obs_dts = C{3}; % how many years since the first visit to the star
obs_dur = C{4}; % how many days in that observation

%%

% axesm ('globe','Grid', 'on');
% view(60,60)
% axis off
figureMap = figure;
axesm('MapProjection','robinson','Grid','on','GLineWidth',2)
% axesm('MapProjection','stereo','MapLatLimit',[-83 -90],'PLineLocation',1,'ParallelLabel','on','Grid','on','GLineWidth',2)
p1 = scatterm(starlats,starlons,'*', 'linewidth', 2,'DisplayName','Stark 2015 targets');
p2 = scatterm(deeplats,deeplons,'rv', 'linewidth', 2);
p3 = scatterm(brightlats,brightlons,'g+', 'linewidth', 2);
legend([p1 p3 p2],{'Stark 2015 targets','Magnitude 2 stars','Hubble/Chandra deep fields'})
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figureMap,'SkyMap_hires.png')

% figureRA = figure;
% plot(starids,(starlons>-5.183)&(starlons<-4.711),starids,(starlons>-3.404)&(starlons<-3.019),...
%     starids,(starlons>4.515)&(starlons<4.982))


%%
figureMinSeps = figure;
plot(closest_to_each_bright)
title('Closest target star to each bright star')
xlabel('Bright star index value (arb.)')
ylabel('Angular distance to nearest target star (deg.)')