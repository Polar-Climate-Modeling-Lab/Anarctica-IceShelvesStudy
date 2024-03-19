% Pranab - 2020

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% takes in iceshelf data
% finds the boundary of the IS
% loads cordex data - transforms rotated to geographic grid
% finds the lon/lat points within the IS boundary
% extracts the temp values for those lon/lat points [and saves in npz
% format]
% DJF season is Dec-centred (Dec/yr1 - Feb/yr2)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Antarctic boundaries, grounding line, and masks from InSAR
% NEEDS antbound toolbox:
% https://in.mathworks.com/matlabcentral/fileexchange/60246-antarctic-boundaries-grounding-line-and-masks-from-insar?focused=7177816&s_tid=gn_loc_drop&tab=example

% NEEDS Antarctic Mapping Tools:
% https://in.mathworks.com/matlabcentral/fileexchange/47638-antarctic-mapping-tools

% NEEDS Rotated grid transform code: 
% https://in.mathworks.com/matlabcentral/fileexchange/43435-rotated-grid-transform
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear; clc; close all

aux_pth = '/home/pranab/Documents/MATLAB-AUX';

addpath (aux_pth)
addpath ([aux_pth,'/bedmap2_toolbox_v4.6.2/'])
addpath ([aux_pth,'/bedmap2_tiff/'])
addpath ([aux_pth,'/seqMK/'])
addpath ([aux_pth,'/antbounds_v3.0.2-InSAR/antbounds/'])
addpath ([aux_pth,'/AntarcticMappingTools_v5.16/AntarcticMappingTools/'])
addpath ([aux_pth,'/hatchpattern'])
addpath ([aux_pth,'/m_map/'])
addpath ([aux_pth,'/Polar_Stereo_transform/'])
addpath ([aux_pth,'/Scientific_colormaps/'])
addpath ([aux_pth,'/Custom-colormaps/'])

addpath /home/pranab/Documents/climate_indices_Swetha/
addpath /media/pranab/'Backup Plus'/Backup_12_05_2021/MetUM_ANDREW_Orr_11deg/Daily_temperature_tas/
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Rotates the MetUM coordinate and prepares the lon lat vectors
ncfile = 'Antarctic_CORDEX_MetUM_0p11deg_3-hourly_tas_1980.nc' ; % nc file name
% To get information about the nc file
% ncdisp(ncfile)
lat = ncread(ncfile,'grid_latitude');
lon = ncread(ncfile,'grid_longitude');
%contourf(x,y,squeeze(tas(:,:,1)))
[x,y] = meshgrid(lon,lat);

%----------------------rotated grid to geographic grid---------------------
% remapping to regular grid
% BEFORE applying inpolygon
% NEED to convert from rotated grid to geographic grid
xx = reshape(x,1,[]);
yy = reshape(y,1,[]);
grid_in = [xx', yy'];
np_lon = 13.08; np_lat = 0.52; % north pole lon/lat provided in the nc file
SP_coor = [(np_lon - 180) -np_lat]; % SP_coor are the coordinates of the South Pole in the rotated grid
[grid_out] = rotated_grid_transform(grid_in, 2, SP_coor); % transforms from rotated to geographic grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The follwing steps are important because the the pdf calculation will
% loop over the indices of lon/lat points; so the indices of
% glon/glat must match that of extracted tas

xx1 = reshape(grid_out(:,1),[length(lat), length(lon)]); % this reshaping gives the correct lon/lat grid/matrix
yy1 = reshape(grid_out(:,2),[length(lat), length(lon)]);
xx1 = xx1'; yy1 = yy1'; % transpose of xx1 & yy1 has size [504, 392] which is same horizontal grid
% as 'tas'  extracted from model
% Now, reshaping of xx1/yy1 and tas will have all lon lat grid points in
% consistent order
xx3 = reshape(xx1,1,[]); yy3 = reshape(yy1,1,[]); % this will correspond to reshaping of tas from model
% xx3 & yy3 shoud be used for finding lon/lat points within the ice shelves
glon = xx3; glat = yy3;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loops over the different ice shelves to compute PDFs

% List the ice shelves (mentioned in Trusel et al., 2015) + LarsenB &
% Thwaites

load Iceshelf_names_ALL.mat
islf = IS;
% islf = {'LarsenC', 'LarsenB', 'Wilkins', 'George VI', 'Venable', 'Abbot', ...
%     'Pine Island', 'Thwaites', 'Getz', 'Ross', 'Shackleton', 'West', 'Amery', ...
%     'Baudouin', 'Fimbul', 'Ronne'};

% islf = {'LarsenC', 'LarsenB', 'Wilkins', 'George VI', 'Venable', 'Abbot', ...
%     'Pine Island', 'Thwaites', 'Getz', 'Sulzberger', 'Ross', 'Rennick', 'Mertz', 'Moscow University', 'Totten', ...
%     'Shackleton', 'West', 'Amery', 'Baudouin', 'Lazarev', 'Fimbul', 'Filchner', 'Ronne'};

dis = [];
for jis = 1:length(islf)
    isst = islf{jis}

    % finds lon/lat of the boundary of an ice shelf
    % Combined Ross IS
    if strcmp (isst, 'Ross') == 1

        [lt,ln] = antbounds_data('Ross West');
        dp = find(ln<0);
        ln(dp) = ln(dp)+360;

        [lt1,ln1] = antbounds_data('Ross East');
        lnd1 = find(ln1<=0);
        ln1(lnd1) = [];
        lt1(lnd1) = [];

        ln = [ln; ln1];
        lt = [lt; lt1];
    else % for rest of the ice shelves
        [lt,ln] = antbounds_data(isst); % ln, lt --> lon/lat of the boundary of the ice shelf
    end

%     hold on
%     plot(ln,lt,'ok')
%
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % check if a specific point is within the polygon
    %in = inpolygon(-74.6, -101.0, lt,ln);

    glon1 = glon; % this keeps the glon unchanged and any modification 
                  % for Ross (the if loop) will not affect actual glon
    if strcmp (isst, 'Ross') == 1
        dg = find(glon1<0);
        glon1(dg) = glon1(dg)+360;
    end

    in = inpolygon(glat, glon1, lt, ln);
    ds = find (in == 1);
    dis = [dis;ds'];
    length(ds)
    
%     hold on
%     scatter(glon(ds), glat(ds))
end


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Extracting temp data for histogram and pdf plots
mx = input ('PLOT PDF for 3-hourly TAS (1) :: PLOT PDF for TMAX (2) :: ');
tasmx = [];
% Intary = []; yr = [];

for jyr = 1980:1981%2017

    % Dec-centred: 1980 DJF means Dec-1980 to Feb-1981
    % year-1: Dec of jyr
    ncfile = ['Antarctic_CORDEX_MetUM_0p11deg_3-hourly_tas_',num2str(jyr),'.nc'];

    if (rem(jyr, 4)==0) || (rem(jyr, 400)== 0)
        ist = 335*8+1; iend = 31*8;
    else
        ist = 334*8+1; iend = 31*8;
    end

    Dec = ncread(ncfile,'tas', [1 1 ist], [Inf Inf iend]);

    % year-2: JF of jyr+1
    ncfile = ['Antarctic_CORDEX_MetUM_0p11deg_3-hourly_tas_',num2str(jyr+1),'.nc'];

    if (rem(jyr+1, 4)==0) || (rem(jyr+1, 400)== 0)
        ist = 1; iend = (31+29)*8;
    else
        ist = 1; iend = (31+28)*8;
    end

    JF = ncread(ncfile,'tas', [1,1,ist], [Inf, Inf, iend]);

    tas = cat(3, Dec, JF); % 3-hourly DJF tas

    [m, n, o] = size(tas);
    lnt = 1:o;


    %%%%% Find Tmax from tas
    if mx == 2 % for daily maximum
        tasx = zeros(m, n, 2);
        for it = 1:8:o
            tas1 = max(tas(:,:,it:it+7),[],3);
            tasx = cat(3, tasx, tas1);
        end
        clear tas; tas = tasx(:,:,3:end);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~

%     [m, n, o] = size(tas);
    tasmx = cat(3, tasmx, tas);
end

%%

thX1 = input ('what is the melting point threshold? ----- :: (1 for 271.15) :: (2 for 273.15) :: ');
if thX1 == 1
    thX = 271.15;
elseif thX1 == 2
    thX = 273.15;
end

[m, n, o] = size(tasmx);
tas1 = reshape(tasmx,m*n, o);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Intar = zeros(length(tas1),1)+NaN;
prc99_1 = zeros(length(tas1),1)+NaN;

for i = 1:length(dis)
    
    tasIS = tas1(dis(i),:);
    if isnan(sum(tasIS))
        Intar(dis(i)) = NaN;
        continue
    end
    
% % %     h = histogram(tasIS,'Normalization','pdf');
% % %     % h.FaceColor = [0 0.5 0.5];
% % %     h.EdgeColor = 'none';% Fit Data

    % fitted distribution
    pd = fitdist(tasIS','Kernel');
%         pd = fitdist(tasIS,'Normal');
%     x_pdf = [min(tasIS)-10:0.1:max(tasIS)+10];
    x_pdf = [210:0.1:290];
%     x_pdf = [min(tasIS):0.1:max(tasIS)];
%     x_pdf = linspace(min(tasIS), max(tasIS), 300);
    y_pdf = pdf(pd,x_pdf);


% finding 99th percentile for the temp dist - finding the pc99 for x-dist
%     ynz = find(y_pdf>0);
%     xnz = x_pdf(ynz);
%     prc99_1(dis(i)) = prctile (xnz, 95); % 99th percentile of the pdf
    
    prc99_1(dis(i)) = prctile (tasIS', 95); % 99th percentile of the
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     hold on
% 	plot(x_pdf, y_pdf)
    
    lims = (x_pdf >= thX); % hresholds 271.15K and 273.15 are tested

    if (length(find(lims==1))<2)
        Intar(dis(i)) = 0; % When the temp max never crosses melting point
    else
        % Integration with trapezoidal method
        mltAR1 = trapz(x_pdf(lims),y_pdf(lims))*100;
        Intar(dis(i)) = mltAR1;
    end
    

%     if (max(x_pdf)<271.15) || (length(find(lims==1))<2)
%         Intar(dis(i)) = 0; % When the temp max never crosses melting point
%     else
%         % Integration with trapezoidal method
%         Intar(dis(i)) = trapz(x_pdf(lims),y_pdf(lims))*100;
%     end
%     Intar(dis(i)) = nanmean(tasIS); % TEST

end


Intar1 = reshape(Intar,[length(lon), length(lat)]);
prc99 = reshape(prc99_1,[length(lon), length(lat)]);
prc99 = prc99 - thX; prc99(prc99<0) = nan;
%%
f1 = figure

pcolorps(yy1,xx1,Intar1); 
% caxis([0 100]); 
caxis([0 50]);
% caxis([0 45]);

% use the following link to set custom maps
% link: https://in.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap?s_tid=FX_rc2_behav&s_tid=mwa_osa_a
% mycolormap = customcolormap(linspace(0,1,11), {'#a60126','#d7302a','#f36e43','#faac5d','#fedf8d','#fcffbf','#d7f08b','#a5d96b','#68bd60','#1a984e','#006936'});
% colorbar('southoutside');
% colormap(mycolormap);
% axis off;


% #008000 green
% #ADFF2F green yellow
% #CCCC00 strong yellow
% #FF8C00 dark orange
% #FF4500 orangered
% #FF0000 red
% #FF1493 deep pink
% #FF69B4 hot pink

% #800080 purple

% Best so far
mycolormap = customcolormap(linspace(0,1,8), {'#008000', '#ADFF2F',  '#CCCC00', '#FF8C00', '#FF4500', '#ff0000',  '#FF1493', '#FF69B4'});

% mycolormap = customcolormap(linspace(0,1,8), {'#9ACD32', '#ADFF2F',  '#CCCC00', '#FF8C00', '#FF4500', '#ff0000',  '#FF1493', '#FF69B4'});


colorbar('southoutside');
colormap(flipud(mycolormap));
axis off;

% HEATWAVE colormap
% cc = colormap(hsv(30));
% k1 = cc(14,:); k2 = flipud (cc(2:6,:)); k3 = flipud(cc(26:29,:)); k = [k1;k2;k3]; colormap(k)

% test
% cc = colormap(autumn(30)); k = flipud(cc); colormap(k)
% cc = colormap(spring(30)); k = flipud(cc); colormap(k)
% cc = colormap(hot(30)); k = flipud(cc); colormap(k)

% from scientific colormap folder
% load buda.mat;colormap(flipud(buda))

% load lajolla.mat;colormap(lajolla)

% load hawaii.mat; a = hawaii; colormap(flipud(hawaii))


colorbar('fontsize',12, 'fontweight', 'b')
box on

% grounding line data from asaid
load asaid_gl
lat = lat(1:200:end);
lon = lon(1:200:end);
hold on

% patchps(lat,lon,[0.9 0.9 0.9], 'Edgecolor', 'none')
patchps(lat,lon,[0.9 0.9 0.9], 'Edgecolor', 'k')


% Bedmap2 elevation data
[lat,lon,sfz] = bedmap2_data('sfz','res','5 km'); % loads topography from bedmap2
% % [C,h] = contourps(lat,lon,sfz, [0, 200, 500], 'k');
% 
hold on
% [C,h] = contourps(lat,lon,sfz, [200, 200], 'k', 'linewidth', 1.3); % plotting contour (with 'ps' extn): AMT plots
% hold on
% [C,h] = contourps(lat,lon,sfz, [500, 500], 'k', 'linewidth', 1);
% [C,h] = contourps(lat,lon,sfz, [1000, 1000], 'k', 'linewidth', 1);


% [C,h] = contourps(lat,lon,sfz, 1000:1000:5000, 'k', 'linewidth', 1.0);
% clabel(C,h,'LabelSpacing',300,'fontsize',8)
% if mx == 1
%     caxis([0 50])
% end

%% plots for individual ice shelf: AP+Ronne --> Baudouin

f2 = figure
% mapzoomps('Pine Island Glacier','mapwidth',1000,'inset','se')
% mapzoomps('Abbot ice shelf','mapwidth',1000,'inset','se')
% mapzoomps('Dronning Maud Land','mapwidth',3000,'inset','se')

mapzoomps(-77,-21,'size',[4100 2500],'se','frame','off') % AP+Ronne --> Baudouin

pcolorps(yy1,xx1,Intar1); caxis([0 0.05]);caxis([0 100]);
% HEATWAVE colormap
cc = colormap(hsv(30));
k1 = cc(14,:); k2 = flipud (cc(2:6,:)); k3 = flipud(cc(26:29,:)); k = [k1;k2;k3]; colormap(k)
%colorbar('fontweight', 'b')

load asaid_gl
lat = lat(1:200:end);
lon = lon(1:200:end);
hold on

patchps(lat,lon,[0.7 0.7 0.7],'Edgecolor','none')
[lat,lon,sfz] = bedmap2_data('sfz','res','5 km'); % loads topography from bedmap2
% 
% hold on
% [C,h] = contourps(lat,lon,sfz, [200, 200], 'k', 'linewidth', 1.3);
box on
% if mx == 1
%     caxis([0 50])
% end
scalebarps('location','nw')
scarlabel({'Antarctic Peninsula'},'fontangle','bold','fontsize',8, 'edgecolor','none',...
    'fontangle','italic')
scarlabel({'Ronne ice shelf'},'fontangle','italic','fontsize',8)
scarlabel({'Dronning Maud Land'},'fontangle','bold','fontsize',9)


%% plots for individual ice shelf: WA+Ross

f3 = figure

mapzoomps(-81.0,-132,'size',[3000 1600],'ne','frame','off') % WA+Ross

pcolorps(yy1,xx1,Intar1); caxis([0 100]);
% HEATWAVE colormap
cc = colormap(hsv(30));
k1 = cc(14,:); k2 = flipud (cc(2:6,:)); k3 = flipud(cc(26:29,:)); k = [k1;k2;k3]; colormap(k)
%colorbar('fontweight', 'b')

load asaid_gl
lat = lat(1:200:end);
lon = lon(1:200:end);
hold on

patchps(lat,lon,[0.7 0.7 0.7],'Edgecolor','none')
[lat,lon,sfz] = bedmap2_data('sfz','res','5 km'); % loads topography from bedmap2
% 
% hold on
% [C,h] = contourps(lat,lon,sfz, [200, 200], 'k', 'linewidth', 1.3);
box on
% if mx == 1
%     caxis([0 50])
% end
scalebarps
scarlabel({'West Antarctica'},'fontangle','bold','fontsize',10)
scarlabel({'Ross ice shelf'},'fontangle','italic','fontsize',8)
% 'Pine Island Glacier','Getz ice shelf',
%% plots for individual ice shelf: WA+Ross

f4 = figure

mapzoomps(-74.0,108,'size',[2500 3600],'nw','frame','off') % WA+Ross

pcolorps(yy1,xx1,Intar1); caxis([0 100]);
% HEATWAVE colormap
cc = colormap(hsv(30));
k1 = cc(14,:); k2 = flipud (cc(2:6,:)); k3 = flipud(cc(26:29,:)); k = [k1;k2;k3]; colormap(k)
colorbar('fontweight', 'b')

load asaid_gl
lat = lat(1:200:end);
lon = lon(1:200:end);
hold on

patchps(lat,lon,[0.7 0.7 0.7],'Edgecolor','none')
[lat,lon,sfz] = bedmap2_data('sfz','res','5 km'); % loads topography from bedmap2
% 
% hold on
% [C,h] = contourps(lat,lon,sfz, [200, 200], 'k', 'linewidth', 1.3);
box on
% if mx == 1
%     caxis([0 50])
% end

scalebarps('location','se')
scarlabel({'East Antarctica'},'fontangle','bold','fontsize',10)
% scarlabel({'Amery ice shelf'},'fontangle','italic','fontsize',8)

%%
f5 = figure

pcolorps(yy1,xx1,prc99); caxis([0.0 5.0]); %caxis([10.0 15.0]);

% HEATWAVE colormap
cc = colormap(hsv(30));
k1 = cc(14,:); k2 = flipud (cc(2:6,:)); k3 = flipud(cc(26:29,:)); k = [k1;k2;k3]; colormap(k)


colorbar('fontsize',12,'fontweight', 'b')
box on

% grounding line data from asaid
load asaid_gl
lat = lat(1:200:end);
lon = lon(1:200:end);
hold on

patchps(lat,lon,[0.9 0.9 0.9], 'Edgecolor', 'k')


% Bedmap2 elevation data
[lat,lon,sfz] = bedmap2_data('sfz','res','5 km'); % loads topography from bedmap2


%%
% figure
% m_proj('stereographic','lat',-90,'long',0,'radius',31);
% 
% % m_contourf(xx1, yy1, Intar1)
% m_pcolor(xx1, yy1, Intar1)
% 
% % m_contourf(lon, lat, Intar1)
% 
% % m_contourf(xx1, yy1, squeeze(Dec(:,:,1)'))
% 
% hold on
% 
% m_grid('xtick',7,'tickdir','out','ytick',[0 -70 -80],'linest','--', 'fontsize',5,'color',[0.5 0.4 0.7])%,'color','b');
% 
% 
% m_coast('color','r');
% 
% caxis([0 100])
% % 
% %  
% c = colormap(jet(28));
% % if (idat == 3 && idc == 4)
% %     c = colormap(bluewhitered(28));
% % end
% c(11:18,:) = [];
% 
% colormap(c)
%% SAVE FIGURES
% ['SAVING PLOTS']
% if mx == 1
%     print (f1, '-r600', 'SpatialMap_TmeltProb_1979_2018_MetUM_T3hr', '-dpng')
%     print (f2, '-r600', 'IsMap_TmeltProb_1979_2018_MetUM_T3hr', '-dpng')
% elseif mx == 2
%     print (f1, '-r600', 'SpatialMap_TmeltProb_1979_2018_MetUM_Tmax', '-dpng')
%     print (f2, '-r600', 'IsMap_TmeltProb_1979_2018_MetUM_Tmax', '-dpng')
% end
