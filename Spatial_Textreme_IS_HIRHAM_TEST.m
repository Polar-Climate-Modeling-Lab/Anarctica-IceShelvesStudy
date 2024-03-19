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
addpath ([aux_pth,'/seqMK/'])
addpath ([aux_pth,'/antbounds_v3.0.2-InSAR/antbounds/'])
addpath ([aux_pth,'/AntarcticMappingTools_v5.16/AntarcticMappingTools/'])
addpath ([aux_pth,'/hatchpattern'])
addpath ([aux_pth,'/m_map/'])
addpath ([aux_pth,'/Polar_Stereo_transform/'])

addpath /home/pranab/Documents/climate_indices_Swetha/
addpath /media/pranab/'Backup Plus'/Backup_12_05_2021/HIRHAM_data/
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%~~~~~~~~~~~~~~~~~~~~ T distribution for ice shelves ~~~~~~~~~~~~~~~~~~
% Rotates the MetUM coordinate and prepares the lon lat vectors
ncfile = 'DMI-HIRHAM5_AT6_tas_197912.nc'; % nc file name
% To get information about the nc file
% ncdisp(ncfile)
lat = ncread(ncfile,'rlat');
lon = ncread(ncfile,'rlon');
%contourf(x,y,squeeze(tas(:,:,1)))
[x,y] = meshgrid(lon,lat);

%----------------------rotated grid to geographic grid---------------------
% remapping to regular grid
% BEFORE applying inpolygon
% NEED to convert from rotated grid to geographic grid
xx = reshape(x,1,[]);
yy = reshape(y,1,[]);
grid_in = [xx', yy'];
np_lon = -166.92; np_lat = 6.08; % north pole lon/lat provided in the nc file
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

islf = {'LarsenC', 'LarsenB', 'Wilkins', 'George VI', 'Venable', 'Abbot', ...
    'Pine Island', 'Thwaites', 'Getz', 'Ross', 'Shackleton', 'West', 'Amery', ...
    'Baudouin', 'Fimbul', 'Ronne'};
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

length(dis)
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Extracting temp data for histogram and pdf plots
mx = input ('PLOT PDF for 3-hourly T2m (1) :: PLOT PDF for TMAX (2) :: ');
tasmx = [];
% Intary = []; yr = [];



%%
Intar = zeros(length(glon),1)+NaN;

for i = 1:length(dis)
    i
    tasIS = [];
    
    for jyr = 1979:1980;%2017

        % Dec-centred: 1980 DJF means Dec-1980 to Feb-1981
        % year-1: Dec of jyr
        ncfile = ['DMI-HIRHAM5_AT6_tas_',num2str(jyr),'12.nc'];
        Dec = ncread(ncfile,'tas'); Dec = squeeze(Dec); [m, n, o] = size(Dec); % hourly
        if mx == 1 % this converts hourly data to 3-hourly
            Dec1 = Dec(:,:,1:3:o); clear Dec;
            Dec = Dec1; clear Dec1 % 3-hourly
        end

        % year-2: J& F of jyr+1
        ncfile = ['DMI-HIRHAM5_AT6_tas_',num2str(jyr+1),'01.nc'];
        Jan = ncread(ncfile,'tas'); Jan = squeeze(Jan); [m, n, o] = size(Jan);
        if mx == 1 % this converts hourly data to 3-hourly
            Jan1 = Jan(:,:,1:3:o); clear Jan;
            Jan = Jan1; clear Jan1
        end

        ncfile = ['DMI-HIRHAM5_AT6_tas_',num2str(jyr+1),'02.nc'];
        Feb = ncread(ncfile,'tas'); Feb = squeeze(Feb); [m, n, o] = size(Feb);
        if mx == 1 % this converts hourly data to 3-hourly
            Feb1 = Feb(:,:,1:3:o); clear Feb;
            Feb = Feb1; clear Feb1
        end

        tas = cat(3, Dec, Jan, Feb); % 3-hourly DJF tas

        [m, n, o] = size(tas);
        lnt = 1:o;

        if mx == 2 % for daily maximum
            %%%%% Find Tmax from tas
            tasx = zeros(m, n, 2);
            for it = 1:24:o % skips by 24 since HIRHAM output is hourly
                tas1 = max(tas(:,:,it:it+23),[],3);
                tasx = cat(3, tasx, tas1);
            end
            clear tas; tas = tasx(:,:,3:end); clear tasx
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~
        tasmx = cat(3, tasmx, tas);
    
        [m, n, o] = size(tasmx);
        tas1 = reshape(tasmx,m*n, o);
        tasIS = [tasIS; tas1(dis(i),:)'];
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isnan(sum(tasIS))
        Intar(dis(i)) = NaN;
        continue
    end

    % fitted distribution
    pd = fitdist(tasIS,'Kernel');
%         pd = fitdist(tasIS,'Normal');
    x_pdf = [min(tasIS)-10:0.1:max(tasIS)+10];
    y_pdf = pdf(pd,x_pdf);

    lims = (x_pdf >= 271.15); % 271.15K is taken because melting starts at -2 degC
    % Integration with trapezoidal method
    Intar(dis(i)) = trapz(x_pdf(lims),y_pdf(lims))*100;

%     Intar(dis(i)) = nanmean(tasIS); % TEST

end




Intar1 = reshape(Intar,[length(lon), length(lat)]);
%%
f1 = figure

pcolorps(yy1,xx1,Intar1); caxis([0 100]);
cc = colormap(jet(30));
cc(11:20,:) = [];
colormap(cc)
% caxis([0 100]); cc = colormap(jet(30)); 
% caxis([0 100]); cc = colormap(flipud(hot(30))); 
% c(12:20,:) = [];
% 
% c = colormap(jet(28));
% c(11:18,:) = [];
colorbar('fontweight', 'b')
box on

% grounding line data from asaid
load asaid_gl
lat = lat(1:200:end);
lon = lon(1:200:end);
hold on

patchps(lat,lon,[0.9 0.9 0.9])


% Bedmap2 elevation data
[lat,lon,sfz] = bedmap2_data('sfz','res','5 km'); % loads topography from bedmap2
% % [C,h] = contourps(lat,lon,sfz, [0, 200, 500], 'k');
% 
hold on
[C,h] = contourps(lat,lon,sfz, [200, 200], 'k', 'linewidth', 1.3); % plotting contour (with 'ps' extn): AMT plots
% hold on
% [C,h] = contourps(lat,lon,sfz, [500, 500], 'k', 'linewidth', 1);
% [C,h] = contourps(lat,lon,sfz, [1000, 1000], 'k', 'linewidth', 1);


% [C,h] = contourps(lat,lon,sfz, 1000:1000:5000, 'k', 'linewidth', 1.0);
% clabel(C,h,'LabelSpacing',300,'fontsize',8)
if mx == 1
    caxis([0 50])
end
%% plots for individual ice shlf

f2 = figure
% mapzoomps('Pine Island Glacier','mapwidth',1000,'inset','se')
mapzoomps('Larsen ice shelf','mapwidth',1000,'inset','se')
% mapzoomps('Dronning Maud Land','mapwidth',3000,'inset','se')

pcolorps(yy1,xx1,Intar1); caxis([0 100]);
cc = colormap(jet(30));
cc(11:20,:) = [];
colormap(cc)
% caxis([0 100]); cc = colormap(jet(30)); 
% caxis([0 100]); cc = colormap(flipud(hot(30))); 
% c(12:20,:) = [];
% 
% c = colormap(jet(28));
% c(11:18,:) = [];
colorbar('fontweight', 'b')
load asaid_gl
lat = lat(1:200:end);
lon = lon(1:200:end);
hold on

patchps(lat,lon,[0.9 0.9 0.9])
[lat,lon,sfz] = bedmap2_data('sfz','res','5 km'); % loads topography from bedmap2
% % [C,h] = contourps(lat,lon,sfz, [0, 200, 500], 'k');
% 
hold on
[C,h] = contourps(lat,lon,sfz, [200, 200], 'k', 'linewidth', 1.3);
box on
if mx == 1
    caxis([0 50])
end
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
% if mx == 1
%     print (f1, '-r600', 'SpatialMap_TmeltProb_1979_2018_HIRHAM_T3hr', '-dpng')
%     print (f2, '-r600', 'IsMap_TmeltProb_1979_2018_HIRHAM_T3hr', '-dpng')
% elseif mx == 1
%     print (f1, '-r600', 'SpatialMap_TmeltProb_1979_2018_HIRHAM_Tmax', '-dpng')
%     print (f2, '-r600', 'IsMap_TmeltProb_1979_2018_HIRHAM_Tmax', '-dpng')
% end
% 
% print (f1, '-r600', 'SpatialMap_TmeltProb_1979_2018_HIRHAM', '-dpng')
% print (f2, '-r600', 'IsMap_TmeltProb_1979_2018_HIRHAM', '-dpng')
%  
