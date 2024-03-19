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
addpath /media/pranab/'Backup Plus'/Backup_12_05_2021/NHM-SMAPv1_1979-1980/

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% lon/lat will be extracted from the sample file, as it will be in degrees
% instead of 'm' as in actual data files

ncfileG = 'NHM-SMAP_v1_1979-1980_3hr.nc';
x = ncread(ncfileG,'LON'); y = ncread(ncfileG,'LAT'); % x & y has same size as the horizontal grid of 'tas'

% tas = ncread(ncfileG,'T2m');
% m_proj('stereographic','lat',-90,'long',0,'radius',31);
% m_pcolor(x, y, squeeze(tas(:,:,1)))
% m_grid('xtick',7,'tickdir','out','ytick',[0 -70 -80],'linest','--', 'fontsize',5,'color',[0.5 0.4 0.7])%,'color','b');
% m_coast('color','r');

% making a vector so as to apply inpolygon
xx = reshape(x,1,[]); glon = xx; % this will correspond to reshaping of tas from model
yy = reshape(y,1,[]); glat = yy;

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


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Extracting temp data for histogram and pdf plots
mx = input (' PLOT PDF for 3-hourly T2m (1) :: PLOT PDF for TMAX (2) :: ');

tasmx = [];
% Intary = []; yr = [];

for jyr = 1979:1980;%2017

    % Dec-centred: 1980 DJF means Dec-1980 to Feb-1981

    % year-1: Dec ----> year-2: Feb
    ncfile = ['NHM-SMAP_v1_tas_3hr_',num2str(jyr),'120103-',num2str(jyr+1),'030100.nc']; ncdisp(ncfile)
    tas = ncread(ncfile,'T2m'); tas = tas + 273.15; % DJF is combined here @ 3-hourly interval

    [m, n, o] = size(tas);
    lnt = 1:o;

    if mx == 2 % for daily maximum
        %%%%% Find Tmax from tas
        tasx = zeros(m, n, 2);
        for it = 1:8:o % 8 because it is a 3-hourly data
            tas1 = max(tas(:,:,it:it+7),[],3);
            tasx = cat(3, tasx, tas1);
        end
        clear tas; tas = tasx(:,:,3:end); clear tasx
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~
    tasmx = cat(3, tasmx, tas);
end

%%
[m, n, o] = size(tasmx);
tas1 = reshape(tasmx,m*n, o);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Intar = zeros(length(tas1),1)+NaN;

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
    x_pdf = [min(tasIS)-10:0.1:max(tasIS)+10];
    y_pdf = pdf(pd,x_pdf);
    
%     hold on
% 	plot(x_pdf, y_pdf)
    
%     lims = (x_pdf >= 271.15) & (x_pdf <= 280.00); % 271.15K is taken because melting starts at -2 degC
    lims = (x_pdf >= 271.15);
    % Integration with trapezoidal method
    Intar(dis(i)) = trapz(x_pdf(lims),y_pdf(lims))*100;

%     Intar(dis(i)) = nanmean(tasIS); % TEST

end


Intar1 = reshape(Intar,[m, n]);
%%
f1 = figure;

pcolorps(y,x,Intar1); caxis([0 100]);
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

f2 = figure;
% mapzoomps('Pine Island Glacier','mapwidth',1000,'inset','se')
mapzoomps('Larsen ice shelf','mapwidth',1000,'inset','se')
% mapzoomps('Dronning Maud Land','mapwidth',3000,'inset','se')

pcolorps(y,x,Intar1); caxis([0 100]);
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
if mx == 1
    print (f1, '-r600', 'SpatialMap_TmeltProb_1979_2018_NHMSMAP_T3hr', '-dpng')
    print (f2, '-r600', 'IsMap_TmeltProb_1979_2018_NHMSMAP_T3hr', '-dpng')
elseif mx == 1
    print (f1, '-r600', 'SpatialMap_TmeltProb_1979_2018_NHMSMAP_Tmax', '-dpng')
    print (f2, '-r600', 'IsMap_TmeltProb_1979_2018_NHMSMAP_Tmax', '-dpng')
end

