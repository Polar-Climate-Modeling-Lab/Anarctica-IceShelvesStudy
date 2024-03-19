% Pranab - 2020

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% takes in iceshelf data
% finds the boundary of the IS
% loads cordex data - transforms rotated to geographic grid
% finds the lon/lat points within the IS boundary
% extracts the temp values for those lon/lat points [and saves in npz
% format]
% DJF season is Dec-centred (Dec/yr1 - Feb/yr2)

% %%%%%%% I M P - NHM-SMAP grid coordinates are in stereo format
%  difficult to work with
%  so extract lon/lat grid from sample file: 'NHM-SMAP_v1_1979-1980_3hr.nc'
%  this is in degree N unit - easier to work with

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
addpath ([aux_pth,'/seqMK/'])
addpath ([aux_pth,'/antbounds_v3.0.2-InSAR/antbounds/'])
addpath ([aux_pth,'/AntarcticMappingTools_v5.16/AntarcticMappingTools/'])
addpath ([aux_pth,'/hatchpattern'])
addpath ([aux_pth,'/m_map/'])
addpath ([aux_pth,'/Polar_Stereo_transform/'])

addpath /home/pranab/Documents/climate_indices_Swetha/
addpath /media/pranab/'Backup Plus'/Backup_12_05_2021/NHM-SMAPv1_1979-1980/
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
yr1 = input ('start year :: ');
yr2 = input('end year :: ');


% ***** CHECK IF INDICES ARE DEC-CENTRED OR NOT ***** %

% Loading Climate indices data
% ASL
load asl_summer_index.mat
yr = asl_summer(:,1);
iyr = find (yr>=yr1 & yr<=yr2);
asl1 = asl_summer(iyr,3); % lon
asl2 = asl_summer(iyr,4); % actual pressure
asl3 = asl_summer(iyr,6); % rel. central pressure

% ENSO
load enso_seasonal_index.mat
yr = enso(:,1);
iyr = find (yr>=yr1 & yr<=yr2);
enso = enso(iyr,2); % summer

% SAM
load sam_seasonal_index.mat
yr = sam_seasonal(:,1);
iyr = find (yr>=yr1 & yr<=yr2);
sam = sam_seasonal(iyr,6); % summer

clearvars -except asl* enso sam yr1 yr2

%%
% Extract the lon/lat from this file as it is in geo degree unit
% (since the complete/actual file sets have stereographic lon/lat, difficuly to work with)
ncfileG = 'NHM-SMAP_v1_1979-1980_3hr.nc';
x = ncread(ncfileG,'LON'); y = ncread(ncfileG,'LAT'); [m,n] = size(x);

% tas = ncread(ncfileG,'T2m');
% m_proj('stereographic','lat',-90,'long',0,'radius',31);
% m_pcolor(x, y, squeeze(tas(:,:,1)))
% m_grid('xtick',7,'tickdir','out','ytick',[0 -70 -80],'linest','--', 'fontsize',5,'color',[0.5 0.4 0.7])%,'color','b');
% m_coast('color','r');

% making a vector so as to apply inpolygon
xx = reshape(x,1,[]); glon = xx;
yy = reshape(y,1,[]); glat = yy;
%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Loops over the different ice shelves to compute PDFs
mx = input ('PLOT PDF for 3-hourly T2m (1) :: PLOT PDF for TMAX (2) :: ');

% List the ice shelves (mentioned in Trusel et al., 2015) + LarsenB &
% Thwaites

islf = {'LarsenC', 'LarsenB', 'Wilkins', 'George VI', 'Venable', 'Abbot', ...
    'Pine Island', 'Thwaites', 'Getz', 'Ross', 'Shackleton', 'West', 'Amery', ...
    'Baudouin', 'Fimbul', 'Ronne'};

for jis = 1:length(islf)
    isst = islf{jis};

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
        [lt,ln] = antbounds_data(isst);
    end

%     plot(ln,lt,'ok')

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

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Extracting temp data for histogram and pdf plots
    tasIS = [];
    Intary = []; yr = [];
    for jyr = yr1:yr2

        % Dec-centred: 1980 DJF means Dec-1980 to Feb-1981

        % year-1: Dec ----> year-2: Feb
        ncfile = ['NHM-SMAP_v1_tas_3hr_',num2str(jyr),'120103-',num2str(jyr+1),'030100.nc']; ncdisp(ncfile)
        tas = ncread(ncfile,'T2m'); tas = tas + 273.15; % DJF is combined here @ 3-hourly interval

        [m, n, o] = size(tas);
        lnt = 1:o;
             
        %%%%% Find Tmax from tas               
        if mx == 2
            tasx = zeros(m, n, 2);
            for it = 1:8:o % 8 --> because of 3-hourly data
                tas1 = max(tas(:,:,it:it+7),[],3);
                tasx = cat(3, tasx, tas1);
            end
            clear tas; tas = tasx(:,:,3:end); clear tasx
        end

        %~~~~~~~~~~~~~~~~~~~~~~~~
        [m, n, o] = size(tas);
        tasyr = [];
        for i = 1: o
            tas1 = reshape(tas(:,:,i),1,[]);
            tasIS = [tasIS; tas1(ds)'];
            tasyr = [tasyr; tas1(ds)'];
        end

        % plotting PDFs of each DJF season
        pd = fitdist(tasyr,'Kernel');

        x_pdf = [min(tasyr)-10:0.1:max(tasyr)+10];
        y = pdf(pd,x_pdf);

        lims = (x_pdf >= 271.15) & (x_pdf <= 280.00);
        Intar1 = trapz(x_pdf(lims),y(lims))*100;
        Intary = [Intary; Intar1];
        yr = [yr; jyr];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % correlation of melt potential (Intary) and climate indices
    
    [ra1, pa1] = corrcoef(asl1,Intary); [ra2, pa2] = corrcoef(asl2,Intary); [ra3, pa3] = corrcoef(asl3,Intary);
    [re, pe] = corrcoef(enso,Intary); [rs, ps] = corrcoef(sam,Intary);
    re = re(1,2); ra1 = ra1(1,2); ra2 = ra2(1,2); ra3 = ra3(1,2); rs = rs(1,2);
    %*********************************************
    f3 = figure(3); % correlation bar chart for the ice shelf
    subplot(4,4,jis)
    %*********************************************
    barh(1, rs, 0.35)
    if ps(1,2)<=0.05
        hold on
        scatter(0, 1, 25, 'ow', 'filled')
        hold on; scatter(0, 1, 25, 'xk')
        hold on; scatter(0, 1, 25, 'ok')
    end
    if rs<=0
        text(0+0.085, 1, 'SAM', 'color', 'b', 'fontsize', 5, 'fontweight', 'b')
    else
        text(rs+0.085, 1, 'SAM', 'color', 'b', 'fontsize', 5, 'fontweight', 'b')
    end
    
    hold on; barh(2,ra1,0.35,'r')
    if pa1(1,2)<=0.05
        scatter(0, 2, 25, 'ow', 'filled')
        hold on; scatter(0, 2, 25, 'xk')
        hold on; scatter(0, 2, 25, 'ok')
    end
    if ra1<=0
        text(0+0.085, 2, 'ASL-lon', 'color', 'r', 'fontsize', 5, 'fontweight', 'b')
    else
        text(ra1+0.085, 2, 'ASL-lon', 'color', 'r', 'fontsize', 5, 'fontweight', 'b')
    end
    
    hold on; barh(3,ra2,0.35,'g')
    if pa2(1,2)<=0.05
        scatter(0, 3, 25, 'ow', 'filled')
        hold on; scatter(0, 3, 25, 'xk')
        hold on; scatter(0, 3, 25, 'ok')
    end
    if ra2<=0
        text(0+0.085, 3, 'ASL-cen', 'color', 'g', 'fontsize', 5, 'fontweight', 'b')
    else
        text(ra2+0.085, 3, 'ASL-cen', 'color', 'g', 'fontsize', 5, 'fontweight', 'b')
    end
    
    hold on; barh(4,ra3,0.35,'m')
    if pa3(1,2)<=0.05
        scatter(0, 4, 25, 'ow', 'filled')
        hold on; scatter(0, 4, 25, 'xk')
        hold on; scatter(0, 4, 25, 'ok')
    end
    if ra3<=0
        text(0+0.085, 4, 'ASL-rel', 'color', 'm', 'fontsize', 5, 'fontweight', 'b')
    else
        text(ra3+0.085, 4, 'ASL-rel', 'color', 'm', 'fontsize', 5, 'fontweight', 'b')
   end
    
    hold on; barh(5,re,0.35,'k')
    if pe(1,2)<=0.05
        scatter(0, 5, 25, 'ow', 'filled')
        hold on; scatter(0, 5, 25, 'xk')
        hold on; scatter(0, 5, 25, 'ok')
    end
    if re<=0
        text(0+0.085, 5, 'ENSO', 'color', 'k', 'fontsize', 5, 'fontweight', 'b')
    else
        text(re+0.085, 5, 'ENSO', 'color', 'k', 'fontsize', 5, 'fontweight', 'b')
    end

    xlim([-1.0 1.0])
    ylim([0.5 5.5])
    set(gca, 'fontsize', 6, 'fontweight', 'b')
    title(isst, 'fontweight', 'b')
    if jis > 12
        xlabel('Correlation', 'fontweight', 'b')
    end
    
    if (jis == 1) || (jis == 5) || (jis == 9) || (jis == 13)
        ylabel('Indices', 'fontweight', 'b')
    end
    yticks(1:5)
    box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting interannual changes in melt potential 
    % (in terms of area under the pdf curve for T>273.15)
    % Fitting trend line to Data
    b = polyfit(yr,Intary, 1);
    fr = polyval(b, yr);
    
    % trend for 1979-1997 & 1998-2018
    if length(yr) >= 38
        ind1 = find(yr == 1997);
        b80 = polyfit(yr(1:ind1),Intary(1:ind1), 1);
        fr80 = polyval(b80, yr(1:ind1));
        b80 = sprintf('%0.2f', b80(1));

        b99 = polyfit(yr(ind1+2:end),Intary(ind1+2:end), 1);
        fr99 = polyval(b99, yr(ind1+1:end));
        b99 = sprintf('%0.2f', b99(1));
    end
    
    % Plot Data & Linear Fit
    %*********************************************
    f1 = figure(1);
    subplot(4,4,jis)
    %*********************************************
%    plot(yr,Intary, 'g', 'linewidth', 2.0)
    plot(yr,Intary, 'color', [0.47 0.67 0.19], 'linewidth', 2.0)
    set(gca, 'fontsize', 6, 'fontweight', 'b')
    title(isst, 'fontweight', 'b')
    %     grid on
%     text(1981, 57.5, num2str(b(1)), 'fontsize', 6, 'fontweight', 'b')
    hold on
    plot(yr, fr, '--k')
    
    if length(yr) >= 38
        plot(yr(ind1+2:end), fr99, 'b', 'linewidth',1.2)
        plot(yr(1:ind1), fr80, 'r', 'linewidth',1.2)
        hold on
        plot((zeros(5,1) + 1999), linspace(0,60,5)', 'k', 'linewidth',1.2)
        text(1981, 56, [b80, ' (yr1-98)'], 'color', 'r', 'fontsize', 5, 'fontweight', 'b')
        text(2001, 56, [b99, ' (99-)'], 'color', 'b', 'fontsize', 5, 'fontweight', 'b')
    end
    if mx == 1
        ylim([0 70]); 
    else
        ylim([0 100]); 
    end
    xlim([yr(1) yr(end)])
    if jis > 12
        xlabel('Years')
    end
    if (jis == 1) || (jis == 5) || (jis == 9) || (jis == 13)
        ylabel('Tmelt Prob.(%)')
    end
    box on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PDF plot for all years DJF
    pd = fitdist(tasIS,'Kernel');
    x_pdf = [min(tasIS)-10:0.1:max(tasIS)+10];
    y = pdf(pd,x_pdf);
    %*********************************************
    f2 = figure(2);
    subplot(4,4,jis)
    %*********************************************
    hold on
    h = histogram(tasIS,'Normalization','pdf');
    % h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'none';% Fit Data
    hold on
    plot(x_pdf,y, 'k')
    
    lims = (x_pdf >= 271.15);% & (x_pdf <= 280.00);
    Intar = trapz(x_pdf(lims),y(lims));

    plot((zeros(5,1) + 271.15), linspace(0,0.5,5)', 'k', 'linewidth',1.2)
    xlim([255.0 280.0])
    ylim([0, 0.45])
    set(gca, 'fontsize', 6, 'fontweight', 'b')
    %title(isst,  'fontsize', 6, 'fontweight', 'b');
    
    text(257,0.42, isst, 'fontsize', 6, 'fontweight', 'b')
%     text(257,0.39, ['Melt Prob: ', num2str(Intar*100), ' %'], 'fontsize', 6, 'fontweight', 'b')
    text(257,0.36, [num2str(Intar*100), '%'], 'fontsize', 6, 'fontweight', 'b')
    if jis > 12
        xlabel('1-hourly tas', 'fontweight', 'b')
    end
    if (jis == 1) || (jis == 5) || (jis == 9) || (jis == 13)
        ylabel('frequency', 'fontweight', 'b')
    end
    box on
 
end
%%
% figure
% seqMK ([yr Intary]);

%     dn = find (in ~= 1);
%     tas1(dn) = nan;
%     
%     
%     %tasP = reshape(tas1,[392,504]);
%     tasP = reshape(tas1,[392,504]);
% 
% % figure
% m_proj('stereographic','lat',-90,'long',0,'radius',31);
% 
% % m_contourf(xx1, yy1, tasP)
% 
% m_contourf(xx1, yy1, squeeze(Dec(:,:,1)'))
% 
% hold on
% 
% m_grid('xtick',7,'tickdir','out','ytick',[0 -70 -80],'linest','--', 'fontsize',5,'color',[0.5 0.4 0.7])%,'color','b');
% 
% 
% m_coast('color','r');
% % 
%  
%%
% SAVE FIGURES
['SAVING PLOTS']
if mx == 1
    print (f1, '-r600', 'Timeseries_TmeltProb_1979_2018_NHMSMAP_T3hr', '-dpng')
    print (f2, '-r600', 'PdfDist_TmeltProb_1979_2018_NHMSMAP_T3hr', '-dpng')
    print (f3, '-r600', 'Corr_indices_TmeltProb_1979_2018_NHMSMAP_T3hr', '-dpng')
elseif mx == 2
    print (f1, '-r600', 'Timeseries_TmeltProb_1979_2018_NHMSMAP_Tmax', '-dpng')
    print (f2, '-r600', 'PdfDist_TmeltProb_1979_2018_NHMSMAP_Tmax', '-dpng')
    print (f3, '-r600', 'Corr_indices_TmeltProb_1979_2018_NHMSMAP_Tmax', '-dpng')
end
