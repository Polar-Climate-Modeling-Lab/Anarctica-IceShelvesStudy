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
addpath /media/pranab/'Backup Plus'/Backup_12_05_2021/SeaIce_Tony_daily2seasonal/
% addpath /media/pranab/'Backup Plus'/Backup_12_05_2021/SeaIce_Tony_seasonal/
addpath /media/pranab/'Backup Plus'/Backup_12_05_2021/MetUM_ANDREW_Orr_11deg/Daily_temperature_tas/
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

yr1 = 1979; % input ('start year :: ');
yr2 = 2018; % input('end year :: ');

sk = input ('season for sea ice (1 for DJF :: 2 for JJA :: 3 for SON) :: ');
% path_si = "/media/pranab/Backup Plus/Backup_12_05_2021/SeaIce_Tony_daily2seasonal/";
% path_si = "/media/pranab/Backup Plus/Backup_12_05_2021/SeaIce_Tony_seasonal/";

% Loading sea ice data 
% Weddell
fl = dir("/media/pranab/Backup Plus/Backup_12_05_2021/SeaIce_Tony_daily2seasonal/weddell_*.csv");
% fl = dir("/media/pranab/Backup Plus/Backup_12_05_2021/SeaIce_Tony_seasonal/weddell_*.csv");
wdl = [];
a = readtable(fl(sk).name); fl(sk).name
yr = a.Var1; iyr = find (yr>=yr1 & yr<=yr2); % year index to be used for extracting the req year range
if sk == 1
    ds = find(yr == 1987); iyr(ds) = []; % removing the summer of year 1987 since the sea ice data is discontinuous for 1987-88
end
yr = yr(iyr);
% enso = enso(iyr
s1 = a.Var5;
k = 1; 
for i = 1:length(s1)
    wdl(k) = str2double(s1{i}(4:end)); k = k+1;
end
wdl = wdl(iyr);

% ABS
fl = dir("/media/pranab/Backup Plus/Backup_12_05_2021/SeaIce_Tony_daily2seasonal/abs_*.csv");
% fl = dir("/media/pranab/Backup Plus/Backup_12_05_2021/SeaIce_Tony_seasonal/abs_*.csv");

abs = [];
a = readtable(fl(sk).name); fl(sk).name
s1 = a.Var5;
k = 1; 
for i = 1:length(s1)
    abs(k) = str2double(s1{i}(4:end)); k = k+1;
end
abs = abs(iyr);

% Indian
fl = dir("/media/pranab/Backup Plus/Backup_12_05_2021/SeaIce_Tony_daily2seasonal/indian_*.csv");
% fl = dir("/media/pranab/Backup Plus/Backup_12_05_2021/SeaIce_Tony_seasonal/indian_*.csv");

indian = [];
a = readtable(fl(sk).name); fl(sk).name
s1 = a.Var5;
k = 1; 
for i = 1:length(s1)
    indian(k) = str2double(s1{i}(4:end)); k = k+1;
end
indian = indian(iyr);

% Ross
fl = dir("/media/pranab/Backup Plus/Backup_12_05_2021/SeaIce_Tony_daily2seasonal/ross_*.csv");
% fl = dir("/media/pranab/Backup Plus/Backup_12_05_2021/SeaIce_Tony_seasonal/ross_*.csv");

ross = [];
a = readtable(fl(sk).name); fl(sk).name
s1 = a.Var5;
k = 1; 
for i = 1:length(s1)
    ross(k) = str2double(s1{i}(4:end)); k = k+1;
end
ross = ross (iyr);

% West Pacific
fl = dir("/media/pranab/Backup Plus/Backup_12_05_2021/SeaIce_Tony_daily2seasonal/west_pacific_*.csv");
% fl = dir("/media/pranab/Backup Plus/Backup_12_05_2021/SeaIce_Tony_seasonal/west_pacific_*.csv");

wp = [];
a = readtable(fl(sk).name); fl(sk).name
s1 = a.Var5;
k = 1; 
for i = 1:length(s1)
    wp(k) = str2double(s1{i}(4:end)); k = k+1;
end
wp = wp(iyr);

figure(4)
plot(yr, abs, '-sr')
hold on
plot(yr, wdl, '-ok')
plot(yr, indian, '-*b')
plot(yr, ross, '-og')
plot(yr, wp, '-om')
legend('abs','wdl','ind','ross','wp')
clearvars -except wdl abs indian wp ross yr1 yr2 yr iyr
%% load melt potential data
% yr = yr1:yr2;
load HIRHAM_timeseries_corr_pdf.mat; mltH = mlt_p(iyr,:);
load MetUM_timeseries_corr_pdf.mat; mltM = mlt_p(iyr,:);

%%

for jis = 1:length(islf)
    isst = islf{jis};
    IntaryM = mltM(:,jis ); 
    IntaryH = mltH(:,jis );
      
    
    
    %**************************
    f2 = figure(2);
    %**************************
       
    
        % correlation of melt potential (Intary) and sea ice indices
    
        % MetUM
    [wdl1, pwdlM] = corrcoef(wdl,IntaryM); rwdlM = wdl1(1,2);
    [abs1, pabsM] = corrcoef(abs,IntaryM); rabsM = abs1(1,2);
    [indian1, pindianM] = corrcoef(indian,IntaryM); rindianM = indian1(1,2);
    [wp1, pwpM] = corrcoef(wp,IntaryM); rwpM = wp1(1,2);
    [ross1, prossM] = corrcoef(ross,IntaryM); rrossM = ross1(1,2);
    
        % HIRHAM
    [wdl1, pwdlH] = corrcoef(wdl,IntaryH); rwdlH = wdl1(1,2);
    [abs1, pabsH] = corrcoef(abs,IntaryH); rabsH = abs1(1,2);
    [indian1, pindianH] = corrcoef(indian,IntaryH); rindianH = indian1(1,2);
    [wp1, pwpH] = corrcoef(wp,IntaryH); rwpH = wp1(1,2);
    [ross1, prossH] = corrcoef(ross,IntaryH); rrossH = ross1(1,2);
    
    
%     
%     
    rwdl = [rwdlM rwdlH]; pwdl = [pwdlM(1, 2) pwdlH(1, 2)]; % WDL
    rabs = [rabsM rabsH]; pabs = [pabsM(1, 2) pabsH(1, 2)]; % ABS
    rindian = [rindianM rindianH]; pindian = [pindianM(1, 2) pindianH(1, 2)]; % INDIAN
    rwp = [rwpM rwpH]; pwp = [pwpM(1, 2) pwpH(1, 2)]; % WP
    rross = [rrossM rrossH]; pross = [prossM(1, 2) prossH(1, 2)]; % ROSS

    cr = [rwdl;rabs;rindian;rwp;rross];
    
    subplot(4,4,jis)
    hold on
    b = barh(1:5,cr,0.65);
    b(2).FaceColor = 'g'; b(2).EdgeColor = 'none';
    b(1).FaceColor = 'b'; b(1).EdgeColor = 'none';
    
    
    % WDL
    if pwdl(1)<=0.05
        if rwdl(1)>=0
            scatter(0.1+rwdl(1), 0.85, 15, 'ok', 'filled')
        else
            scatter(-0.1+rwdl(1), 0.85, 15, 'ok', 'filled')
        end
    end
    
    if pwdl(2)<=0.05
        if rwdl(2)>=0
            scatter(0.1+rwdl(2), 1.15, 15, 'ok', 'filled')
        else
            scatter(-0.1+rwdl(2), 1.15, 15, 'ok', 'filled')
        end
    end
    

    % ABS
    if pabs(1)<=0.05
        if rabs(1)>=0
            scatter(0.1+rabs(1), 1.85, 15, 'ok', 'filled')
        else
            scatter(-0.1+rabs(1), 1.85, 15, 'ok', 'filled')
        end
    end
    
    if pabs(2)<=0.05
        if rabs(2)>=0
            scatter(0.1+rabs(2), 2.15, 15, 'ok', 'filled')
        else
            scatter(-0.1+rabs(2), 2.15, 15, 'ok', 'filled')
        end
    end
    
    % INDIAN
    if pindian(1)<=0.05
        if rindian(1)>=0
            scatter(0.1+rindian(1), 2.85, 15, 'ok', 'filled')
        else
            scatter(-0.1+rindian(1), 2.85, 15, 'ok', 'filled')
        end
    end
    
    if pindian(2)<=0.05
        if rindian(2)>=0
            scatter(0.1+rindian(2), 3.15, 15, 'ok', 'filled')
        else
            scatter(-0.1+rindian(2), 3.15, 15, 'ok', 'filled')
        end
    end
 
    % WP
    if pwp(1)<=0.05
        if rwp(1)>=0
            scatter(0.1+rwp(1), 3.85, 15, 'ok', 'filled')
        else
            scatter(-0.1+rwp(1), 3.85, 15, 'ok', 'filled')
        end
    end
    
    if pwp(2)<=0.05
        if rwp(2)>=0
            scatter(0.1+rwp(2), 4.15, 15, 'ok', 'filled')
        else
            scatter(-0.1+rwp(2), 4.15, 15, 'ok', 'filled')
        end
    end
    
    % ROSS
    if pross(1)<=0.05
        if rross(1)>=0
            scatter(0.1+rross(1), 4.85, 15, 'ok', 'filled')
        else
            scatter(-0.1+rross(1), 4.85, 15, 'ok', 'filled')
        end
    end
    
    if pross(2)<=0.05
        if rross(2)>=0
            scatter(0.1+rross(2), 5.15, 15, 'ok', 'filled')
        else
            scatter(-0.1+rross(2), 5.15, 15, 'ok', 'filled')
        end
    end
    
    
    xlim([-0.80 0.80])
    ylim([0.5 5.5])
    set(gca, 'fontsize', 6, 'fontweight', 'b')
    title(isst, 'fontweight', 'b')
    if jis > 12
        xlabel('Correlation', 'fontweight', 'b')
    end
    
    if (jis == 1) || (jis == 5) || (jis == 9) || (jis == 13)
        ylabel('Sea ice sectors', 'fontweight', 'b')
    end
    
    yticks(1:5)
    yticklabels({'WDL', 'ABS', 'IND', 'WP', 'ROSS'})
    
    box on
    grid on
   
end

%% SAVE FIGURES
% ['SAVING PLOTS']
% if mx == 1
%     print (f1, '-r600', 'Timeseries_TmeltProb_1979_2018_MetUM_T3hr', '-dpng')
%     print (f2, '-r600', 'PdfDist_TmeltProb_1979_2018_MetUM_T3hr', '-dpng')
%     print (f3, '-r600', 'Corr_indices_TmeltProb_1979_2018_MetUM_T3hr', '-dpng')
% elseif mx == 2
%     print (f1, '-r600', 'Timeseries_TmeltProb_1979_2018_MetUM_Tmax', '-dpng')
%     print (f2, '-r600', 'PdfDist_TmeltProb_1979_2018_MetUM_Tmax', '-dpng')
%     print (f3, '-r600', 'Corr_indices_TmeltProb_1979_2018_MetUM_Tmax', '-dpng')
% end  


