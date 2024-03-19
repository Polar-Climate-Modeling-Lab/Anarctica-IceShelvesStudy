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
addpath ([aux_pth,'/Butterworth Filters/'])

addpath /home/pranab/Documents/climate_indices_Swetha/
addpath /media/pranab/'Backup Plus'/Backup_12_05_2021/MetUM_ANDREW_Orr_11deg/Daily_temperature_tas/
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

yr1 = 1979; %input ('start year :: ');
yr2 = 2018; % input('end year :: ');
yr = yr1:yr2;

% load HIRHAM_timeseries_corr.mat;
% corH = cor; pvH = pv; mltH = mlt_p;

% load MetUM_timeseries_corr.mat; % variables loaded: cor, pv, mlt_p, islf
% corM = cor; pvM = pv; mltM = mlt_p; 

load MetUM_timeseries_corr_pdf.mat
corM = cor; pvM = pv; mltM = mlt_p; pdfm = pdfM; xpdfm = x_pdfM; mltARm = mltAR;

load HIRHAM_timeseries_corr_pdf.mat
corH = cor; pvH = pv; mltH = mlt_p; pdfh = pdfM; xpdfh = x_pdfM; mltARh = mltAR;

%%
for jis = 1:1%length(islf)
    
    isst = islf{jis};
    IntaryM = mltM(:,jis ); 
    IntaryH = mltH(:,jis );

    
    yrM = []; yrH = [];
    tipM = seqMK([yr' smooth(IntaryM)]); itM = find(tipM(:,3) == 1);
    tipH = seqMK([yr' smooth(IntaryH)]); itH = find(tipH(:,3) == 1);
    
%     tipM = seqMK([yr' IntaryM]); itM = find(tipM(:,3) == 1);
%     tipH = seqMK([yr' IntaryH]); itH = find(tipH(:,3) == 1);
    if ~isempty(itM)
        yrM = yr(itM);
    end
    if ~isempty(itH)
        yrH = yr(itH);
    end
    
    % Plotting interannual changes in melt potential 
    % (in terms of area under the pdf curve for T>271.15)
    %***********************************
    f1 = figure(1)
    %***********************************
    % Fitting trend line to Data
     ind1 = find(yr == 1997);
    
    % MetUM
    % Fitting trend line to Data
 
    b = polyfit(yr,IntaryM', 1); fr = polyval(b, yr); 
    b80 = polyfit(yr(1:ind1),IntaryM(1:ind1)', 1); fr80 = polyval(b80, yr(1:ind1)); b80 = sprintf('%0.2f', b80(1));
    b99 = polyfit(yr(ind1+2:end),IntaryM(ind1+2:end)', 1); fr99 = polyval(b99, yr(ind1+2:end)); b99 = sprintf('%0.2f', b99(1));   
    
    % MetUM
    subplot(4,4,jis)
    plot(yr,IntaryM, 'color', 'b', 'linewidth', 2.0)
    hold on
%     plot(yr,IntaryM, 'color', [0.47 0.67 0.19], 'linewidth', 2.0)

	if ((jis >= 1) & (jis <= 3))|((jis >= 7) & (jis <= 10))|((jis >= 14) & (jis <= 15))  %%% picking the only ice shelves that show trend reversal
        if ~isempty(yrM)
            for j = 1:length(yrM)
%                 plot((zeros(5,1) + yrM(1)), linspace(0,100,5)', 'b', 'linewidth',1.2)
%                 if (yrM(j) >= 1994) & (yrM(j) <= 2002)
                if (yrM(j) >= 1980) & (yrM(j) <= 2015)
                    plot((zeros(5,1) + yrM(j)), linspace(0,100,5)', 'b', 'linewidth',1.2)
                end
            end
        end
    end
    
    hold on
    
    % trend line plots
%     plot(yr, fr, '--k')
%     plot(yr(ind1+2:end), fr99, '--k', 'linewidth',1.2)
%     plot(yr(1:ind1), fr80, 'k', 'linewidth',1.2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HIRHAM
    % Fitting trend line to Data
    b = polyfit(yr,IntaryH', 1); fr = polyval(b, yr);
    b80 = polyfit(yr(1:ind1),IntaryH(1:ind1)', 1); fr80 = polyval(b80, yr(1:ind1)); b80 = sprintf('%0.2f', b80(1));
    b99 = polyfit(yr(ind1+2:end),IntaryH(ind1+2:end)', 1); fr99 = polyval(b99, yr(ind1+2:end)); b99 = sprintf('%0.2f', b99(1));

    % HIRHAM
    subplot(4,4,jis)
    hold on
    plot(yr,IntaryH, 'color', 'g', 'linewidth', 2.0)
    hold on
    
    if ((jis >= 1) & (jis <= 3))|((jis >= 7) & (jis <= 10))|((jis >= 14) & (jis <= 15))
        if ~isempty(yrH)
            for j = 1:length(yrH)
%                 plot((zeros(5,1) + yrM(1)), linspace(0,100,5)', 'g', 'linewidth',1.2)
%                 if (yrH(j) >= 1994) & (yrH(j) <= 2002)
                if (yrH(j) >= 1980) & (yrH(j) <= 2015)
                    plot((zeros(5,1) + yrH(j)), linspace(0,100,5)', 'g', 'linewidth',1.2)
                end
            end
        end
    end



%     if (jis >= 7) & (jis <= 10)
%         plot((zeros(5,1) + 1998), linspace(0,100,5)', '--k', 'linewidth',1.2)
%     end
        % trend line plots
%     plot(yr(ind1+2:end), fr99, '--k', 'linewidth',1.2)
%     plot(yr(1:ind1), fr80, 'k', 'linewidth',1.2)

    set(gca, 'fontsize', 6, 'fontweight', 'b')
    title(isst, 'fontweight', 'b')
    ylim([0 100])

    xlim([yr(1) yr(end)])
    if jis > 12
        xlabel('Years')
    end
    
    
%     text(1981, 56, [b80, ' (80-98)'], 'color', 'r', 'fontsize', 5, 'fontweight', 'b')
%     text(2001, 56, [b99, ' (99-)'], 'color', 'b', 'fontsize', 5, 'fontweight', 'b')

    
    if (jis == 1) || (jis == 5) || (jis == 9) || (jis == 13)
        ylabel('Melt Pot.(%)')
    end
    box on
%     grid on
    
    
    %**************************
    f2 = figure(2)
    %**************************
       
    re = [corM(1, jis) corH(1, jis)]; pe = [pvM(1, jis) pvH(1, jis)]; % ENSO
    ra1 = [corM(2, jis) corH(2, jis)]; pa1 = [pvM(2, jis) pvH(2, jis)]; % ASL-lon
    ra3 = [corM(3, jis) corH(3, jis)]; pa3 = [pvM(3, jis) pvH(3, jis)]; % ASL-rel
    rs = [corM(4, jis) corH(4, jis)]; ps = [pvM(4, jis) pvH(4, jis)]; % SAM
    
    cr = [rs;ra1;ra3;re];
    
    subplot(4,4,jis)
    hold on
    b = barh(1:4,cr,0.65);
    b(2).FaceColor = 'g'; b(2).EdgeColor = 'none';
    b(1).FaceColor = 'b'; b(1).EdgeColor = 'none';
    % SAM
   
    if ps(1)<=0.05
        if rs(1)>=0
            scatter(0.1+rs(1), 0.85, 15, 'ok', 'filled')
        else
            scatter(-0.1+rs(1), 0.85, 15, 'ok', 'filled')
        end
    end
    
    if ps(2)<=0.05
        if rs(2)>=0
            scatter(0.1+rs(2), 1.15, 15, 'ok', 'filled')
        else
            scatter(-0.1+rs(2), 1.15, 15, 'ok', 'filled')
        end
    end
    

    % asl-rel

    if pa1(1)<=0.05
        if ra1(1)>=0
            scatter(0.1+ra1(1), 1.85, 15, 'ok', 'filled')
        else
            scatter(-0.1+ra1(1), 1.85, 15, 'ok', 'filled')
        end
    end
    
    if pa1(2)<=0.05
        if ra1(2)>=0
            scatter(0.1+ra1(2), 2.15, 15, 'ok', 'filled')
        else
            scatter(-0.1+ra1(2), 2.15, 15, 'ok', 'filled')
        end
    end
    
    % ASL-lon
    if pa3(1)<=0.05
        if ra3(1)>=0
            scatter(0.1+ra3(1), 2.85, 15, 'ok', 'filled')
        else
            scatter(-0.1+ra3(1), 2.85, 15, 'ok', 'filled')
        end
    end
    
    if pa3(2)<=0.05
        if ra3(2)>=0
            scatter(0.1+ra3(2), 3.15, 15, 'ok', 'filled')
        else
            scatter(-0.1+ra3(2), 3.15, 15, 'ok', 'filled')
        end
    end
    
    xlim([-0.80 0.80])
    ylim([0.5 4.5])
    set(gca, 'fontsize', 6, 'fontweight', 'b')
    title(isst, 'fontweight', 'b')
    if jis > 12
        xlabel('Correlation', 'fontweight', 'b')
    end
    
    if (jis == 1) || (jis == 5) || (jis == 9) || (jis == 13)
        ylabel('Indices', 'fontweight', 'b')
    end
    yticks(1:4)
    box on
    grid on
    
    
    %***********************************
    f3 = figure(3)
    %***********************************
    subplot(4,4,jis)
    plot(xpdfm, pdfm(jis,:), 'color', 'b', 'linewidth', 2.0)
    hold on
    plot(xpdfh, pdfh(jis,:), 'color', 'g', 'linewidth', 2.0)
    plot((zeros(5,1) + 271.15), linspace(0,1.0,5)', 'k', 'linewidth',1.2)
    
    text(250, 0.35, [num2str(mltARm(jis),'%.2f') '%'], 'color', 'b', 'fontsize', 7, 'fontweight', 'b')
    text(250, 0.25, [num2str(mltARh(jis),'%.2f') '%'], 'color', 'g', 'fontsize', 7, 'fontweight', 'b')
    
    xlim([245 280])
    ylim([0.0 0.525])
    
        set(gca, 'fontsize', 6, 'fontweight', 'b')
    title(isst, 'fontweight', 'b')
    if jis > 12
        xlabel('Daily t-max', 'fontweight', 'b')
    end
    
    if (jis == 1) || (jis == 5) || (jis == 9) || (jis == 13)
        ylabel('frequency', 'fontweight', 'b')
    end
    yticks(1:4)
    box on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

%% SAVE FIGURES
% ['SAVING PLOTS']

% print (f1, '-r600', 'Timeseries_MeltP_Tmax', '-dpng')
% print (f2, '-r600', 'Corr_indices_TmeltProb_Tmax', '-dpng')



