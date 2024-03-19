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

thX = input ('melt threshold :: ');

if thX == 271.15
    load MetUM_timeseries_corr_pdf_27115.mat
elseif thX == 273.15
    load MetUM_timeseries_corr_pdf_27315.mat
end

corM = fliplr(cor); pvM = fliplr(pv); mltM = fliplr(mlt_p); 
pdfm = flipud(pdfM); xpdfm = x_pdfM; mltARm = flipud(mltAR); prcM = fliplr(prc_is); prc95M_IS = flipud(prc95_IS)-thX;
load MetUM_timeseries_corr_pdf_3hr.mat
pdfm3 = flipud(pdfM);

if thX == 271.15
    load HIRHAM_timeseries_corr_pdf_27115.mat
elseif thX == 273.15
    load HIRHAM_timeseries_corr_pdf_27315.mat
end
corH = fliplr(cor); pvH = fliplr(pv); mltH = fliplr(mlt_p); 
pdfh = flipud(pdfM); xpdfh = x_pdfM; mltARh = flipud(mltAR); prcH = fliplr(prc_is); prc95H_IS = flipud(prc95_IS)- thX;
load HIRHAM_timeseries_corr_pdf_3hr.mat
pdfh3 = flipud(pdfM);

% 
islf = fliplr(islf); islf{1} = 'Ronne-Filchner';
%%
tb = cell(16,3); tb2 = cell(16,3); 
for jis = 1:length(islf)
    
    isst = islf{jis};
 
   % Plotting interannual changes in melt potential 
    % (in terms of area under the pdf curve for T>thX)
    %***********************************
    f1 = figure(1);
    %***********************************
    
    yrM = [];
    
    IntaryM = mltM(:,jis );     
    IntaryM1 = IntaryM; % for plotting the unfiltered timeseries later
    [b,filb,fila]=lopass_butterworth(IntaryM,1/20,1,4); IntaryM = b; 
    tipM = seqMK([yr' IntaryM]); itM = find(tipM(:,3) == 1);
    
    % Trend of lp-filtered time series (full length) - MetUM
    [trM,intM] = trend(IntaryM); [hM,p_vM] = Mann_Kendall(IntaryM',0.05); % h1==1 ---> there is trend present
    % Fitting trend line to Data
    bM1 = polyfit(yr,IntaryM', 1); fr1M = polyval(bM1, yr); trM = sprintf('%0.2f', trM(1));
    
    if ~isempty(itM)
        yrM = yr(itM);
    end
    
     % MetUM
    subplot(4,4,jis)
    
    plot(yr,IntaryM1, 'color', 'b', 'linewidth', 2.0)
    hold on
    
    cr3M = 0;
    
	if ((jis >= 14) && (jis <= 16))||((jis >= 7) && (jis <= 10))||((jis >= 2) && (jis <= 3))  %%% picking the only ice shelves that show trend reversal
        
        if ~isempty(yrM)
            for j = 1:length(yrM)
                
                ind1 = find (yr == yrM(j));
                
                if (length(yr) - ind1) <= 1
                    continue
                end
                
                yrd1 = yrM(j) - yr1; yrd2 = yr2- yrM(j);

                [tr1,int1] = trend(IntaryM(1:ind1-1)); [H1,p_value1] = Mann_Kendall(IntaryM(1:ind1)',0.05); % H1==1 ---> there is trend present
                [tr2,int2] = trend(IntaryM(ind1+1:end)); [H2,p_value2] = Mann_Kendall(IntaryM(1:ind1)',0.05);

                if (H1 == 1) && (H2 == 1)
                    
%                     yrd1, yrd2
%                     scatter(yrM(j),0,'ob','filled') % plotting all trend reversal points
                    
                    if (yrd1 >= 15) && (yrd2 >= 15)
                        
                        scatter(yrM(j),0,'ob','filled')
                        
                        cr3M = 1;

                        % Fitting trend line to Data
                        b1 = polyfit(yr(1:ind1-1),IntaryM(1:ind1-1)', 1); fr1 = polyval(b1, yr(1:ind1-1)); b1 = sprintf('%0.2f', b1(1));
                        b2 = polyfit(yr(ind1+1:end),IntaryM(ind1+1:end)', 1); fr2 = polyval(b2, yr(ind1+1:end)); b2 = sprintf('%0.2f', b2(1));  

                        % trend line plots
                        hold on
                        
                        if (jis == 2 && yrM(j) == 1999) || (jis == 3 && yrM(j) == 1995) % Only for the CP with max trend difference
                            plot(yr(ind1+1:end), fr2, 'k', 'linewidth',1.5); fr2MU = fr2; yr2MU = yr(ind1+1:end);
                            plot(yr(1:ind1-1), fr1, 'k', 'linewidth',1.5); fr1MU = fr1; yr1MU = yr(1:ind1-1);
    %                         whos fr1 fr2
                        end
 
                        [isst,'--','MetUM','---',num2str(yrM(j)),'::',num2str(tr1),'::',num2str(tr2)]
%                         [yrM(j) tr1 tr2]
                    end
                      
                end
             
            end
        end

    end
    
%     'cr3M::', cr3M
    
	% Plotting the whole trend only when there is no shift (satisfying all 3 criteria)

    if cr3M == 0
        if hM == 1
            plot(yr, fr1M, 'k', 'linewidth',1.5)
    %             text(1981, 90, [num2str(trM,'%.2f')], 'color', 'b', 'fontsize', 7, 'fontweight', 'b')
        else
            plot(yr, fr1M, '--k', 'linewidth',0.6)
        end
    end
    
    hold on
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HIRHAM
    yrH = [];
     
    IntaryH = mltH(:,jis );
    IntaryH1 = IntaryH; % for plotting the unfiltered timeseries later
    [b,filtb,filta]=lopass_butterworth(IntaryH,1/20,1,4); IntaryH = b;  
    tipH = seqMK([yr' IntaryH]); itH = find(tipH(:,3) == 1);
    
    if ~isempty(itH)
        yrH = yr(itH);
    end
    
    % HIRHAM
    subplot(4,4,jis)
    hold on

    % Trend of lp-filtered time series (full length) - HIRHAM
    [trH,intH] = trend(IntaryH); [hH,p_vH] = Mann_Kendall(IntaryH',0.05); % h1 ---> there is trend present
    % Fitting trend line to Data
    bH1 = polyfit(yr,IntaryH', 1); fr1H = polyval(bH1, yr); trH = sprintf('%0.2f', trH(1));
    
    plot(yr,IntaryH1, 'color', 'g', 'linewidth', 2.0); grid on
    hold on
    cr3 = 0;
    

    if ((jis >= 14) && (jis <= 16))||((jis >= 7) && (jis <= 10))||((jis >= 2) && (jis <= 3))
        
        if ~isempty(yrH)
            
            
            for j = 1:length(yrH)
                
                ind1 = find (yr == yrH(j));
                
%                 if (length(yr) - ind1) <= 1
%                     continue
%                 end
                
                yrd1 = yrH(j) - yr1; yrd2 = yr2- yrH(j);

                [tr1,int1] = trend(IntaryH(1:ind1-1)); [H1,p_value1] = Mann_Kendall(IntaryH(1:ind1)',0.05); % H1 ---> there is trend present
                [tr2,int2] = trend(IntaryH(ind1+1:end)); [H2,p_value2] = Mann_Kendall(IntaryH(1:ind1)',0.05);

                if (H1 == 1) && (H2 == 1)
%                     yrd1, yrd2
%                     scatter(yrH(j),0,'og','filled') % plotting all trend reversal points
                    
                    if (yrd1 >= 15) && (yrd2 >= 15)
                        cr3 = 1;
                        scatter(yrH(j),0,'og','filled')

                        % Fitting trend line to Data
                        b1 = polyfit(yr(1:ind1-1),IntaryH(1:ind1-1)', 1); fr1 = polyval(b1, yr(1:ind1-1)); b1 = sprintf('%0.2f', b1(1));
                        b2 = polyfit(yr(ind1+1:end),IntaryH(ind1+1:end)', 1); fr2 = polyval(b2, yr(ind1+1:end)); b2 = sprintf('%0.2f', b2(1));  

                        % trend line plots
                        hold on
                        if (jis == 2 && yrH(j) == 2002) || (jis == 3 && yrH(j) == 2003)
                            plot(yr(ind1+1:end), fr2, 'r', 'linewidth',1.5)
                            plot(yr(1:ind1-1), fr1, 'r', 'linewidth',1.5)
                        end
                        
                        [isst,'--','HIRHAM','---',num2str(yrH(j)),'::',num2str(tr1),'::',num2str(tr2)]
%                         [yrH(j) tr1 tr2]
                    end
                    
                end
            end
            
            
        end
    end

    % Plotting the whole trend only when there is no shift (satisfying all 3 criteria)
    if cr3 == 0
        if hH == 1
            plot(yr, fr1H, 'r', 'linewidth',1.5)
        else
            plot(yr, fr1H, '--r', 'linewidth',0.6)
        end
    end

    set(gca, 'fontsize', 6, 'fontweight', 'b')
    title(isst, 'fontweight', 'b')
%     ylim([0 100])
    ylim([0 60])

    xlim([yr(1) yr(end)])
    if jis > 12
        xlabel('Years')
    end
    
    
%     text(1981, 56, [b80, ' (80-98)'], 'color', 'r', 'fontsize', 5, 'fontweight', 'b')
%     text(2001, 56, [b99, ' (99-)'], 'color', 'b', 'fontsize', 5, 'fontweight', 'b')

    
    if (jis == 1) || (jis == 5) || (jis == 9) || (jis == 13)
        ylabel('MPI - freq (%)')
    end
    box on
%     grid on

    if jis == 1
        text(1983, 70, ['MetUM'], 'color', 'b', 'fontsize', 10, 'fontweight', 'b')
        text(2010, 70, ['HIRHAM5'], 'color', 'g', 'fontsize', 10, 'fontweight', 'b')
    end


%%%%% JUST to overlay metUM trend

    if cr3M == 0
        if hM == 1
            plot(yr, fr1M, 'k', 'linewidth',1.5)
        else
            plot(yr, fr1M, '--k', 'linewidth',0.6)
        end
    elseif cr3M == 1
        plot(yr2MU, fr2MU, 'k', 'linewidth',1.5)
        plot(yr1MU, fr1MU, 'k', 'linewidth',1.5)        
    end
    
    
    % writing trend values to a cell
    tb{jis,1} = islf{jis}; 
    if hM == 1
        tb{jis,2} = [num2str(trM),'*']; 
    else
        tb{jis,2} = num2str(trM); 
    end
    if hH == 1
        tb{jis,3} = [num2str(trH),'*']; 
    else
        tb{jis,3} = num2str(trH);
    end
    %**************************
    f2 = figure(2);
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
    f3 = figure(3);
    %***********************************
    subplot(4,4,jis)
    plot(xpdfm, pdfm(jis,:), 'color', 'b', 'linewidth', 1.2)
    hold on
    ynz = find(pdfm(jis,:)>0); xnz = xpdfm(ynz); prc95 = prctile (xnz, 95); prcm95 = prc95 - thX;
    plot(xpdfm, pdfm3(jis,:), 'color', 'b', 'linewidth', 1.2, 'linestyle', '--')
    
    plot(xpdfh, pdfh(jis,:), 'color', 'g', 'linewidth', 1.2)
    ynz2 = find(pdfh(jis,:)>0); xnz2 = xpdfm(ynz2); prc95 = prctile (xnz2, 95); prch95 = prc95 - thX;
    plot(xpdfm, pdfh3(jis,:), 'color', 'g', 'linewidth', 1.2, 'linestyle', '--')
    
    
    plot((zeros(5,1) + thX), linspace(0,1.0,5)', 'k', 'linewidth',1.2)
    
    plot((zeros(5,1) + 271.15), linspace(0,1.0,5)', '--k', 'linewidth',0.7)
    
    text(248, 0.45, [num2str(mltARm(jis),'%.2f') '%'], 'color', 'b', 'fontsize', 6, 'fontweight', 'b')
    text(248, 0.25, [num2str(mltARh(jis),'%.2f') '%'], 'color', 'g', 'fontsize', 6, 'fontweight', 'b')
    
    text(248, 0.38, ['(',num2str(prc95M_IS(jis),'%.2f'), 'K',')'], 'color', 'b', 'fontsize', 5, 'fontweight', 'b')
    text(248, 0.18, ['(',num2str(prc95H_IS(jis),'%.2f'), 'K',')'], 'color', 'g', 'fontsize', 5, 'fontweight', 'b')
    
%     
%     text(248, 0.38, ['(',num2str(prcm95,'%.2f'), 'K',')'], 'color', 'b', 'fontsize', 5, 'fontweight', 'b')
%     text(248, 0.18, ['(',num2str(prch95,'%.2f'), 'K',')'], 'color', 'g', 'fontsize', 5, 'fontweight', 'b')
    
    xlim([245 280])
%     ylim([0.0 0.525])
    ylim([0.0 0.225])
    
    set(gca, 'fontsize', 6, 'fontweight', 'b')
    title(isst, 'fontweight', 'b')
    if jis > 12
        xlabel('Temperature', 'fontweight', 'b')
    end
    
    if (jis == 1) || (jis == 5) || (jis == 9) || (jis == 13)
        ylabel('frequency', 'fontweight', 'b')
    end
%     yticks(1:4)
    ylim([0.0 0.6])
    box on
    
    if jis == 1
        text(245, 0.775, ['MetUM'], 'color', 'b', 'fontsize', 10, 'fontweight', 'b')
        text(280, 0.775, ['HIRHAM5'], 'color', 'g', 'fontsize', 10, 'fontweight', 'b')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        
%   % Plotting interannual changes in 99th percentile 
%   %***********************************
    f4 = figure(4);
    %***********************************

     yrM = [];
    
    prc99M = prcM(:,jis);     
    prc99M = prc99M - thX;
     
    % Trend of lp-filtered time series (full length) - MetUM
    [trM,intM] = trend(prc99M); [hM,p_vM] = Mann_Kendall(prc99M',0.05); % h1==1 ---> there is trend present
    % Fitting trend line to Data
    bM1 = polyfit(yr,prc99M', 1); fr1M = polyval(bM1, yr); %trM = sprintf('%0.2f', trM(1));
    
    if ~isempty(itM)
        yrM = yr(itM);
    end

%     if ((jis ==1) || (jis == 8))
%         prcl1 = -5; prcl2 = 1;
%     else
%         prcl1 = -2; prcl2 = 4;
%     end

    prcl1 = -4.25; prcl2 = 3.1;
    
    
    % MetUM
    subplot(4,4,jis)

    plot(yr,prc99M, 'color', 'b', 'linewidth', 2.0)
%     bar(yr,prc99M1-273.15); %, 'color', 'b', 'linewidth', 2.0)

    hold on
    if hM == 1
        plot(yr, fr1M, 'k', 'linewidth',1.5)
    else
        plot(yr, fr1M, '--k', 'linewidth',0.6)
%         if ((jis >= 1 && jis <= 4))
%         text(2003, tloc, [num2str(trM,'%.3f')], 'color', 'b', 'fontsize', 7, 'fontweight', 'b')
    end
  
    
    set(gca, 'fontsize', 6, 'fontweight', 'b')
    title(isst, 'fontweight', 'b')

    xlim([yr(1) yr(end)])
    if jis > 12
        xlabel('Years')
    end
  
    % HIRHAM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    yrH = [];
    prc99H = prcH(:,jis);     
    prc99H = prc99H - thX;

    
    % Trend of lp-filtered time series (full length) - MetUM
    [trH,intH] = trend(prc99H); [hH,p_vH] = Mann_Kendall(prc99H',0.05); % h1==1 ---> there is trend present
    % Fitting trend line to Data
    bH1 = polyfit(yr,prc99H', 1); fr1H = polyval(bH1, yr); %trM = sprintf('%0.2f', trM(1));

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(yr,prc99H, 'color', 'g', 'linewidth', 2.0)
%     bar(yr,prc99H1-273.15); %, 'color', 'b', 'linewidth', 2.0)

    hold on
    if hH == 1
        plot(yr, fr1H, 'r', 'linewidth',1.5)
%         if ((jis >= 1 && jis <= 4))
%         text(1985, tloc, [num2str(trH,'%.3f')], 'color', 'g', 'fontsize', 7, 'fontweight', 'b')
%         else
%             text(1981, 15, [num2str(trM,'%.2f')], 'color', 'b', 'fontsize', 7, 'fontweight', 'b')
%         end
    else
        plot(yr, fr1H, '--r', 'linewidth',0.6)
    end
  
    
    set(gca, 'fontsize', 6, 'fontweight', 'b')
    title(isst, 'fontweight', 'b')

    xlim([yr(1) yr(end)])
    if jis > 12
        xlabel('Years')
    end
      
    if (jis == 1) || (jis == 5) || (jis == 9) || (jis == 13)
        ylabel('MPI - int (K)')
    end
    
    %%%%% this is just to overlay MetUM trend line on top
    if hM == 1
        plot(yr, fr1M, 'k', 'linewidth',1.5)
    else
        plot(yr, fr1M, '--k', 'linewidth',0.6)
%         if ((jis >= 1 && jis <= 4))
%         text(2003, tloc, [num2str(trM,'%.3f')], 'color', 'b', 'fontsize', 7, 'fontweight', 'b')
    end
    grid on
    hold on
    
    ylim([prcl1 prcl2])

    if jis == 1
        text(1983, 2.5, ['MetUM'], 'color', 'b', 'fontsize', 10, 'fontweight', 'b')
        text(2010, 2.5, ['HIRHAM5'], 'color', 'g', 'fontsize', 10, 'fontweight', 'b')
    end
    
    tb2{jis,1} = islf{jis}; 
    if hM == 1
        tb2{jis,2} = [num2str(trM),'*']; 
    else
        tb2{jis,2} = num2str(trM); 
    end
    if hH == 1
        tb2{jis,3} = [num2str(trH),'*']; 
    else
        tb2{jis,3} = num2str(trH);
    end
end

if thX == 273.15
    writecell(tb,'MP_freq_Trend_27315.txt','Delimiter','tab');
    writecell(tb2,'MP_Int_Trend_27315.txt','Delimiter','tab');
elseif thX == 271.15
    writecell(tb,'MP_freq_Trend_27115.txt','Delimiter','tab');
    writecell(tb2,'MP_Int_Trend_27115.txt','Delimiter','tab');
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



