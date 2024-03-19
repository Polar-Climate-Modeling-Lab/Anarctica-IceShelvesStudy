% correlation values between MP-freq and MP-int for MetUM and HIRHAM

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

load MetUM_timeseries_corr_pdf_27315.mat
mltM = fliplr(mlt_p); prcM = fliplr(prc_is) - 271.15;


load HIRHAM_timeseries_corr_pdf_27315.mat
mltH = fliplr(mlt_p); prcH = fliplr(prc_is) - 271.15;

islf = fliplr(islf); islf{1} = 'Ronne-Filchner';
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm = []; rh = [];
cor = cell(16,3);

for jis = 1:16
    
    cor{jis,1} = islf{jis};
    
    a = mltM(:,jis); b = prcM(:,jis); [r, p] = corr(a, b);
    if p <= 0.05
        cor{jis,2} = [num2str(r),'*'];
    else
        cor{jis,2} = num2str(r);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ah = mltH(:,jis); bh = prcH(:,jis); [r, p] = corr(ah, bh);
    if p <= 0.05
        cor{jis,3} = [num2str(r),'*'];
    else
        cor{jis,3} = num2str(r);
    end
   
    cor{jis,1} = islf{jis};
    
end

writecell(cor,'MP_freqVSint_corr.txt','Delimiter','tab');
