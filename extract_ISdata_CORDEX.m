% Pranab - 2020

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% takes in iceshelf data
% finds the boundary of the IS
% loads cordex data - transforms rotated to geographic grid
% finds the lon/lat points within the IS boundary
% extracts the temp values for those lon/lat points [and saves in npz
% format]

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Antarctic boundaries, grounding line, and masks from InSAR
% NEEDS antbound toolbox:
% https://in.mathworks.com/matlabcentral/fileexchange/60246-antarctic-boundaries-grounding-line-and-masks-from-insar?focused=7177816&s_tid=gn_loc_drop&tab=example

% NEEDS Antarctic Mapping Tools:
% https://in.mathworks.com/matlabcentral/fileexchange/47638-antarctic-mapping-tools

% NEEDS Rotated grid transform code: 
% https://in.mathworks.com/matlabcentral/fileexchange/43435-rotated-grid-transform
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear; clc;

addpath G:\MATLAB-aux\
addpath G:\MATLAB-aux\antbounds_v3.0.2-InSAR\antbounds\
addpath G:\MATLAB-aux\AntarcticMappingTools_v5.16\AntarcticMappingTools\
addpath I:\MetUM_ANDREW_Orr_11deg\'Daily temperature-tas'\
addpath G:\MY-WORK-LIVE\RealProj_Jclim\codes-data_MATLAB\MATLAB_Aux_files\m_map\

% finds lon/lat of the boundary of an ice shelf
[lt,ln] = antbounds_data('Amery');
plot(ln,lt,'ok')

% check if a specific point is within the polygon
%in = inpolygon(-74.6, -101.0, lt,ln);

ncfile = 'Antarctic_CORDEX_MetUM_0p11deg_3-hourly_tas_1980.nc' ; % nc file name
% To get information about the nc file
ncdisp(ncfile)

lat = ncread(ncfile,'grid_latitude');
lon = ncread(ncfile,'grid_longitude');


%contourf(x,y,squeeze(tas(:,:,1)))

[x,y] = meshgrid(lon,lat);

% remapping to regular grid
% BEFORE applying inpolygon
% NEED to convert from rotated grid to geographic grid

xx = reshape(x,1,[]);
yy = reshape(y,1,[]);
grid_in = [xx', yy'];

np_lon = 13.08; np_lat = 0.52;
SP_coor = [(np_lon - 180) -np_lat];

[grid_out] = rotated_grid_transform(grid_in, 2, SP_coor);
 
xx1 = reshape(grid_out(:,1),[392,504]);
yy1 = reshape(grid_out(:,2),[392,504]);

glon = grid_out(:,1);
glat = grid_out(:,2);

in = inpolygon(glat, glon, lt, ln);
ds = find (in == 1);

% extracting temp data for histogram and pdf plots

for jyr = 1980:2017
%     if jyr == 1
%         yr = 1980;
%     elseif jy == 2
%         yr = 1997;
%     end
%         
    
    ncfile = ['Antarctic_CORDEX_MetUM_0p11deg_3-hourly_tas_',num2str(jyr),'.nc'];

    tas = ncread(ncfile,'tas');

    tasIS = [];
    for i = 1: length(tas)
        tas1 = reshape(tas(:,:,i)',1,[]);
        tasIS = [tasIS; tas1(ds)'];
    end

    dn = find (in ~= 1);
    tas1(dn) = nan;

    %tasP = reshape(tas1,[392,504]);
    tasP = reshape(tas1,[392,504]);

    % Plotting histogram and PDF
    % pd = fitdist(tasIS,'Normal');
    pd = fitdist(tasIS,'Kernel');

    x_pdf = [210:0.1:290];
    y = pdf(pd,x_pdf);

    figure(1)
    
    if jy == 1
        cl = 'b';
    elseif jy == 2
        cl = 'r';
    end
    
    hold on
    h = histogram(tasIS,'Normalization','pdf');
    % h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'none';
    hold on
    plot(x_pdf,y, cl)
    xlim([200.0 300.0])
 

end

%%
figure(2)
m_proj('stereographic','lat',-90,'long',0,'radius',31);

m_contourf(xx1, yy1, tasP)

hold on
% 
% m_contour(lon,lat1,stp1',-100:10:100,'LineColor','k','linewi',1.0)
% 
% m_contour(lon,lat1,stp2',-100:10:100,'LineColor','w','linewi',1.0)
% 
% m_contour(lon,lat1,stp1p',-100:10:100, 'LineColor','k','linewi',1.0)
% % 
% m_contour(lon,lat1,stp1n',-100:10:100,'LineColor','k','linewi',1.0)


%m_contourf(lon,lat1,hgtenso',[200,-100:10:100,200], 'LineColor','none')

m_grid('xtick',7,'tickdir','out','ytick',[0 -70 -80],'linest','--', 'fontsize',5,'color',[0.5 0.4 0.7])%,'color','b');


m_coast('color','r');

 




