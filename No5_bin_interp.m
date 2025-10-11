%% No5_bin_interp.m
% NSE
%Removal of noisy ensembles and creates bin interpolated matrices (25m)

% Applies draft offset to depths.
% Separates cropped individual transects by time (interim cell array), 
% Remove noisy ensembles (Plot velocity and depths to determine if 
% individual ensembles need to be removed; lengthy section.)

% Flags internal, surface (if any) and full-NaN ensembles and removes.
% (keeps bottom NaNs)

% Interpolates vel data to grid. (x=25m;z=0.25m) 
% data where 0 is left bank looking up-estuary; create 3d matrices where 
% the third dimension is the transect number.
% **Plots to check bin-avg velocity and see if points are skewed by outliers.

% Deliverables: separate mat file for each line and survey date.
% BI_adcp_L#_mmddyy_cropped_grid_rotate_bininterp.mat: 
% u_bininterp, v_bininterp, along_bininterp, cross_bininterp, time_bininterp,
% h_bininterp(offset added to depth), start (starting index for each transect), 
% last (ending index for each transect), X (bin-averaged projected distance), 
% bins (from concatenated file, cropped to match individual survey dimensions, 
% includes added offset).

%% Run each line / survey date separately.

%% load mat files - Eli's date specific bin file
load BI_bins.mat

%% manually choose line

% line= 1;
% load BI_reference_line_L1.mat

% line= 2;
% load BI_reference_line_L2.mat

line= 3;
load BI_reference_line_L3.mat

% line= 4;
% load BI_reference_line_L4.mat

%% MANUAL choose one survey: EDIT FILE NAMES (#)!!
%% load correct file and CHECK LINE NUMBER!!
% offset corresponds to water line -> ducer

% date= "042023";
% off = 0.4;
% load BI_adcp_L3_042023_cropped_grid_rotate.mat
% load BI_adcp_L3_042023_cropped.mat

% date= "062823";
% off = 0.3;
% load BI_adcp_L3_062823_cropped_grid_rotate.mat
% load BI_adcp_L3_062823_cropped.mat

% date= "071323";
% off = 0.6;
% load BI_adcp_L3_071323_cropped_grid_rotate.mat
% load BI_adcp_L3_071323_cropped.mat

% date= "031524";
% off = 0.6;
% load BI_adcp_L3_031524_cropped_grid_rotate.mat
% load BI_adcp_L3_031524_cropped.mat
% 
date= "031824";
off = 0.6;
load BI_adcp_L3_031824_cropped_grid_rotate.mat
load BI_adcp_L3_031824_cropped.mat

%% redefine 'bins' and 'depth' adding offset

bins = bins + off;
depth = depth + off;

%% use time to find start and end indices

dttime= datetime(time,'ConvertFrom','datenum'); 

%can use diff instead of loop %time(2:end)-time(1:end-1) **later
diffsec= zeros(1,length(time)-1); % one < length
for i = 1:(length(time) - 1)
    diffsec(i) = (time(i+1) - time(i)) * 24 * 3600; %convert to sec, bin 
    % index corresponds to the last data point in line
end
% these are differences in seconds between each bin, now find difference

lastm1= find(diff(diffsec) > 200); % one < last data point for line
lastint= lastm1+1; % last data point in line
last= [lastint, length(time)]; % adding the final transect end
first= lastint+1;
start= [1, first]; % start indices of line (added one)

%% separate; cell array option
% set up cells

tt = 1:length(last); % vector 1xnumber of transects
mt = max(tt); % number of transects

u_t = cell(1, mt);  
v_t = cell(1, mt);
time_t = cell(1, mt);
depth_t = cell(1, mt);
lat_t = cell(1, mt);
lon_t = cell(1, mt);
along_t = cell(1, mt);
cross_t = cell(1, mt);
proj_distance_t = cell(1, mt);

for j= 1:length(start)
    u_t{j} = u(:,start(j):last(j));
    v_t{j} = v(:,start(j):last(j));
    time_t{j} = time(start(j):last(j));
    depth_t{j} = depth(start(j):last(j));
    lat_t{j} = lat(start(j):last(j));
    lon_t{j} = lon(start(j):last(j));
    along_t{j} = along(:,start(j):last(j));
    cross_t{j} = cross(:,start(j):last(j));
    proj_distance_t{j} = proj_distance(start(j):last(j));
end

%% plot cells (all transects from survey date) on map

% ncfilecoast = 'grid_bbay_v7s.nc';
% lon_coast = ncread(ncfilecoast,'lon_coast');
% lat_coast = ncread(ncfilecoast,'lat_coast');
% 
% % plot transects (if they need to be cropped more, return to script
% % No3_data_crop then rerun the rest of the scripts.
% figure color 'w'
% for in=1:max(tt)
% plot(lon_coast, lat_coast,'m')
% ax=gca;
% set(gca,'xlim',[-74.115,-74.085]) %BI
% ax.YLim=[39.745, 39.775]; 
% lon_range= ax.XLim(2)-ax.XLim(1);
% lat_range= ax.YLim(2)-ax.YLim(1);
% lat_range_km = lat_range*111; 
% lon_range_km = cos(lat_range)*111;
% ax.PlotBoxAspectRatio = [1 lat_range/lon_range 1]
% hold on
% plot(lon_t{in},lat_t{in})
% end

%% pcolor for u, v, along, cross velocities

% % adjusting the rows
% pl = ceil(mt/2); %ceil round up
% % cc = ceil(sqrt(mt)); %used to accomodate # of transects, keeping it square.
% % rr = ceil(mt/cc);
% 
% figure('color','w');
% for xx = 1:length(u_t)
%     subplot(pl,2,xx)
%     pcolor(u_t{xx})
%     set(gca,'ydir','reverse','YLim',[0,60])
%     shading flat
%     colormap(redblue)
%     cb=colorbar; %ylabel(cb,'u velocity (m/s)')
%     caxis([-1.5,1.5])
%     title(strcat('transect = ',num2str(tt(xx))))
% end
% sgtitle('u velocity (m/s)');
% 
% figure('color','w');
% for xx = 1:length(u_t)
%     subplot(pl,2,xx)
%     pcolor(v_t{xx})
%     set(gca,'ydir','reverse','YLim',[0,60])
%     shading flat
%     colormap(redblue)
%     cb=colorbar; %ylabel(cb,'v velocity (m/s)')
%     caxis([-1.5,1.5])
%     title(strcat('transect = ',num2str(tt(xx))))
% end
% sgtitle('v velocity (m/s)');
% 
% figure('color','w');
% for xx = 1:length(u_t)
%     subplot(pl,2,xx)
%     pcolor(along_t{xx})
%     set(gca,'ydir','reverse','YLim',[0,60])
%     shading flat
%     colormap(redblue)
%     cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
%     caxis([-1.5,1.5])
%     title(strcat('transect = ',num2str(tt(xx))))
% end
% sgtitle('along channel velocity (m/s)');
% 
% figure('color','w');
% for xx = 1:length(u_t)
%     subplot(pl,2,xx)
%     pcolor(cross_t{xx})
%     set(gca,'ydir','reverse','YLim',[0,60])
%     shading flat
%     colormap(redblue)
%     cb=colorbar; %ylabel(cb,'cross channel velocity (m/s)')
%     caxis([-1.5,1.5])
%     title(strcat('transect = ',num2str(tt(xx))))
% end
% sgtitle('cross channel velocity (m/s)');

%% crop bins
bins = bins(1:height(u));

%% plot problem transects for manual removal of ensembles
%protected in if statements so can run these regardless of selection
%% LINE 2 042023-----------------------------------------------------------
%--------------------------------------------------------------------------
% L2, 4/20/23, Transect 2--------------------------------------------------
if line==2 && date=="042023"
    figure('color','w');
    xx = 2;
    % pcolor(time_t{xx},bins,along_t{xx})
        pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [1,418; 7,420; 2,459; 4,468; 6,419; 7,419];

    mask_a2 = ones(size(along_t{xx}));
    mask_a2(brushedData(1,1):brushedData(2,1),brushedData(1,2)) = NaN;
    mask_a2(brushedData(1,1):brushedData(2,1),brushedData(2,2)) = NaN;
    mask_a2(brushedData(3,1),brushedData(3,2)) = NaN;
    mask_a2(brushedData(4,1),brushedData(4,2)) = NaN;
    mask_a2(brushedData(5,1):brushedData(6,1),brushedData(5,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a2;
    cross_t{xx} = cross_t{xx}.*mask_a2;
    u_t{xx} = u_t{xx}.*mask_a2;
    v_t{xx} = v_t{xx}.*mask_a2;

        %check: scaled
    figure('color','w');
    xx = 2;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 4/20/23, Transect 8--------------------------------------------------
if line==2 && date=="042023"
    figure('color','w');
    xx = 8;
    % pcolor(time_t{xx},bins,along_t{xx})
        pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [76,1.3550; 115,0.9725; 1,147];

    mask_a8 = ones(size(along_t{xx}));
    mask_d8 = ones(size(depth_t{xx}));

    mask_a8(:,brushedData(1,1):brushedData(2,1)) = NaN;
    mask_d8(:,brushedData(1,1):brushedData(2,1)) = NaN;

    mask_a8(brushedData(3,1),brushedData(3,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a8;
    depth_t{xx} = depth_t{xx}.*mask_d8;
    cross_t{xx} = cross_t{xx}.*mask_a8;
    u_t{xx} = u_t{xx}.*mask_a8;
    v_t{xx} = v_t{xx}.*mask_a8;

        %check: scaled
    figure('color','w');
    xx = 8;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 4/20/23, Transect 9--------------------------------------------------
if line==2 && date=="042023"
    figure('color','w');
    xx = 9;
         % pcolor(time_t{xx},bins,along_t{xx})
  pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export
    brushedData = [52,1.0750; 143,4.1600];

    mask_a9 = ones(size(along_t{xx}));
    mask_d9 = ones(size(depth_t{xx}));

    mask_a9(:,brushedData(1):brushedData(2)) = NaN;
    mask_d9(:,brushedData(1):brushedData(2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a9;
    depth_t{xx} = depth_t{xx}.*mask_d9;
    cross_t{xx} = cross_t{xx}.*mask_a9;
    u_t{xx} = u_t{xx}.*mask_a9;
    v_t{xx} = v_t{xx}.*mask_a9;

        %check: scaled
    figure('color','w');
    xx = 9;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------

%% LINE 2 062823-----------------------------------------------------------
%--------------------------------------------------------------------------
% L2, 6/28/23, Transect 1--------------------------------------------------
if line==2 && date=="062823"
    figure('color','w');
    xx = 1;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [1,156;2,156;3,162;4,164];

    mask_a1 = ones(size(along_t{xx}));
    mask_a1(:,brushedData(1,2)) = NaN;
    mask_a1(brushedData(3,1),brushedData(3,2)) = NaN;
    mask_a1(brushedData(4,1),brushedData(4,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a1;
    cross_t{xx} = cross_t{xx}.*mask_a1;
    u_t{xx} = u_t{xx}.*mask_a1;
    v_t{xx} = v_t{xx}.*mask_a1;

    %check: scaled
    figure('color','w');
    xx = 1;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 6/28/23, Transect 2--------------------------------------------------
if line==2 && date=="062823"
    figure('color','w');
    xx = 2;
    % pcolor(along_t{xx})
         pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    % plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [55,2.5183; 102,3.1775];

    mask_a2 = ones(size(along_t{xx}));
    mask_d2 = ones(size(depth_t{xx}));

    mask_a2(:,brushedData(1):brushedData(2)) = NaN;
    mask_d2(:,brushedData(1):brushedData(2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a2;
    depth_t{xx} = depth_t{xx}.*mask_d2;
    cross_t{xx} = cross_t{xx}.*mask_a2;
    u_t{xx} = u_t{xx}.*mask_a2;
    v_t{xx} = v_t{xx}.*mask_a2;

    %check: scaled
    figure('color','w');
    xx = 2;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 6/28/23, Transect 11-------------------------------------------------
if line==2 && date=="062823"
    figure('color','w');
    xx = 11;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [4,85;5,86];

    mask_a11 = ones(size(along_t{xx}));
    mask_a11(brushedData(1,1),brushedData(1,2)) = NaN;
    mask_a11(brushedData(2,1),brushedData(2,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a11;
    cross_t{xx} = cross_t{xx}.*mask_a11;
    u_t{xx} = u_t{xx}.*mask_a11;
    v_t{xx} = v_t{xx}.*mask_a11;

    %check: scaled
    figure('color','w');
    xx = 11;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 6/28/23, Transect 12-------------------------------------------------
if line==2 && date=="062823"
    figure('color','w');
    xx = 12;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [4,55;2,60;4,66];

    mask_a12 = ones(size(along_t{xx}));
    mask_a12(brushedData(1,1),brushedData(1,2)) = NaN;
    mask_a12(brushedData(2,1):brushedData(3,1),brushedData(2,2)) = NaN;
    mask_a12(brushedData(3,1),brushedData(3,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a12;
    cross_t{xx} = cross_t{xx}.*mask_a12;
    u_t{xx} = u_t{xx}.*mask_a12;
    v_t{xx} = v_t{xx}.*mask_a12;

    %check: scaled
    figure('color','w');
    xx = 12;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 6/28/23, Transect 13-------------------------------------------------
if line==2 && date=="062823"
    figure('color','w');
    xx = 13;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [71,1; 71,6];

    mask_a13 = ones(size(along_t{xx}));
    mask_a13(:,brushedData(1):brushedData(2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a13;
    cross_t{xx} = cross_t{xx}.*mask_a13;
    u_t{xx} = u_t{xx}.*mask_a13;
    v_t{xx} = v_t{xx}.*mask_a13;

    %check: scaled
    figure('color','w');
    xx = 13;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end

%% LINE 2 071323-----------------------------------------------------------
%--------------------------------------------------------------------------
% L2, 7/13/23, Transect 5--------------------------------------------------
if line==2 && date=="071323"
    figure('color','w');
    xx = 5;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [29, 30; 4, 210]; 

    mask_a5 = ones(size(along_t{xx}));
    mask_a5(brushedData(1,1), brushedData(1,2)) = NaN;
    mask_a5(brushedData(2,1), brushedData(2,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a5;
    cross_t{xx} = cross_t{xx}.*mask_a5;
    u_t{xx} = u_t{xx}.*mask_a5;
    v_t{xx} = v_t{xx}.*mask_a5;

    %check: scaled
    figure('color','w');
    xx = 5;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 7/13/23, Transect 6--------------------------------------------------
if line==2 && date=="071323"
    figure('color','w');
    xx = 6;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    % plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [2, 172; 3, 172; 2,182]; 

    mask_a6 = ones(size(along_t{xx}));
    mask_a6(brushedData(1,1):brushedData(2,1), brushedData(1,2)) = NaN;
    mask_a6(brushedData(3,1), brushedData(3,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a6;
    cross_t{xx} = cross_t{xx}.*mask_a6;
    u_t{xx} = u_t{xx}.*mask_a6;
    v_t{xx} = v_t{xx}.*mask_a6;

    %check: scaled
    figure('color','w');
    xx = 6;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 7/13/23, Transect 7--------------------------------------------------
if line==2 && date=="071323"
    figure('color','w');
    xx = 7;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    % plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [19,6]; 

    mask_a7 = ones(size(along_t{xx}));
    mask_a7(brushedData(1,1), brushedData(1,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a7;
    cross_t{xx} = cross_t{xx}.*mask_a7;
    u_t{xx} = u_t{xx}.*mask_a7;
    v_t{xx} = v_t{xx}.*mask_a7;

    %check: scaled
    figure('color','w');
    xx = 7;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 7/13/23, Transect 8--------------------------------------------------
if line==2 && date=="071323"
    figure('color','w');
    xx = 8;
    pcolor(along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [2,191; 3,197; 2,216]; 

    mask_a8 = ones(size(along_t{xx}));
    mask_a8(brushedData(1,1), brushedData(1,2)) = NaN;
    mask_a8(brushedData(2,1), brushedData(2,2)) = NaN;
    mask_a8(brushedData(3,1), brushedData(3,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a8;
    cross_t{xx} = cross_t{xx}.*mask_a8;
    u_t{xx} = u_t{xx}.*mask_a8;
    v_t{xx} = v_t{xx}.*mask_a8;

    % check: scaled
    figure('color','w');
    xx = 8;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 7/13/23, Transect 9--------------------------------------------------
if line==2 && date=="071323"
    figure('color','w');
    xx = 9;
    pcolor(along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [26,34;4,193;1,212;3,215;2,218;3,224]; 

    mask_a9 = ones(size(along_t{xx}));
    mask_a9(brushedData(1,1), brushedData(1,2)) = NaN;
    mask_a9(brushedData(2,1), brushedData(2,2)) = NaN;
    mask_a9(:,brushedData(3,2)) = NaN;
    mask_a9(brushedData(4,1), brushedData(4,2)) = NaN;
    mask_a9(brushedData(5,1), brushedData(5,2)) = NaN;
    mask_a9(brushedData(6,1), brushedData(6,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a9;
    cross_t{xx} = cross_t{xx}.*mask_a9;
    u_t{xx} = u_t{xx}.*mask_a9;
    v_t{xx} = v_t{xx}.*mask_a9;

    % check: scaled
    figure('color','w');
    xx = 9;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 7/13/23, Transect 10-------------------------------------------------
if line==2 && date=="071323"
    figure('color','w');
    xx = 10;
    pcolor(along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [5,234]; 

    mask_a10 = ones(size(along_t{xx}));
    mask_a10(brushedData(1,1), brushedData(1,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a10;
    cross_t{xx} = cross_t{xx}.*mask_a10;
    u_t{xx} = u_t{xx}.*mask_a10;
    v_t{xx} = v_t{xx}.*mask_a10;

    % check: scaled
    figure('color','w');
    xx = 10;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 7/13/23, Transect 14-------------------------------------------------
if line==2 && date=="071323"
    figure('color','w');
    xx = 14;
    pcolor(along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [6,190]; 

    mask_a14 = ones(size(along_t{xx}));
    mask_a14(brushedData(1,1), brushedData(1,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a14;
    cross_t{xx} = cross_t{xx}.*mask_a14;
    u_t{xx} = u_t{xx}.*mask_a14;
    v_t{xx} = v_t{xx}.*mask_a14;

    % check: scaled
    figure('color','w');
    xx = 14;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end

%--------------------------------------------------------------------------
% L2, 7/13/23, Transect 15-------------------------------------------------
if line==2 && date=="071323"
    figure('color','w');
    xx = 15;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [17,174;6,337;4,349;5,349;4,375;5,375;4,379]; 

    mask_a15 = ones(size(along_t{xx}));
    mask_a15(brushedData(1,1), brushedData(1,2)) = NaN;
    mask_a15(brushedData(2,1), brushedData(2,2)) = NaN;
    mask_a15(brushedData(3,1):brushedData(4,1), brushedData(3,2)) = NaN;
    mask_a15(brushedData(5,1):brushedData(6,1), brushedData(5,2)) = NaN;
    mask_a15(brushedData(7,1), brushedData(7,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a15;
    cross_t{xx} = cross_t{xx}.*mask_a15;
    u_t{xx} = u_t{xx}.*mask_a15;
    v_t{xx} = v_t{xx}.*mask_a15;

    %check: scaled
    figure('color','w');
    xx = 14;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 7/13/23, Transect 17-------------------------------------------------
if line==2 && date=="071323"
    figure('color','w');
    xx = 17;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [1,193;7,193]; 

    mask_a17 = ones(size(along_t{xx}));
    mask_a17(brushedData(1,1):brushedData(2,1), brushedData(1,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a17;
    cross_t{xx} = cross_t{xx}.*mask_a17;
    u_t{xx} = u_t{xx}.*mask_a17;
    v_t{xx} = v_t{xx}.*mask_a17;

    %check: scaled
    figure('color','w');
    xx = 17;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end

%% LINE 2 031524-----------------------------------------------------------
%--------------------------------------------------------------------------
% L2, 3/15/24, Transect 9--------------------------------------------------
if line==2 && date=="031524"
    figure('color','w');
    xx = 9;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    % clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [1,39;7,42]; 

    mask_a9 = ones(size(along_t{xx}));
    mask_a9(brushedData(1,1):brushedData(2,1), brushedData(1,2):brushedData(2,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a9;
    cross_t{xx} = cross_t{xx}.*mask_a9;
    u_t{xx} = u_t{xx}.*mask_a9;
    v_t{xx} = v_t{xx}.*mask_a9;

    %check: scaled
    figure('color','w');
    xx = 9;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 3/15/24, Transect 10-------------------------------------------------
if line==2 && date=="031524"
    figure('color','w');
    xx = 10;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [1,33;4,33]; 

    mask_a10 = ones(size(along_t{xx}));
    mask_a10(brushedData(1,1):brushedData(2,1), brushedData(1,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a10;
    cross_t{xx} = cross_t{xx}.*mask_a10;
    u_t{xx} = u_t{xx}.*mask_a10;
    v_t{xx} = v_t{xx}.*mask_a10;

    %check: scaled
    figure('color','w');
    xx = 10;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 3/15/24, Transect 11-------------------------------------------------
if line==2 && date=="031524"
    figure('color','w');
    xx = 11;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [1,26;5,26]; 

    mask_a11 = ones(size(along_t{xx}));
    mask_a11(brushedData(1,1):brushedData(2,1), brushedData(1,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a11;
    cross_t{xx} = cross_t{xx}.*mask_a11;
    u_t{xx} = u_t{xx}.*mask_a11;
    v_t{xx} = v_t{xx}.*mask_a11;

    %     %check: scaled
    figure('color','w');
    xx = 11;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end

%% LINE 2 031824-----------------------------------------------------------
%--------------------------------------------------------------------------
% L2, 3/18/24, Transect 1--------------------------------------------------
if line==2 && date=="031824"
    figure('color','w');
    xx = 1;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [4,28;6,30]; 

    mask_a1 = ones(size(along_t{xx}));
    mask_a1(brushedData(1,1), brushedData(1,2)) = NaN;
    mask_a1(brushedData(2,1), brushedData(2,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a1;
    cross_t{xx} = cross_t{xx}.*mask_a1;
    u_t{xx} = u_t{xx}.*mask_a1;
    v_t{xx} = v_t{xx}.*mask_a1;

    % check: scaled
    figure('color','w');
    xx = 1;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L2, 3/18/24, Transect 8--------------------------------------------------
if line==2 && date=="031824"
    figure('color','w');
    xx = 8;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [7,13;8,39]; 

    mask_a8 = ones(size(along_t{xx}));
    mask_a8(brushedData(1,1), brushedData(1,2)) = NaN;
    mask_a8(brushedData(2,1), brushedData(2,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a8;
    cross_t{xx} = cross_t{xx}.*mask_a8;
    u_t{xx} = u_t{xx}.*mask_a8;
    v_t{xx} = v_t{xx}.*mask_a8;

    %check: scaled
    figure('color','w');
    xx = 8;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end

%% LINE 3 042023-----------------------------------------------------------
%--------------------------------------------------------------------------
% L3, 4/20/23, Transect 2--------------------------------------------------
if line==3 && date=="042023"
    figure('color','w');
    xx = 2;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [2:3,194]; %not actually brushing if depth not needed to be filtered, just identify r/c

    mask_a2 = ones(size(along_t{xx}));
    mask_a2(brushedData(1,1):brushedData(1,2), brushedData(1,3)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a2;
    cross_t{xx} = cross_t{xx}.*mask_a2;
    u_t{xx} = u_t{xx}.*mask_a2;
    v_t{xx} = v_t{xx}.*mask_a2;

    %check: scaled
    figure('color','w');
    xx = 2;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L3, 4/20/23, Transect 3--------------------------------------------------
if line==3 && date=="042023"
    figure('color','w');
    xx = 3;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    brushedData = [6,204;6,225]; 

    mask_a3 = ones(size(along_t{xx}));
    mask_a3(brushedData(1,1), brushedData(1,2)) = NaN;
    mask_a3(brushedData(2,1), brushedData(2,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a3;
    cross_t{xx} = cross_t{xx}.*mask_a3;
    u_t{xx} = u_t{xx}.*mask_a3;
    v_t{xx} = v_t{xx}.*mask_a3;

    %check: scaled
    figure('color','w');
    xx = 3;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L3, 4/20/23, Transect 6--------------------------------------------------
if line==3 && date=="042023"
    figure('color','w');
    xx = 6;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    brushedData = [5,91;8,121]; 

    mask_a6 = ones(size(along_t{xx}));
    mask_a6(brushedData(1,1), brushedData(1,2)) = NaN;
    mask_a6(brushedData(2,1), brushedData(2,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a6;
    cross_t{xx} = cross_t{xx}.*mask_a6;
    u_t{xx} = u_t{xx}.*mask_a6;
    v_t{xx} = v_t{xx}.*mask_a6;

    %check: scaled
    figure('color','w');
    xx = 6;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L3, 4/20/23, Transect 7--------------------------------------------------
if line==3 && date=="042023"
    figure('color','w');
    xx = 7;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    brushedData = [5,57;1,70;2,70;1,71;3,93;5,96]; 

    mask_a7 = ones(size(along_t{xx}));
    mask_a7(brushedData(1,1), brushedData(1,2)) = NaN;
    mask_a7(brushedData(2,1), brushedData(2,2)) = NaN;
    mask_a7(brushedData(3,1), brushedData(3,2)) = NaN;
    mask_a7(brushedData(4,1), brushedData(4,2)) = NaN;
    mask_a7(brushedData(5,1), brushedData(5,2)) = NaN;
    mask_a7(brushedData(6,1), brushedData(6,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a7;
    cross_t{xx} = cross_t{xx}.*mask_a7;
    u_t{xx} = u_t{xx}.*mask_a7;
    v_t{xx} = v_t{xx}.*mask_a7;

    %check: scaled
    figure('color','w');
    xx = 7;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L3, 4/20/23, Transect 8--------------------------------------------------
if line==3 && date=="042023"
    figure('color','w');
    xx = 8;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    brushedData = [1,65;7,65;20,65;2,82]; 

    mask_a8 = ones(size(along_t{xx}));
    mask_d8 = ones(size(depth_t{xx}));
    mask_a8(brushedData(2,1):brushedData(3,1),brushedData(2,2)) = NaN;
    mask_a8(brushedData(4,1), brushedData(4,2)) = NaN;
    mask_d8(:,brushedData(1,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a8;
    cross_t{xx} = cross_t{xx}.*mask_a8;
    depth_t{xx} = depth_t{xx}.*mask_d8;
    u_t{xx} = u_t{xx}.*mask_a8;
    v_t{xx} = v_t{xx}.*mask_a8;

    %check: scaled
    figure('color','w');
    xx = 8;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L3, 4/20/23, Transect 9--------------------------------------------------
if line==3 && date=="042023"
    figure('color','w');
    xx = 9;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    brushedData = [3,86; 1,89; 1,93; 1,91; 1,92]; 
    
    mask_a9 = ones(size(along_t{xx}));
    mask_d9 = ones(size(depth_t{xx}));
    mask_a9(brushedData(1,1),brushedData(1,2)) = NaN;
    mask_a9(brushedData(2,1), brushedData(2,2):brushedData(3,2)) = NaN;
    mask_d9(:,brushedData(4,2):brushedData(5,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a9;
    cross_t{xx} = cross_t{xx}.*mask_a9;
    depth_t{xx} = depth_t{xx}.*mask_d9;
    u_t{xx} = u_t{xx}.*mask_a9;
    v_t{xx} = v_t{xx}.*mask_a9;

    %check: scaled
    figure('color','w');
    xx = 9;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L3, 4/20/23, Transect 14-------------------------------------------------
if line==3 && date=="042023"
    figure('color','w');
    xx = 14;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    brushedData = [1,156; 1,200; 4,147; 6,157; 11,159; 7,160; 13,161; 8,162; 
        16,162; 9,163; 60,167; 7,168; 8,169; 11,170; 11,171; 11,173; 11,174;  
     12,175; 13,176; 13,181; 15,182; 15,183; 16,184; 16,189; 17,190; 17,193; 18,194; 18,199; 13,172]; 
    mask_a14 = ones(size(along_t{xx}));
    mask_d14 = ones(size(depth_t{xx}));
    mask_d14(:,brushedData(1,2):brushedData(2,2)) = NaN;
    mask_a14(brushedData(3,1),brushedData(3,2)) = NaN;
    mask_a14(brushedData(4,1):brushedData(5,1), brushedData(4,2):brushedData(5,2)) = NaN;
    mask_a14(brushedData(5,1):brushedData(6,1), brushedData(5,2):brushedData(5,2)) = NaN;
    mask_a14(brushedData(6,1):brushedData(7,1), brushedData(6,2):brushedData(7,2)) = NaN;
    mask_a14(brushedData(8,1):brushedData(9,1), brushedData(8,2):brushedData(9,2)) = NaN;
    mask_a14(brushedData(10,1):brushedData(11,1), brushedData(10,2):brushedData(11,2)) = NaN;
    mask_a14(brushedData(12,1):brushedData(11,1), brushedData(12,2)) = NaN;
    mask_a14(brushedData(13,1):brushedData(11,1), brushedData(13,2)) = NaN;
    mask_a14(brushedData(14,1):brushedData(11,1), brushedData(14,2):brushedData(15,2)) = NaN;
    mask_a14(brushedData(16,1):brushedData(11,1), brushedData(16,2):brushedData(17,2)) = NaN;
    mask_a14(brushedData(18,1):brushedData(11,1), brushedData(18,2)) = NaN;
    mask_a14(brushedData(19,1):brushedData(11,1), brushedData(19,2):brushedData(20,2)) = NaN;
    mask_a14(brushedData(21,1):brushedData(11,1), brushedData(21,2):brushedData(22,2)) = NaN;
    mask_a14(brushedData(23,1):brushedData(11,1), brushedData(23,2):brushedData(24,2)) = NaN;
    mask_a14(brushedData(25,1):brushedData(11,1), brushedData(25,2):brushedData(26,2)) = NaN;
    mask_a14(brushedData(27,1):brushedData(11,1), brushedData(27,2):brushedData(28,2)) = NaN;
    mask_a14(brushedData(29,1):brushedData(11,1), brushedData(29,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a14;
    cross_t{xx} = cross_t{xx}.*mask_a14;
    depth_t{xx} = depth_t{xx}.*mask_d14;
    u_t{xx} = u_t{xx}.*mask_a14;
    v_t{xx} = v_t{xx}.*mask_a14;

%interpolate depth gap (when vel data exists)
depth_int = depth_t{xx};  
x = [155 201];
y = [2.09072550867928 5.51500000000000];
x_interp = 156:200;
y_interp = interp1(x, y, x_interp);
for i = 1:length(x_interp)
    idx = x_interp(i);
    if isnan(depth_int(idx))
        depth_int(idx) = y_interp(i);
    end
end
depth_t{xx} = depth_int;

    %check: scaled
    figure('color','w');
    xx = 14;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end

%% LINE 3 062823-----------------------------------------------------------
%--------------------------------------------------------------------------
% L3, 6/28/23, Transect 1--------------------------------------------------
if line==3 && date=="062823"
    figure('color','w');
    xx = 1;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

   brushedData = [1,147; 5,149; 6,153; 7,153]; 

    mask_a1 = ones(size(along_t{xx}));
    mask_a1(brushedData(1,1):brushedData(2,1), brushedData(1,2):brushedData(2,2)) = NaN;
    mask_a1(brushedData(3,1):brushedData(4,1), brushedData(3,2):brushedData(4,2)) = NaN;
   
    along_t{xx} = along_t{xx}.*mask_a1;
    cross_t{xx} = cross_t{xx}.*mask_a1;
    u_t{xx} = u_t{xx}.*mask_a1;
    v_t{xx} = v_t{xx}.*mask_a1;

    %check: scaled
    figure('color','w');
    xx = 1;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L3, 6/28/23, Transect 2--------------------------------------------------
if line==3 && date=="062823"
    figure('color','w');
    xx = 2;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

   brushedData = [1,88; 5,91; 7,73; 8,73; 6,74; 4,85; 7,97; 9,100; 10,102; 13,118]; 

    mask_a2 = ones(size(along_t{xx}));
    mask_a2(brushedData(1,1):brushedData(2,1), brushedData(1,2):brushedData(2,2)) = NaN;
    mask_a2(brushedData(3,1):brushedData(4,1), brushedData(3,2)) = NaN;
    mask_a2(brushedData(5,1), brushedData(5,2)) = NaN;
    mask_a2(brushedData(6,1), brushedData(6,2)) = NaN;
    mask_a2(brushedData(7,1), brushedData(7,2)) = NaN;
    mask_a2(brushedData(8,1):brushedData(9,1), brushedData(8,2):brushedData(9,2)) = NaN;
    mask_a2(brushedData(10,1), brushedData(10,2)) = NaN;
    
    along_t{xx} = along_t{xx}.*mask_a2;
    cross_t{xx} = cross_t{xx}.*mask_a2;
    u_t{xx} = u_t{xx}.*mask_a2;
    v_t{xx} = v_t{xx}.*mask_a2;

    %check: scaled
    figure('color','w');
    xx = 2;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L3, 6/28/23, Transect 3--------------------------------------------------
if line==3 && date=="062823"
    figure('color','w');
    xx = 3;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

   brushedData = [3,80; 4,80; 19,163]; 

    mask_a3 = ones(size(along_t{xx}));
    mask_a3(brushedData(1,1):brushedData(2,1), brushedData(1,2):brushedData(2,2)) = NaN;
    mask_a3(brushedData(3,1), brushedData(3,2)) = NaN;
   
    along_t{xx} = along_t{xx}.*mask_a3;
    cross_t{xx} = cross_t{xx}.*mask_a3;
    u_t{xx} = u_t{xx}.*mask_a3;
    v_t{xx} = v_t{xx}.*mask_a3;

    %check: scaled
    figure('color','w');
    xx = 3;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L3, 6/28/23, Transect 4--------------------------------------------------
if line==3 && date=="062823"
    figure('color','w');
    xx = 4;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

   brushedData = [4,67; 4,68; 5,69; 8,83; 6,93]; 

    mask_a4 = ones(size(along_t{xx}));
    mask_a4(brushedData(1,1), brushedData(1,2):brushedData(2,2)) = NaN;
    mask_a4(brushedData(3,1), brushedData(3,2)) = NaN;
    mask_a4(brushedData(4,1), brushedData(4,2)) = NaN;
    mask_a4(brushedData(5,1), brushedData(5,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a4;
    cross_t{xx} = cross_t{xx}.*mask_a4;
    u_t{xx} = u_t{xx}.*mask_a4;
    v_t{xx} = v_t{xx}.*mask_a4;

    %check: scaled
    figure('color','w');
    xx = 4;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end

%% LINE 3 031524-----------------------------------------------------------
%--------------------------------------------------------------------------
% L3, 3/15/24, Transect 10-------------------------------------------------
if line==3 && date=="031524"
    figure('color','w');
    xx = 10;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

   brushedData = [17,14; 15,26; 15,27]; 

    mask_a2 = ones(size(along_t{xx}));
    mask_a2(brushedData(1,1), brushedData(1,2)) = NaN;
    mask_a2(brushedData(2,1), brushedData(2,2):brushedData(3,2)) = NaN;
   
    along_t{xx} = along_t{xx}.*mask_a2;
    cross_t{xx} = cross_t{xx}.*mask_a2;
    u_t{xx} = u_t{xx}.*mask_a2;
    v_t{xx} = v_t{xx}.*mask_a2;

    %check: scaled
    figure('color','w');
    xx = 10;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end

%% LINE 4 042023-----------------------------------------------------------
%--------------------------------------------------------------------------
% L4, 4/20/23, Transect 5--------------------------------------------------
if line==4 && date=="042023"
    figure('color','w');
    xx = 5;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [26,107; 15,113]; %not actually brushing if depth not needed to be filtered, just identify r/c
%can manually open matrix to confirm
    mask_a5 = ones(size(along_t{xx}));
    mask_a5(brushedData(1,1) , brushedData(1,2)) = NaN;
    mask_a5(brushedData(2,1) , brushedData(2,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a5;
    cross_t{xx} = cross_t{xx}.*mask_a5;
    u_t{xx} = u_t{xx}.*mask_a5;
    v_t{xx} = v_t{xx}.*mask_a5;

    %check: scaled
    figure('color','w');
    xx = 5;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L4, 4/20/23, Transect 7--------------------------------------------------
if line==4 && date=="042023"
    figure('color','w');
    xx = 7;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [15,108]; 

    mask_a10 = ones(size(along_t{xx}));
    mask_a10(brushedData(1,1), brushedData(1,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a10;
    cross_t{xx} = cross_t{xx}.*mask_a10;
    u_t{xx} = u_t{xx}.*mask_a10;
    v_t{xx} = v_t{xx}.*mask_a10;

    %check: scaled
    figure('color','w');
    xx = 7;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L4, 4/20/23, Transect 10-------------------------------------------------
if line==4 && date=="042023"
    figure('color','w');
    xx = 10;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [21, 73; 10, 60]; 

    mask_a10 = ones(size(along_t{xx}));
    mask_a10(brushedData(1,1), brushedData(1,2)) = NaN;
    mask_a10(brushedData(2,1), brushedData(2,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a10;
    cross_t{xx} = cross_t{xx}.*mask_a10;
    u_t{xx} = u_t{xx}.*mask_a10;
    v_t{xx} = v_t{xx}.*mask_a10;

    %check: scaled
    figure('color','w');
    xx = 10;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end

%% LINE 4 062823-----------------------------------------------------------
%--------------------------------------------------------------------------
% L4, 6/28/23, Transect 7--------------------------------------------------
if line==4 && date=="062823"
    figure('color','w');
    xx = 7;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [10, 14]; 

    mask_a7 = ones(size(along_t{xx}));
    mask_a7(brushedData(1,1), brushedData(1,2)) = NaN;
        along_t{xx} = along_t{xx}.*mask_a7;
    cross_t{xx} = cross_t{xx}.*mask_a7;
    u_t{xx} = u_t{xx}.*mask_a7;
    v_t{xx} = v_t{xx}.*mask_a7;

    %check: scaled
    figure('color','w');
    xx = 7;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L4, 6/28/23, Transect 15-------------------------------------------------
if line==4 && date=="062823"
    figure('color','w');
    xx = 15;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [5, 147; 7, 147]; 

    mask_a15 = ones(size(along_t{xx}));
      mask_a15(brushedData(1,1):brushedData(2,1), brushedData(1,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a15;
    cross_t{xx} = cross_t{xx}.*mask_a15;
    u_t{xx} = u_t{xx}.*mask_a15;
    v_t{xx} = v_t{xx}.*mask_a15;

    %check: scaled
    figure('color','w');
    xx = 15;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end

%% LINE 4 071323-----------------------------------------------------------
%--------------------------------------------------------------------------
% L4, 7/13/23, Transect 4--------------------------------------------------
if line==4 && date=="071323"
    figure('color','w');
    xx = 4;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [1, 144; 22, 166]; 

    mask_a4 = ones(size(along_t{xx}));
    mask_a4(brushedData(1,1):brushedData(2,1), brushedData(1,2):brushedData(2,2)) = NaN;

    along_t{xx} = along_t{xx}.*mask_a4;
    cross_t{xx} = cross_t{xx}.*mask_a4;
    u_t{xx} = u_t{xx}.*mask_a4;
    v_t{xx} = v_t{xx}.*mask_a4;

    %check: scaled
    figure('color','w');
    xx = 4;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L4, 7/13/23, Transect 7--------------------------------------------------
if line==4 && date=="071323"
    figure('color','w');
    xx = 7;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [18,120; 11, 93]; 

    mask_a7 = ones(size(along_t{xx}));
    mask_a7(brushedData(1,1), brushedData(1,2)) = NaN;
    mask_a7(brushedData(2,1), brushedData(2,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a7;
    cross_t{xx} = cross_t{xx}.*mask_a7;
    u_t{xx} = u_t{xx}.*mask_a7;
    v_t{xx} = v_t{xx}.*mask_a7;

    %check: scaled
    figure('color','w');
    xx = 7;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L4, 7/13/23, Transect 9--------------------------------------------------
if line==4 && date=="071323"
    figure('color','w');
    xx = 9;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [13,51]; 

    mask_a9 = ones(size(along_t{xx}));
    mask_a9(brushedData(1,1), brushedData(1,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a9;
    cross_t{xx} = cross_t{xx}.*mask_a9;
    u_t{xx} = u_t{xx}.*mask_a9;
    v_t{xx} = v_t{xx}.*mask_a9;

    %check: scaled
    figure('color','w');
    xx = 9;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end

%% LINE 4 031524-----------------------------------------------------------
%--------------------------------------------------------------------------
% L4, 3/15/24, Transect 1--------------------------------------------------
if line==4 && date=="031524"
    figure('color','w');
    xx = 1;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [1,33;27,33;16,76]; 

    mask_a1 = ones(size(along_t{xx}));
    mask_a1(brushedData(1,1):brushedData(2,1), brushedData(1,2)) = NaN;
    mask_a1(brushedData(1,1):brushedData(3,1), brushedData(3,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a1;
    cross_t{xx} = cross_t{xx}.*mask_a1;
    u_t{xx} = u_t{xx}.*mask_a1;
    v_t{xx} = v_t{xx}.*mask_a1;

    %check: scaled
    figure('color','w');
    xx = 1;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%--------------------------------------------------------------------------
% L4, 7/13/23, Transect 9--------------------------------------------------
if line==4 && date=="031524"
    figure('color','w');
    xx = 9;
    pcolor(along_t{xx})
         % pcolor(time_t{xx},bins,along_t{xx})
      set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
        % plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    plot(depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))

    % data brush: select start point and end point > export (for first run)
    % copied below for subsequent runs
    brushedData = [1,102;17,109]; 

    mask_a9 = ones(size(along_t{xx}));
    mask_a9(brushedData(1,1):brushedData(2,1), brushedData(1,2)) = NaN;
    mask_a9(brushedData(1,1):brushedData(2,1), brushedData(2,2)) = NaN;
    along_t{xx} = along_t{xx}.*mask_a9;
    cross_t{xx} = cross_t{xx}.*mask_a9;
    u_t{xx} = u_t{xx}.*mask_a9;
    v_t{xx} = v_t{xx}.*mask_a9;

    %check: scaled
    figure('color','w');
    xx = 9;
    pcolor(time_t{xx},bins(1:73),along_t{xx})
        % pcolor(along_t{xx})
    set(gca,'ydir','reverse','YLim',[0,40])
    shading flat
    hold on
    plot(time_t{xx},depth_t{xx},'k','LineWidth',2)
    colormap(redblue)
    % cb=colorbar; %ylabel(cb,'along channel velocity (m/s)')
    clim([-1.5,1.5])
    title('along channel velocity (m/s)',strcat('transect = ',num2str(tt(xx))))
end
%% ------------------------------------------------------------------------
%% crop bins (if skipped above)
bins = bins(1:height(u));

%% transect number vector

transect_num = size(time);
for j= 1:length(start) 
    transect_num(start(j):last(j)) = j;
end

num_transects = max(transect_num);  % total number of transects
%just another way to say length(start) 

%% horizontal grid (25m) 
%note: keep dref and X unique to each line.

dx = 0.025; %horizontal bin size (25m)
% x-coordinate for common grid
% 0 should be the left bank looking up-estuary, dref should be the total
% distance to the right bank 
xref=[xo,xw];
yref=[yo,yw];

dref= m_lldist(xref,yref); 
X = 0+dx/2:dx:dref-dx/2; %"proj distance bin avg"
% edges = 0:dx:dref; %if needed...

%% vertical

Z = bins;%1.25:0.25:18;%:1.25:0.25:20;

%% preallocate matrices

along_bininterp = NaN(size(Z,2), length(X), num_transects);  % Z along rows, X along columns
cross_bininterp = NaN(size(Z,2), length(X), num_transects);  
u_bininterp = NaN(size(Z,2), length(X), num_transects); 
v_bininterp = NaN(size(Z,2), length(X), num_transects);  
h_bininterp = NaN(1,length(X), num_transects);
time_bininterp = NaN(1,length(X), num_transects);
Nx_all = zeros(1,length(X), num_transects);
Nx_good = zeros(1,length(X), num_transects); 

%% interpolate internal nans (either this or next section)
% % xx=1;
% % m=340;
% 
% for xx = 1:size(along_t, 2)
%     a_temp = along_t{xx};
%     nbins = size(a_temp, 1);
%     nens = size(a_temp, 2);
% 
%     for m = 1:nens
%         col = a_temp(:, m);
%         idx_valid = find(~isnan(col));
% 
%         if numel(idx_valid) < 2
%             continue  % not enough data to interpolate
%         end
%         idx_start = idx_valid(1);
%         idx_end = idx_valid(end);
%         %or:
%         % idx_start_interp = find(~isnan(col),1,'first');
%         % idx_end_interp = find(~isnan(col),1,'last');
% 
%         interp_range = idx_start:idx_end;
%         interp_vals = interp1(idx_valid, col(idx_valid), interp_range);
% 
%         col(idx_start:idx_end) = interp_vals;
% 
%         a_temp(:, m) = col;
%     end
%     along_t{xx} = a_temp;
% end

%% flag nans (internal, surface, optional full-nan ensembles)
%bottom nans are not touched
flag = -999;

for xx = 1:size(along_t, 2)
    a_temp = along_t{xx};
    nbins = size(a_temp, 1);
    nens = size(a_temp, 2);

    for m = 1:nens
        col = a_temp(:, m);
        idx_valid = find(~isnan(col));

        % flag full nan ensembles
        if isempty(idx_valid)
            col(:) = flag;
            a_temp(:, m) = col;
            continue
        end

        idx_start = idx_valid(1);
        idx_end = idx_valid(end);

        % flag internal nans
        internal_range = idx_start:idx_end;
        internal_nan_idx = internal_range(isnan(col(internal_range)));
        col(internal_nan_idx) = flag;

        % flag surface nans
        surface_range = 1:(idx_start-1);
        surface_nan_idx = surface_range(isnan(col(surface_range)));
        col(surface_nan_idx) = flag; 

        a_temp(:, m) = col;
    end

    along_t{xx} = a_temp;
end

%% INTERPOLANT
%add cross, u, v

for m= 1:num_transects
    along_1=along_t{m};
    cross_1=cross_t{m};
    u_1=u_t{m};
    v_1=v_t{m};
    proj_1 = proj_distance_t{m};
    h_1 = depth_t{m};
    t_1 = time_t{m};

        [x,y] = meshgrid(proj_1,bins);
        x2=x(:);
        y2=y(:);
        along_12=along_1(:);
        cross_12=cross_1(:);
        u_12=u_1(:);
        v_12=v_1(:);
        
        % %remove flagged data (surface/internal/edge nans; just not removing 
        % %nans below last good data)
        x2(along_1==-999)=[];
        y2(along_1==-999)=[];
        along_12(along_1==-999)=[];
        cross_12(along_1==-999)=[];
        u_12(along_1==-999)=[];
        v_12(along_1==-999)=[];

        % %remove all nans
        % b=find(isnan(along_12));
        % along_12(b) = [];
        % x2(b) = [];
        % y2(b) = [];

        F = scatteredInterpolant(x2, y2, along_12);
        F.Method = 'natural';
        F.ExtrapolationMethod = 'none';
        Fc = scatteredInterpolant(x2, y2, cross_12);
        Fc.Method = 'natural';
        Fc.ExtrapolationMethod = 'none';
        Fu = scatteredInterpolant(x2, y2, u_12);
        Fu.Method = 'natural';
        Fu.ExtrapolationMethod = 'none';
        Fv = scatteredInterpolant(x2, y2, v_12);
        Fv.Method = 'natural';
        Fv.ExtrapolationMethod = 'none';

        % crop horiz range X by index...
        xmin = min(proj_1);
        xmax = max(proj_1);
        X_idx = find(X >= xmin & X <= xmax);
        X_crop = X(X_idx); %crop x by index of range

        [X_grd, Z_grd] = meshgrid(X_crop, Z);
 
       h_bininterp_crop = interp1(proj_1, h_1, X_crop);%, 'linear', NaN);
       t_bininterp_crop = interp1(proj_1, t_1, X_crop);
       along_bininterp_crop = F(X_grd, Z_grd);
       cross_bininterp_crop = Fc(X_grd, Z_grd);
       u_bininterp_crop = Fu(X_grd, Z_grd);
       v_bininterp_crop = Fv(X_grd, Z_grd);

along_bininterp_temp = nan(size(Z,2), length(X));%;, num_transects);  % Z along rows, X along columns
cross_bininterp_temp = nan(size(Z,2), length(X));%;, num_transects);  % Z along rows, X along columns
u_bininterp_temp = nan(size(Z,2), length(X));%;, num_transects);  % Z along rows, X along columns
v_bininterp_temp = nan(size(Z,2), length(X));%;, num_transects);  % Z along rows, X along columns
h_bininterp_temp = nan(1, length(X));
t_bininterp_temp = nan(1, length(X));

along_bininterp_temp(:, X_idx) = along_bininterp_crop; %insert interpolated values in the correct X columns
cross_bininterp_temp(:, X_idx) = cross_bininterp_crop; %insert interpolated values in the correct X columns
u_bininterp_temp(:, X_idx) = u_bininterp_crop; %insert interpolated values in the correct X columns
v_bininterp_temp(:, X_idx) = v_bininterp_crop; %insert interpolated values in the correct X columns
h_bininterp_temp(1,X_idx) = h_bininterp_crop;
t_bininterp_temp(1,X_idx) = t_bininterp_crop;

%new; could remove
% shp = alphaShape(x2,y2); %
% ca = criticalAlpha(shp,'one-region');
% shp.Alpha = ca;
% along_grd(~inShape(shp,X_grd(:),Z_grd(:)))=NaN;

along_bininterp(:,:,m) = along_bininterp_temp;
cross_bininterp(:,:,m) = cross_bininterp_temp;
u_bininterp(:,:,m) = u_bininterp_temp;
v_bininterp(:,:,m) = v_bininterp_temp;
h_bininterp(1,:,m) = h_bininterp_temp;
time_bininterp(1,:,m) = t_bininterp_temp;

% Nx_all(:,n,m) = length(a); %use to find number of ensembles (including bad ens.)
% Nx_good(:,n,m) = length(find(~isnan(mean(u_t{m}(:,a),1,'omitnan')))); %check # of good velocity points

clear x y %along_1 proj_1 x y F xmin xmax X_idx X_crop X_grd Z_grd along_grd_crop
end

% figure;
% imagesc(X, Z, a1a);
% set(gca, 'YDir', 'reverse');
% xlabel('Distance (m)');
% ylabel('Depth (m)');
% title('Interpolated Along-channel Velocity');
% colorbar;

%% Plot depth-averaged, bin-interped velocity to check for points skewed by
% outliers. Use bar graph to check data coverage.
% need to adjust for interpolated nx_good/nx_all

% for xx = 1:size(along_bininterp,3)
% figure('color','w')
% subplot(4,1,1)
% plot(proj_distance_t{xx},mean(u_t{xx},1,'omitnan'))
% hold on
% plot(X,mean(u_bininterp(:,:,xx),1,'omitnan'),'ro--')
% ylabel('Depth-averaged u')
% title(sprintf('dx = %.1f m\n',dx*1e3))
% sgtitle(strcat('transect = ',num2str(tt(xx))))
% legend('original','bin-interpolated')
% 
% subplot(4,1,2)
% plot(proj_distance_t{xx},mean(v_t{xx},1,'omitnan'))
% hold on
% plot(X,mean(v_bininterp(:,:,xx),1,'omitnan'),'ro--')
% ylabel('Depth-averaged v')
% % title(sprintf('dx = %.1f m\n',dx*1e3))
% 
% subplot(4,1,3)
% plot(proj_distance_t{xx},depth_t{xx},'k'); hold on
% plot(X,h_bininterp(:,:,xx),'m*:')
% ylabel('bin-avg depth')
% set(gca,'ydir','reverse','YLim',[0,max(depth)])
% 
% subplot(4,1,4)
% bar(X,Nx_good(:,:,1))
% ylabel('# good ensembles per horiz bin')
% xlabel('Projected distance')
% 
%     % export_fig(fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\fig\Bin-Avg\L2\031524\skewed', ...
%     % sprintf('transect%d.png', xx)), '-png','-m2');
% end

%% load bin-avg data for comparison

% avg = load('BI_adcp_L3_062823_cropped_grid_rotate_binavg.mat');
% avg2 = load('BI_adcp_L3_062823_cropped_grid_rotate_binavg_extrap.mat');

% avg = load('BI_adcp_L3_042023_cropped_grid_rotate_binavg.mat');
% avg2 = load('BI_adcp_L3_042023_cropped_grid_rotate_binavg_extrap.mat');

%% figures to compare: clean data, bin-avg data, bin-interp data
% % mplot=1;
% for mplot = 1:num_transects  
% figure('color','w');
% subplot(1,3,1);
% pcolorjw(proj_distance_t{mplot},bins,along_t{mplot})%73x359
% % shading flat
% hold on
% plot(proj_distance_t{mplot},depth_t{mplot},'k','LineWidth',2)
% set(gca, 'YDir', 'reverse');
% xlabel('ens');
% ylabel('bin');
% title('Along Velocity');
% colormap('redblue')
% colorbar;
% clim([-1.5,1.5])
% xlim([0,0.4])
% 
% subplot(1,3,2);
% pcolorjw(X, bins, avg.along_binavg(:,:,mplot));  % from original averaging
% % shading flat
% hold on 
% plot(X,avg.h_binavg(:,:,mplot),'k','LineWidth',2)
% set(gca, 'YDir', 'reverse');
% xlabel('Distance (m)');
% ylabel('Depth (m)');
% title('Bin-Averaged Along Velocity');
% colormap('redblue')
% colorbar;
% clim([-1.5,1.5])
% xlim([0,0.4])
% 
% subplot(1,3,3);
% pcolorjw(X, Z, along_bininterp(:,:,mplot))  % from interpolated
% % shading flat
% hold on
% plot(X,h_bininterp(1,:,mplot),'k','LineWidth',2)
% set(gca, 'YDir', 'reverse');
% xlabel('Distance (m)');
% ylabel('Depth (m)');
% title('Interpolated Along Velocity');
% colormap('redblue')
% colorbar;
% clim([-1.5,1.5])
% xlim([0,0.4])
% sgtitle(sprintf('Transect %d', mplot)); 
% 
% % filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No5_bin_interp\QA', ...
% %     sprintf('Bin-Interp_QA_%s_line%d_transect%d', date, line, mplot));
% % export_fig([filename, '.png'], '-m2');
% % savefig([filename, '.fig']);
% end


%% figures to compare: clean data and bin-interp data
% % mplot=1;
% for mplot = 1:num_transects  
% figure('color','w');
% subplot(1,2,1);
% pcolorjw(proj_distance_t{mplot},bins,along_t{mplot})%73x359
% % shading flat
% hold on
% plot(proj_distance_t{mplot},depth_t{mplot},'k','LineWidth',2)
% set(gca, 'YDir', 'reverse');
% xlabel('ens');
% ylabel('bin');
% title('Along Velocity');
% colormap('redblue')
% colorbar;
% clim([-1.5,1.5])
% xlim([0,0.4])
% 
% subplot(1,2,2);
% pcolorjw(X, Z, along_bininterp(:,:,mplot))  % from interpolated
% % shading flat
% hold on
% plot(X,h_bininterp(1,:,mplot),'k','LineWidth',2)
% set(gca, 'YDir', 'reverse');
% xlabel('Distance (m)');
% ylabel('Depth (m)');
% title('Interpolated Along Velocity');
% colormap('redblue')
% colorbar;
% clim([-1.5,1.5])
% xlim([0,0.4])
% sgtitle(sprintf('Transect %d', mplot)); 
% 
% filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No5_bin_interp', ...
%     sprintf('RawvsInterp_%s_line%d_transect%d', date, line, mplot));
% export_fig([filename, '.png'], '-m2');
% savefig([filename, '.fig']);
% 
% end

%% figures to compare CROSS
% % mplot=1;

% for mplot = 1:num_transects  
% figure('color','w');
% subplot(1,3,1);
% pcolorjw(proj_distance_t{mplot},bins,cross_t{mplot})%73x359
% % shading flat
% hold on
% plot(proj_distance_t{mplot},depth_t{mplot},'k','LineWidth',2)
% set(gca, 'YDir', 'reverse');
% xlabel('ens');
% ylabel('bin');
% title('Cross Velocity');
% colormap('redblue')
% colorbar;
% clim([-1.5,1.5])
% xlim([0,0.4])
% 
% subplot(1,3,2);
% pcolorjw(X, bins, og.cross_binavg(:,:,mplot));  % from original averaging
% % shading flat
% hold on 
% plot(X,og.h_binavg(:,:,mplot),'k','LineWidth',2)
% set(gca, 'YDir', 'reverse');
% xlabel('Distance (m)');
% ylabel('Depth (m)');
% title('Bin-Averaged Cross Velocity');
% colormap('redblue')
% colorbar;
% clim([-1.5,1.5])
% xlim([0,0.4])
% 
% subplot(1,3,3);
% pcolorjw(X, Z, cross_bininterp(:,:,mplot))  % from interpolated
% % shading flat
% hold on
% plot(X,h_bininterp(1,:,mplot),'k','LineWidth',2)
% set(gca, 'YDir', 'reverse');
% xlabel('Distance (m)');
% ylabel('Depth (m)');
% title('Interpolated Cross Velocity');
% colormap('redblue')
% colorbar;
% clim([-1.5,1.5])
% xlim([0,0.4])
% sgtitle(sprintf('Transect %d', mplot)); 
% end


%% save 

% save(strcat('BI_','adcp_','L',num2str(line),'_',...
% date,'_cropped_','grid_','rotate_','bininterp','.mat'),'along_bininterp','cross_bininterp',...
% 'u_bininterp','v_bininterp','h_bininterp','time_bininterp','start','last','X','bins')
