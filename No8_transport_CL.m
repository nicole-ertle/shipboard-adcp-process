%% No8_transport_CL.m
% NSE
% calculates transport from residual velocities obtained from tide fitting


%% SELECT LINE MANUALLY
line=2;
%% ------------------------------------------------------------------------
%% CHANGE THE L# in MAT FILES BELOW:
%% 04-20-23
date = "042023";
load('BI_adcp_L2_042023_tide_fit_bin-avg_sigma_CL.mat')
%% ------------------------------------------------------------------------
%% 06-28-23
% date = "062823";
% load('BI_adcp_L4_062823_tide_fit_bin-avg_sigma_CL.mat')
%% ------------------------------------------------------------------------
%% 07-13-23
% date = "071323";
% load('BI_adcp_L4_071323_tide_fit_bin-avg_sigma_CL.mat')
%% ------------------------------------------------------------------------
%% 03-15-24
% date = "031524";
% load('BI_adcp_L4_031524_tide_fit_bin-avg_sigma_CL.mat')
%% ------------------------------------------------------------------------
%% 03-18-24
% date = "031824";
% load('BI_adcp_L4_031824_tide_fit_bin-avg_sigma_CL.mat')
%% ------------------------------------------------------------------------
%% extract residual velocity
vel_resid = squeeze(T_fit_result_all3(1, :, :));

%% transport = velocity * area of bin

l = (Xgrid(1,2) - Xgrid(1,1)) * 1000; %km to m
w_real = h_mean/9; %avg depth by #sig levels -1 
a_real = repmat(l * w_real, height(Xgrid), 1);
transport = vel_resid .* a_real;
transport_sum = nansum(nansum(transport));

%%
%excel each line, each date - transport through
%try pcolor
figure('color','w')
pcolor(Xgrid, Z_real, transport); shading flat 
hold on
plot(Xgrid(1,:),h_mean,'k','LineWidth',2)
colormap(redblue)
    cb=colorbar; caxis([-5,5])
xlabel('Distance (km)')
ylabel(cb,'Transport (m^3/s)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
xlim([0,Xgrid(end)+0.025]) 
% title('Transport')
title(strcat('Line=', num2str(line),'  Date=',(num2str(date))))

 filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No8_transport\constant-linear', ...
    sprintf('Transport_%s_line%d', date, line));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);


%% SAVE variables
save(strcat('BI_','adcp_','L',num2str(line),'_',date,...
    '_transport_bin-avg_CL-extrap_sigma-tide-fit','.mat'),'vel_resid','transport','transport_sum')





