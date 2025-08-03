% No7_bin_interp_sigma_tide_fit.m
% NSE
% performs tide fit the the x-bin interpolated, z-extrapolated velocity data
% (all sigma levels)

%Deliverables: BI_adcp_L#_DATE_tide_fit_bin-interp_sigma.mat
%T_fit_result_all=tide fit from M2
%T_fit_result_all2=tide fit from M2+M4
%T_fit_result_all3=tide fit from M2+M4+M6
%other related variables; see save in last section.
%% ------------------------------------------------------------------------
%% SELECT LINE MANUALLY
line=3;
%% ------------------------------------------------------------------------
%% CHANGE THE L# in MAT FILES BELOW:
%% 04-20-23
date = "042023";
load('BI_adcp_L3_042023_cropped_grid_rotate_bininterp_extrap.mat')
%% ------------------------------------------------------------------------
%% 06-28-23
% date = "062823";
% load('BI_adcp_L3_062823_cropped_grid_rotate_bininterp_extrap.mat')
%% ------------------------------------------------------------------------
%% 07-13-23
% date = "071323";
% load('BI_adcp_L3_071323_cropped_grid_rotate_bininterp_extrap.mat')
%% ------------------------------------------------------------------------
%% 03-15-24
% date = "031524";
% load('BI_adcp_L3_031524_cropped_grid_rotate_bininterp_extrap.mat')
%% ------------------------------------------------------------------------
%% 03-18-24
% date = "031824";
% load('BI_adcp_L3_031824_cropped_grid_rotate_bininterp_extrap.mat')
%% ------------------------------------------------------------------------
%% plots--need to adapt this,,,
%test: and prep: for loop graph time and horiz vel transect number assign to time
% tstar and depth interp vel and see if fit ok
% 
% for xx = 1:size(along_extrap,3)
% figure('color', 'white')
%     subplot(2,1,1)
% plot(time_bininterp(1,:,xx),along_avg(1,:,xx),'k*-')
% xlabel('time')
% ylabel('velocity')
% title(strcat('Transect =', num2str(xx)));
%     subplot(2,1,2)
% plot(X,along_avg(1,:,xx),'k*-')
% xlabel('distance')
% ylabel('velocity')
% end
%% ------------------------------------------------------------------------
%% plots all transects (see how many pts per section)
% % tt=1:length(X);
% figure('color', 'white')
% for xx = 1:size(along_sigma,3)
% plot(X,along_sigma(1,:,xx))
% % plot(time_bininterp(1,:,xx),along_sigma(:,:,xx))
% hold on
% xlabel('distance')
% ylabel('velocity')
% title('All Transects')
% legend
% % legend(strcat('transect = ',num2str(tt(xx))))
% end
%% ------------------------------------------------------------------------
%% low tide times using VERIFIED water level data from station Atlantic City (NOAA)

if date== "042023"
LT_time_datetime = datetime([2023,4,20,6,24,00; 2023,4,20,18,36,00; 2023,4,21,7,12,00])';
LT_time = datenum(LT_time_datetime);
end

if date== "062823"
LT_time_datetime = datetime([2023,6,28,1,42,00; 2023,6,28,13,24,00; 2023,6,29,2,30,00])';
LT_time = datenum(LT_time_datetime);
end

if date== "071323"
LT_time_datetime = datetime([2023,7,13,2,54,00; 2023,7,13,14,48,00; 2023,7,14,3,54,00])';
LT_time = datenum(LT_time_datetime);
end

if date== "031524"
LT_time_datetime = datetime([2024,3,15,10,18,00; 2024,3,15,22,12,00; 2024,3,16,11,30,00])';
LT_time = datenum(LT_time_datetime);
end

if date== "031824"
LT_time_datetime = datetime([2024,3,18,00,24,00; 2024,3,18,13,42,00; 2024,3,19,1,36,00])';
LT_time = datenum(LT_time_datetime);
end
%% ------------------------------------------------------------------------
% pcolor n=7
% figure('color','w')
% pcolorjw(time_bininterp(1,7,:),-z,squeeze(along_extrap(:,7,:)))
%% ------------------------------------------------------------------------
%% Convert time to fraction of tide cycle (Tstar)

Tstar=zeros(size(time_bininterp));

for xx = 1:size(along_extrap, 3)
    for n=1:length(X)
        prior= find(LT_time<=time_bininterp(1,n,xx),1,'last');
        Tstar(1,n,xx) = (time_bininterp(1,n,xx)-LT_time(prior))/(LT_time(prior+1)-LT_time(prior));
    end
end

%Turn 0s back into NaNs
Tstar(Tstar == 0) = NaN;

Tstar1= Tstar(:,:,1);

% reshape
Tstar_re = squeeze(Tstar)';
%% ------------------------------------------------------------------------
%% create min criteria for transects
n_bins = size(along_sigma, 1); % vertical bins
n_ens = size(along_sigma, 2);  % ensembles
n_tran = size(along_sigma, 3); % transects

along_sigma_clean = along_sigma; % preallocate "cleaned" arrays; 
% these are with the data removed for transects below min criteria.
Tstar_clean = Tstar;
h_bininterp_clean = h_bininterp;
time_bininterp_clean = time_bininterp;

for n = 1:n_ens  % loop over ensembles (X locations)
    valid_transects = false(n_tran, 1);  % track which transects have any valid data for this ensemble
                                           % logical 0s
    for t = 1:n_tran
        velocity_column = along_sigma(:, n, t);  % all bins at this ensemble + transect

        if any(~isnan(velocity_column))
            valid_transects(t) = true;  % this transect has data for this ensemble
        end % transect # for this ensemble
    end

    num_valid_transects = sum(valid_transects);

    if num_valid_transects < 7
        % remove all data for this ensemble (across all bins and transects)
        along_sigma_clean(:,n,:) = NaN;
        Tstar_clean(1,n,:) = NaN;
        h_bininterp_clean(1,n,:) = NaN;
        time_bininterp_clean(1,n,:) = NaN;

        fprintf('Removed ensemble %d (only %d valid transects: [%s])\n', ...
            n, num_valid_transects, num2str(find(valid_transects)'));
    end
end
%% ------------------------------------------------------------------------
%% 1st run: M2
T_fit_result_all = NaN(3, n_bins, n_ens);  
T_fit_rms_all = NaN(n_bins, n_ens);

Tc_names = [{'M2'}];
Tc_periods = [12.4208]; 
Mx_freq = 24*(LT_time(2)-LT_time(1))./Tc_periods; %per tidal cycle (calculated by elapsed time between low tides divided by the tidal period of each harmonic (M2=12.4208 hours, M4=6.2103 hours, M6=4.1402 hours)

for z = 1:n_bins
    for n = 1:n_ens
        % extract data across transects for this bin/ensemble
        vel = squeeze(along_sigma_clean(z,n,:));     % transects x 1
        tstar = squeeze(Tstar_clean(1,n,:));    % transects x 1
    
        Pars = [nanmean(vel); 0; (max(vel) - min(vel)) / 2]; % initial guess: [residual; phase; amplitude]
        if isnan(Pars(1))
            Pars(2)=NaN; % turn 0 to NaNs in phasing 
        end
        if Pars(1:3)==0
            Pars(1:3)=NaN; %change the bottom from 0 to NaNs; comment out if you want to show bottom 0s, then adjust next Tc fits
        end
        T_guess_all(:, z, n) = Pars;

        r= find(isnan(vel));
            vel(r) = [];
            tstar(r) = [];
        rr= find(isnan(Pars));
            Pars(rr) = [];
            if isempty(Pars)
                continue
            else

        fit = fminsearch(@adcp_tide_fit, Pars, [], Mx_freq, tstar, vel); % run tide fit

        % adjust for negative amplitude
        if fit(3) < 0
            fit(3) = -fit(3);
            fit(2) = fit(2) + pi;
        end

        T_fit_result_all(:,z,n) = fit; % store results

        % compute RMS error and sum of squares
        pred = fit(1) + fit(3)*sin(2*pi*Mx_freq*tstar - fit(2));
        T_fit_rms_all(z,n) = sqrt(nanmean((vel - pred).^2));
        ss_all(z,n) = nansum((vel - pred).^2);
    end
    end
end
%% ------------------------------------------------------------------------
%% export m2
[num_params, n_bins, n_ens] = size(T_fit_result_all);
results_flat = reshape(permute(T_fit_result_all, [3, 2, 1]), [], num_params);
ens_list = repmat((1:n_ens)', n_bins, 1);%1–16, repeated 10x
bin_list = repelem((1:n_bins)', n_ens); %1repeated16x,2repeated16x,etc
fit_table_M2 = array2table([ens_list, bin_list, results_flat], ...
    'VariableNames', {'Ensemble', 'Bin', 'Residual', 'Phase_M2', 'Amp_M2'});
sheetName = date + "_M2";
base = 'C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\excel';
filename = 'tidefitresult.xlsx';
output_path = fullfile(base, filename);
writetable(fit_table_M2, output_path, 'Sheet', sheetName);
%% ------------------------------------------------------------------------
%% export m2 rms
rms_flat_M2 = T_fit_rms_all(:);  % [n_bins*n_ens x 1]
rms_table_M2 = table(ens_list, bin_list, rms_flat_M2, ...
    'VariableNames', {'Ensemble', 'Bin', 'RMS_M2'});
sheetName = date + "_RMS_M2";
base = 'C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\excel';
filename = 'rms.xlsx';
output_path = fullfile(base, filename);
writetable(rms_table_M2, output_path, 'Sheet', sheetName);
%% ------------------------------------------------------------------------
%% Tc plot single ens/bin M2
ee = 8;
zz = 7;

Pars = T_guess_all(:, zz, ee);
T_fit_result = T_fit_result_all(:, zz, ee);  
Tstar_re = squeeze(Tstar_clean(1, ee, :));
along_re = squeeze(along_sigma_clean(zz, ee, :));

% think about removing this:
% valid = ~isnan(Tstar_re) & ~isnan(along_re);
% Tstar_re = Tstar_re(valid);
% along_re = along_re(valid);
% if isempty(Tstar_re) || any(isnan(T_fit_result))
%     disp('No valid data or complete fit result at this bin/ensemble.')
%     return
% end
%
N_tc = length(Mx_freq);
tide_time = 0:0.01:1;

figure('color','w')
qp = NaN(N_tc, length(tide_time));
qp_in = NaN(N_tc, length(tide_time));
for m = 1:N_tc
    subplot(N_tc+1, 1, m)
    % fit result
    amp_fit = T_fit_result(2*m+1);
    phase_fit = T_fit_result(2*m);
    % input parameters
    amp_in = Pars(2*m+1);
    phase_in = Pars(2*m);
qp_in(m,:) = amp_in * sin(2 * pi * Mx_freq(m) * tide_time - phase_in);
qp(m,:) = amp_fit * sin(2 * pi * Mx_freq(m) * tide_time - phase_fit);
plot(tide_time, qp_in(m,:), 'Color', [0.5, 0.5, 0.5]); hold on
plot(tide_time, qp(m,:), 'b')
ylabel('Amplitude (m)'); xlabel('T*'); title(Tc_names{m})
    %may not need:
    % if isnan(amp_fit) || isnan(phase_fit)
    %     title(sprintf('%s (NaN)', Tc_names{m}))
    %     continue
    % end
    %
end

% composite signal
subplot(N_tc+1, 1, m+1)
model_sum = sum(qp, 1) + T_fit_result(1);
input_sum = sum(qp_in, 1) + Pars(1);

plot(tide_time, input_sum, 'Color', [0.5, 0.5, 0.5]); hold on
plot(tide_time, model_sum, 'b')
plot(Tstar_re, along_re, 'r*')
ylabel('Amplitude (m)'); xlabel('T*'); title('Sum of Tidal Constituents')
legend('Input parameters','Model Fit (fminsearch)','Observed Data')
sgtitle(sprintf('Ensemble = %d, Bin = %d', ee, zz))

filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\constituents', ...
    sprintf('M2_%s_line%d_ens%d_bin%d', date, line,ee,zz));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% pcolor sigma resid vel M2

residuals = squeeze(T_fit_result_all(1, :, :));  %10x16
[Xgrid, Zgrid] = meshgrid(X, sigma_levels); 

figure('color','w')
pcolorjw(Xgrid, Zgrid, residuals);  
% shading interp
    colormap(redblue)
    cb=colorbar; 
    ylabel(cb,'Residual Velocity (m/s)')
    caxis([-0.5,0.5])
xlabel('Distance (km)')
ylabel('Sigma Levels')
xlim([0,X(end)+0.025])
% xlim([0,max(X)])
title(strcat('Line=', num2str(line),'  Date=',(num2str(date))))
% 
filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit', ...
    sprintf('residual_velocity_M2_%s_line%d_sigma', date, line));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% pcolor resid vel vs average bathy M2

h_mean = squeeze(nanmean(h_bininterp_clean,1)); %[sections, #transects]
h_mean = nanmean(h_mean, 2); %added 2, since want mean of all transects at this ensemble (not mean of single transect at all ens).
h_mean=h_mean';
N_levels = length(sigma_levels);
n_ens = length(h_mean);

Z_real = NaN(N_levels, n_ens);  
for n = 1:n_ens
    h = h_mean(n);
    if isnan(h), continue, end
    Z_real(:, n) = sigma_levels * h; 
end
Z_real = flipud(Z_real); %orient from surface to bottom to match other matrices

figure('color','w')
pcolor(X, Z_real, residuals)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
    colormap(redblue)
    cb=colorbar
    ylabel(cb,'Residual Velocity (m/s)')
    caxis([-0.5,0.5])
xlim([0,X(end)+0.025])
xlabel('Distance (km)')
ylabel('Depth (m)')
title(strcat('Line=', num2str(line),'  Date=',(num2str(date))))
set(gca, 'YDir', 'reverse')  
% 
filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\depth contour', ...
    sprintf('residual_velocity_M2_%s_line%d_depth', date, line));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% pcolor RMS M2

figure('color','w')
pcolorjw(Xgrid, Zgrid, T_fit_rms_all)
% shading flat
    cb=colorbar
    caxis([0,1])
xlim([0,X(end)+0.025])
xlabel('Distance (km)')
ylabel('Sigma Levels')
title(strcat('Line=', num2str(line),'  Date=',(num2str(date))))

filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\rmse', ...
    sprintf('RMSE_M2_%s_line%d_depth', date, line));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
%% 2nd run: M2+M4
Tc_names = {'M2', 'M4'};
Tc_periods = [12.4208, 6.2103];
Mx_freq = 24*(LT_time(2)-LT_time(1)) ./ Tc_periods;
N_tc = length(Mx_freq);

% preallocate storage for second fit
T_fit_result_all2 = NaN(2*N_tc+1, n_bins, n_ens);  
T_fit_rms_all2 = NaN(n_bins, n_ens);
ss_all2 = NaN(n_bins, n_ens);
T_guess_all2 = NaN(2*N_tc + 1, n_bins, n_ens);

for z = 1:n_bins
    for n = 1:n_ens
        fit_M2 = T_fit_result_all(:,z,n);
        % if any(isnan(fit_M2)), continue, end

        % extract data across transects
        vel = squeeze(along_sigma_clean(z,n,:));
        tstar = squeeze(Tstar_clean(1,n,:));

        % construct initial guess for M2 + M4
        Pars2 = NaN(2*N_tc+1, 1);
        Pars2(1) = fit_M2(1);  % residual
        Pars2(2) = fit_M2(2);  % M2 phase
        Pars2(3) = fit_M2(3);  % M2 amp
        Pars2(4) = fit_M2(2) - (-0.16888); % M4 phase (offset based on?)
        Pars2(5) = fit_M2(3) * 0.0125; % M4 amp

        if all(Pars2(1:3) == 0) 
           Pars2(1:5) = NaN; %not sure if this is needed anymore
        end

        T_guess_all2(:, z, n) = Pars2;

        r= find(isnan(vel));
            vel(r) = [];
            tstar(r) = [];
        rr= find(isnan(Pars2));
            Pars2(rr) = [];
            if isempty(Pars2)
                continue
            else

        fit2 = fminsearch(@adcp_tide_fit, Pars2, [], Mx_freq, tstar, vel);  % second fit: M2 + M4 

            if fit2(3) < 0 % correct negative amplitudes
                fit2(3) = -fit2(3);
                fit2(2) = fit2(2) + pi;
            end
            if fit2(5) < 0
                fit2(5) = -fit2(5);
                fit2(4) = fit2(4) + pi;
            end

        T_fit_result_all2(:,z,n) = fit2;

        total = NaN(N_tc, length(tstar)); % reconstruct model and compute SS, RMS
        for m = 1:N_tc
            total(m,:) = fit2(2*m+1) * sin(2*pi*Mx_freq(m)*tstar' - fit2(2*m));
        end
        model_sum = sum(total,1,'omitnan') + fit2(1);
        residuals2 = vel' - model_sum;

        T_fit_rms_all2(z,n) = sqrt(nanmean(residuals2.^2));
        ss_all2(z,n) = nansum(residuals2.^2);
    end
    end
end
%% ------------------------------------------------------------------------
%% export m4
[num_params2, n_bins, n_ens] = size(T_fit_result_all2);
results_flat2 = reshape(permute(T_fit_result_all2, [3, 2, 1]), [], num_params2);
% ens_list = repmat((1:n_ens)', n_bins, 1);%1–16, repeated 10x
% bin_list = repelem((1:n_bins)', n_ens); %1repeated16x,2repeated16x,etc
fit_table_M2M4 = array2table([ens_list, bin_list, results_flat2], ...
    'VariableNames', {'Ensemble', 'Bin', ...
                      'Residual', 'Phase_M2', 'Amp_M2', 'Phase_M4', 'Amp_M4'});
sheetName = date + "_M2M4";
base = 'C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\excel';
filename = 'tidefitresult.xlsx';
output_path = fullfile(base, filename);
writetable(fit_table_M2M4, output_path, 'Sheet', sheetName);
%% ------------------------------------------------------------------------
%% export m4 rms2
rms_flat_M2M4 = T_fit_rms_all2(:); 
rms_table_M2M4 = table(ens_list, bin_list, rms_flat_M2M4, ...
    'VariableNames', {'Ensemble', 'Bin', 'RMS_M2M4'});
sheetName = date + "_RMS_M2M4";
base = 'C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\excel';
filename = 'rms.xlsx';
output_path = fullfile(base, filename);
writetable(rms_table_M2M4, output_path, 'Sheet', sheetName);
%% ------------------------------------------------------------------------
%% Tc plot: choose ens and bin depth M2+M4
% ee = 7;  
% zz = 5;   
T_fit_result = T_fit_result_all2(:, zz, ee);
Pars = T_guess_all2(:, zz, ee);  % input guess

Tstar_re = squeeze(Tstar_clean(1, ee, :)); %extract
along_re = squeeze(along_sigma_clean(zz, ee, :));

% Tc_names = {'M2', 'M4'};
% Tc_periods = [12.4208, 6.2103];
% Mx_freq = 24*(LT_time(2)-LT_time(1)) ./ Tc_periods;
N_tc = length(Mx_freq);
tide_time = 0:0.01:1;

qp = NaN(N_tc, length(tide_time));     % fit
qp_in = NaN(N_tc, length(tide_time));  % input guess
figure('color','w')
for m = 1:N_tc
    subplot(N_tc + 1, 1, m)
    amp_fit = T_fit_result(2*m + 1); % extract fitted and input params
    phase_fit = T_fit_result(2*m);
    amp_in = Pars(2*m + 1);
    phase_in = Pars(2*m);
    qp(m,:) = amp_fit * sin(2 * pi * Mx_freq(m) * tide_time - phase_fit); % compute signals
    qp_in(m,:) = amp_in * sin(2 * pi * Mx_freq(m) * tide_time - phase_in);

    plot(tide_time, qp_in(m,:), 'Color', [0.5, 0.5, 0.5]); hold on
    plot(tide_time, qp(m,:), 'b')
    ylabel('Amp'); xlabel('T*'); title(Tc_names{m})
end

% composite prediction
subplot(N_tc + 1, 1, N_tc + 1)
input_sum = sum(qp_in, 1) + Pars(1);       % residual from guess
model_sum = sum(qp, 1) + T_fit_result(1);  % residual from fit

plot(tide_time, input_sum, 'Color', [0.5, 0.5, 0.5]); hold on
plot(tide_time, model_sum, 'b')
plot(Tstar_re, along_re, 'r*')

ylabel('Amplitude (m)'); xlabel('T*'); title('Sum of tidal constituents')
legend('Input parameters','Model Fit','Observed Data')
sgtitle(sprintf('Ensemble = %d, Bin = %d', ee, zz))

filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\constituents', ...
    sprintf('M2+M4_%s_line%d_ens%d_bin%d', date, line,ee,zz));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% pcolor sigma resid vel M2+M4
residuals2 = squeeze(T_fit_result_all2(1, :, :));  %10x16

figure('color','w')
pcolorjw(Xgrid, Zgrid, residuals2);  
% shading interp
    colormap(redblue)
    cb=colorbar; 
    ylabel(cb,'Residual Velocity (m/s)')
    caxis([-0.5,0.5])
xlabel('Distance (km)')
ylabel('Sigma Levels')
xlim([0,X(end)+0.025])
% xlim([0,max(X)])
title(strcat('Line=', num2str(line),'  Date=',(num2str(date))))
% 
filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit', ...
    sprintf('residual_velocity_M2+M4_%s_line%d_sigma', date, line));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% depth contour with resid vel M2M4
figure('color','w')
pcolor(X, Z_real, residuals2)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
    colormap(redblue)
    cb=colorbar
    ylabel(cb,'Residual Velocity (m/s)')
    caxis([-0.5,0.5])
xlim([0,X(end)+0.025])
xlabel('Distance (km)')
ylabel('Depth (m)')
title(strcat('Line=', num2str(line),'  Date=',(num2str(date))))
set(gca, 'YDir', 'reverse')  

filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\depth contour', ...
    sprintf('residudal_velocity_M2+M4_%s_line%d_depth', date, line));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% pcolor RMS M2+M4
figure('color','w')
pcolorjw(Xgrid, Zgrid, T_fit_rms_all2)
% shading flat
    cb=colorbar
    caxis([0,1])
xlim([0,X(end)+0.025])
xlabel('Distance (km)')
ylabel('Sigma Levels')
title(strcat('Line=', num2str(line),'  Date=',(num2str(date))))

filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\rmse', ...
    sprintf('RMSE_M2+M4_%s_line%d_depth', date, line));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% subplot DEPTH M2M4
phase=cmocean('phase');

figure('color','w')
ax(1)=subplot(2,3,1);
pcolor(Xgrid, Z_real, residuals); shading flat 
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(1),redblue)
    cb=colorbar; caxis([-0.5,0.5])
xlabel('Distance (km)')
ylabel(cb,'Residual Velocity (m/s)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
xlim([0,X(end)+0.025]) %line3
title('Residual Velocity')
sgtitle(strcat('Line=', num2str(line),'  Date=',(num2str(date))))

ax(2)=subplot(2,3,4);
pcolor(Xgrid, Z_real, T_fit_rms_all2) 
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(2),parula)
    cb=colorbar;
    caxis([0,1])
xlim([0,X(end)+0.025])
xlabel('Distance (km)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
title('RMSE')

m2amp = squeeze(T_fit_result_all2(3, :, :));  
m2phase = squeeze(T_fit_result_all2(2, :, :));  
    m2phase_deg= m2phase*180/pi;
m4amp = squeeze(T_fit_result_all2(5, :, :));  
m4phase = squeeze(T_fit_result_all2(4, :, :));  
    m4phase_deg= m4phase*180/pi;

ax(3)=subplot(2,3,2);
pcolor(Xgrid, Z_real, m2amp)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(3),redblue)
cb=colorbar;
ylabel(cb,'Amplitude (m)')
ma=max(max(m2amp));
caxis([0,ma])
xlim([0,X(end)+0.025])
xlabel('Distance (km)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
title('M2 Amplitude')

ax(4)=subplot(2,3,5);
pcolor(Xgrid, Z_real, m2phase_deg)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(4),phase)
cb=colorbar;
ylabel(cb,'Degrees')
caxis([0,360])
xlim([0,X(end)+0.025])
xlabel('Distance (km)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
title('M2 Phase')

ax(5)=subplot(2,3,3);
pcolor(Xgrid, Z_real, m4amp)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(5),redblue)
cb=colorbar;
caxis([0,ma * 1/2])
ylabel(cb,'Amplitude (m)')
xlim([0,X(end)+0.025])
xlabel('Distance (km)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
title('M4 Amplitude')

ax(6)=subplot(2,3,6);
pcolor(Xgrid, Z_real, m4phase_deg)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(6),phase)
cb=colorbar;
ylabel(cb,'Degrees')
caxis([0,360])
xlabel('Distance (km)')
xlim([0,X(end)+0.025])
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
title('M4 Phase')
% 
 filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit', ...
    sprintf('subplot_M2+M4_%s_line%d_sigma', date, line));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% 3rd Run: M2 + M4 + M6
Tc_names = {'M2','M4','M6'};
Tc_periods = [12.4208, 6.2103, 4.1402];
Mx_freq = 24*(LT_time(2)-LT_time(1)) ./ Tc_periods;
N_tc = length(Mx_freq);

% preallocate storage for third fit
T_fit_result_all3 = NaN(2*N_tc+1, n_bins, n_ens);  
T_fit_rms_all3 = NaN(n_bins, n_ens);
ss_all3 = NaN(n_bins, n_ens);
T_guess_all3 = NaN(2*N_tc + 1, n_bins, n_ens);

for z = 1:n_bins
    for n = 1:n_ens
        fit_M2M4 = T_fit_result_all2(:,z,n);

        vel = squeeze(along_sigma_clean(z,n,:));
        tstar = squeeze(Tstar_clean(1,n,:));

        % construct initial guess for M2 + M4 + M6
        Pars3 = NaN(2*N_tc+1, 1);
        Pars3(1) = fit_M2M4(1);  % residual
        Pars3(2) = fit_M2M4(2);  % M2 phase
        Pars3(3) = fit_M2M4(3);  % M2 amp
        Pars3(4) = fit_M2M4(4);  % M4 p
        Pars3(5) = fit_M2M4(5);  % M4 amp
        Pars3(6) = fit_M2M4(2)- (-0.24); % M6 frequency
        Pars3(7) = fit_M2M4(3) * 0.0125; % M6 amp

        if all(Pars3(1:3) == 0)
           Pars3(1:7) = NaN; 
        end

        T_guess_all3(:, z, n) = Pars3;

        r= find(isnan(vel));
            vel(r) = [];
            tstar(r) = [];
        rr= find(isnan(Pars3));
            Pars3(rr) = [];
            if isempty(Pars3)
                continue
            else

        fit3 = fminsearch(@adcp_tide_fit, Pars3, [], Mx_freq, tstar, vel); % third fit: M2 + M4 + M6 

            if fit3(3) < 0 % correct negative amplitudes
                fit3(3) = -fit3(3);
                fit3(2) = fit3(2) + pi;
            end
            if fit3(5) < 0
                fit3(5) = -fit3(5);
                fit3(4) = fit3(4) + pi;
            end
            if fit3(7) <0
                fit3(7) = -fit3(7);
                fit3(6) = fit3(6) + pi;
            end

        T_fit_result_all3(:,z,n) = fit3;

        total = NaN(N_tc, length(tstar)); % reconstruct model and compute SS, RMS
        for m = 1:N_tc
            total(m,:) = fit3(2*m+1) * sin(2*pi*Mx_freq(m)*tstar' - fit3(2*m));
        end
        model_sum = sum(total,1,'omitnan') + fit3(1);
        residuals3 = vel' - model_sum;

        T_fit_rms_all3(z,n) = sqrt(nanmean(residuals3.^2));
        ss_all3(z,n) = nansum(residuals3.^2);
    end
    end
end
%% ------------------------------------------------------------------------
%% Tc plot specify bin/ens M2+M4+M6
% ee = 7;  
% zz = 5;   

T_fit_result = T_fit_result_all3(:, zz, ee);
Pars = T_guess_all3(:, zz, ee);  % input guess

Tstar_re = squeeze(Tstar_clean(1, ee, :)); %extract
along_re = squeeze(along_sigma_clean(zz, ee, :));
valid = ~isnan(Tstar_re) & ~isnan(along_re);
Tstar_re = Tstar_re(valid);
along_re = along_re(valid);

% if isempty(Tstar_re) || any(isnan(T_fit_result))
%     disp('No valid data or fit result for this bin/ensemble.')
%     return
% end

N_tc = length(Mx_freq);
tide_time = 0:0.01:1;

qp = NaN(N_tc, length(tide_time));     % fit
qp_in = NaN(N_tc, length(tide_time));  % input guess

figure('color','w')

for m = 1:N_tc
    subplot(N_tc + 1, 1, m)

    amp_fit = T_fit_result(2*m + 1); % extract fitted and input params
    phase_fit = T_fit_result(2*m);

    amp_in = Pars(2*m + 1);
    phase_in = Pars(2*m);

    qp(m,:) = amp_fit * sin(2 * pi * Mx_freq(m) * tide_time - phase_fit); % compute signals
    qp_in(m,:) = amp_in * sin(2 * pi * Mx_freq(m) * tide_time - phase_in);

    plot(tide_time, qp_in(m,:), 'Color', [0.5, 0.5, 0.5]); hold on
    plot(tide_time, qp(m,:), 'b')
    ylabel('Amp'); xlabel('T*'); title(Tc_names{m})
end

% composite prediction
subplot(N_tc + 1, 1, N_tc + 1)
input_sum = sum(qp_in, 1, 'omitnan') + Pars(1);       % residual from guess
model_sum = sum(qp, 1, 'omitnan') + T_fit_result(1);  % residual from fit

plot(tide_time, input_sum, 'Color', [0.5, 0.5, 0.5]); hold on
plot(tide_time, model_sum, 'b')
plot(Tstar_re, along_re, 'r*')

ylabel('Velocity'); xlabel('T*'); title('Sum of tidal constituents')
legend('Input parameters','Model Fit','Observed Data')
sgtitle(sprintf('Ensemble = %d, Bin = %d', ee, zz))

filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\constituents', ...
    sprintf('M2+M4+M6_%s_line%d_ens%d_bin%d', date, line,ee,zz));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% export m2m4m6
[num_params3, n_bins, n_ens] = size(T_fit_result_all3);
results_flat3 = reshape(permute(T_fit_result_all3, [3, 2, 1]), [], num_params3);
% ens_list = repmat((1:n_ens)', n_bins, 1);%1–16, repeated 10x
% bin_list = repelem((1:n_bins)', n_ens); %1repeated16x,2repeated16x,etc
fit_table_M2M4M6 = array2table([ens_list, bin_list, results_flat3], ...
    'VariableNames', {'Ensemble', 'Bin', ...
                      'Residual', 'Phase_M2', 'Amp_M2', 'Phase_M4', 'Amp_M4','Phase_M6', 'Amp_M6'});
sheetName = date + "_M2M4M6";
base = 'C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\excel';
filename = 'tidefitresult.xlsx';
output_path = fullfile(base, filename);
writetable(fit_table_M2M4M6, output_path, 'Sheet', sheetName);
%% ------------------------------------------------------------------------
%% export m2m4m6 rms3
rms_flat_M2M4M6 = T_fit_rms_all3(:);
rms_table_M2M4M6 = table(ens_list, bin_list, rms_flat_M2M4M6, ...
    'VariableNames', {'Ensemble', 'Bin', 'RMS_M2M4M6'});
sheetName = date + "_RMS_M2M4M6";
base = 'C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\excel';
filename = 'rms.xlsx';
output_path = fullfile(base, filename);
writetable(rms_table_M2M4M6, output_path, 'Sheet', sheetName);
%% ------------------------------------------------------------------------
%% pcolor sigma resid vel M2+M4+M6
residuals3 = squeeze(T_fit_result_all3(1, :, :));  

figure('color','w')
pcolorjw(Xgrid, Zgrid, residuals3);  
% shading interp
    colormap(redblue)
    cb=colorbar; 
    ylabel(cb,'Residual Velocity (m/s)')
    caxis([-0.5,0.5])
xlabel('Distance (km)')
ylabel('Sigma Levels')
xlim([0,X(end)+0.025])
% xlim([0,max(X)])
title(strcat('Line=', num2str(line),'  Date=',(num2str(date))))
% 
filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit', ...
    sprintf('residual_velocity_M2+M4+M6_%s_line%d_sigma', date, line));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% depth contour with resid vel M2+M4+M6
figure('color','w')
pcolor(Xgrid, Z_real, residuals3)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
    colormap(redblue)
    cb=colorbar
    ylabel(cb,'Residual Velocity (m/s)')
    caxis([-0.5,0.5])
xlim([0,X(end)+0.025])
xlabel('Distance (km)')
ylabel('Depth (m)')
title(strcat('Line=', num2str(line),'  Date=',(num2str(date))))
set(gca, 'YDir', 'reverse')  

filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\depth contour', ...
    sprintf('residudal_velocity_M2+M4+M6_%s_line%d_depth', date, line));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% pcolor RMS M2M4M6 
figure('color','w')
pcolorjw(Xgrid, Zgrid, T_fit_rms_all3)
% shading flat
    cb=colorbar
    caxis([0,1])
xlim([0,X(end)+0.025])
xlabel('Distance (km)')
ylabel('Sigma Levels')
title(strcat('Line=', num2str(line),'  Date=',(num2str(date))))
% 
filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit\rmse', ...
    sprintf('RMSE_M2+M4+M6_%s_line%d_depth', date, line));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% subplot all deliverables M2M4M6
phase=cmocean('phase');
figure('color','w')
ax(1)=subplot(2,4,1);
pcolor(X, Z_real, residuals); shading flat 
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(1),redblue)
    cb=colorbar; caxis([-0.5,0.5])
xlabel('Distance (km)')
ylabel(cb,'Residual Velocity (m/s)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
xlim([0,X(end)+0.025]) %line3
title('Residual Velocity')
sgtitle(strcat('Line=', num2str(line),'  Date=',(num2str(date))))

ax(2)=subplot(2,4,5);
pcolor(X, Z_real, T_fit_rms_all2) 
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(2),parula)
    cb=colorbar;
    caxis([0,1])
xlim([0,X(end)+0.025])
xlabel('Distance (km)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
title('RMSE')

m2amp = squeeze(T_fit_result_all2(3, :, :));  %16x10
m2phase = squeeze(T_fit_result_all2(2, :, :));  %16x10
    m2phase_deg= m2phase*180/pi;
m4amp = squeeze(T_fit_result_all2(5, :, :));  %16x10
m4phase = squeeze(T_fit_result_all2(4, :, :));  %16x10
    m4phase_deg= m4phase*180/pi;
m6amp = squeeze(T_fit_result_all3(7, :, :));  %16x10
m6phase = squeeze(T_fit_result_all3(6, :, :));  %16x10
    m6phase_deg= m6phase*180/pi;

ax(3)=subplot(2,4,2);
pcolor(X, Z_real, m2amp)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(3),redblue)
cb=colorbar;
ylabel(cb,'Amplitude (m)')
ma=max(max(m2amp));
caxis([0,ma])
xlabel('Distance (km)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
title('M2 Amplitude')

ax(4)=subplot(2,4,6);
pcolor(X, Z_real, m2phase_deg)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(4),phase)
cb=colorbar;
ylabel(cb,'Degrees')
caxis([0,360])
xlabel('Distance (km)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
title('M2 Phase')

ax(5)=subplot(2,4,3);
pcolor(X, Z_real, m4amp)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(5),redblue)
cb=colorbar;
caxis([0,ma * 1/2])
ylabel(cb,'Amplitude (m)')
xlabel('Distance (km)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
title('M4 Amplitude')

ax(6)=subplot(2,4,7);
pcolor(X, Z_real, m4phase_deg)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(6),phase)
cb=colorbar;
ylabel(cb,'Degrees')
caxis([0,360])
xlabel('Distance (km)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
title('M4 Phase')

ax(7)=subplot(2,4,4);
pcolor(X, Z_real, m6amp)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(7),redblue)
cb=colorbar;
caxis([0,ma * 1/3])
ylabel(cb,'Amplitude (m)')
xlabel('Distance (km)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
title('M6 Amplitude')

ax(8)=subplot(2,4,8);
pcolor(X, Z_real, m6phase_deg)
shading flat
hold on
plot(X,h_mean,'k','LineWidth',2)
colormap(ax(8),phase)
cb=colorbar;
ylabel(cb,'Degrees')
caxis([0,360])
xlabel('Distance (km)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse') 
title('M6 Phase')

filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No7_bin_interp_sigma_tide_fit', ...
    sprintf('subplot_M2+M4+M6_%s_line%d_sigma', date, line));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);
%% ------------------------------------------------------------------------
%% save
save(strcat('BI_','adcp_','L',num2str(line),'_',...
    date,'_tide_','fit_','bin-interp_sigma','.mat'),'T_fit_result_all','T_fit_result_all2','T_fit_result_all3',...
    'Mx_freq','along_sigma_clean','h_bininterp_clean','time_bininterp_clean',...
    'Tc_names','Tc_periods','ss_all','ss_all2','ss_all3','h_mean','Z_real','Zgrid','Xgrid',...
     'T_fit_rms_all','T_fit_rms_all2','T_fit_rms_all3')
