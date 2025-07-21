% No7_bin_avg_depth_avg_tide_fit.m
% NSE
% Tide fitting for horizontally bin averaged, depth averaged velocity data

%% SELECT LINE MANUALLY
line=3;
%% ------------------------------------------------------------------------
%% CHANGE THE L# in MAT FILES BELOW:
%% 04-20-23
% date = "042023";
% load('BI_adcp_L3_042023_cropped_grid_rotate_binavg_extrap.mat')
%% ------------------------------------------------------------------------
%% 06-28-23
date = "062823";
load('BI_adcp_L3_062823_cropped_grid_rotate_binavg_extrap.mat')
%% ------------------------------------------------------------------------
%% 07-13-23
% date = "071323";
% load('BI_adcp_L3_071323_cropped_grid_rotate_binavg_extrap.mat')
%% ------------------------------------------------------------------------
%% 03-15-24
% date = "031524";
% load('BI_adcp_L3_031524_cropped_grid_rotate_binavg_extrap.mat')
%% ------------------------------------------------------------------------
%% 03-18-24
% date = "031824";
% load('BI_adcp_L3_031824_cropped_grid_rotate_binavg_extrap.mat')
%% ------------------------------------------------------------------------
%% average over ensembles (depth-avg)

along_avg = zeros(1,length(X), size(along_extrap,3));

for xx = 1:size(along_extrap, 3)
   for n = 1:length(X)
 along_avg(1,n,xx)= nanmean(along_extrap(:,n,xx));
   end
end
%% ------------------------------------------------------------------------
%% plots
%test: and prep: for loop graph time and horiz vel transect number assign to time
% tstar and depth avg vel and see if fit ok
% 
% for xx = 1:size(along_extrap,3)
% figure('color', 'white')
%     subplot(2,1,1)
% plot(time_binavg(1,:,xx),along_avg(1,:,xx),'k*-')
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
% for xx = 1:size(along_extrap,3)
% plot(X,along_avg(1,:,xx))
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
%% Convert time to fraction of tide cycle (Tstar)

Tstar=zeros(size(time_binavg));

for xx = 1:size(along_extrap, 3)
    for n=1:length(X)
        prior= find(LT_time<=time_binavg(1,n,xx),1,'last');
        Tstar(1,n,xx) = (time_binavg(1,n,xx)-LT_time(prior))/(LT_time(prior+1)-LT_time(prior));
    end
end

%Turn 0s back into NaNs
Tstar(Tstar == 0) = NaN;

Tstar1= Tstar(:,:,1);

% reshape
Tstar_re = squeeze(Tstar)';
%% ------------------------------------------------------------------------
%% create a minimum criteria for transects
% if entire column has 6 or less values, then turn all to NaNs

% Need to loop through ensemble #s:
% Initial values for amplitude and phase
along_avg_re = squeeze(along_avg)'; %reshapes to 16x18;(# of transects x # of ensembles)
h_binavg_re = squeeze(h_binavg)';
time_binavg_re = squeeze(time_binavg)';

   for n=1:length(X)
       if sum(~isnan(along_avg_re(:,n))) <= 6  %count amount of non-NaNs, if <=4 turn to NaNs
          along_avg_re(:,n) = NaN;
          Tstar_re(:,n) = NaN;
          h_binavg_re(:,n) = NaN;
          time_binavg_re(:,n) = NaN;
       end
   end

along_avg_avg = nanmean(along_avg_re); %1x18;mean value (1 x #ensembles)
h_binavg_avg = nanmean(h_binavg_re);
time_binavg_avg = nanmean(time_binavg_re);
% Tstar_re_avg = nanmean(Tstar_re_copy);
%% ------------------------------------------------------------------------
%% 1st Run for M2 Tc 
%Set up tidal constituent parameters
Tc_names = [{'M2'}];
Tc_periods = [12.4208]; 
Mx_freq = 24*(LT_time(2)-LT_time(1))./Tc_periods; %per tidal cycle (calculated by elapsed time between low tides divided by the tidal period of each harmonic (M2=12.4208 hours, M4=6.2103 hours, M6=4.1402 hours)

    for n=1:length(X)
     % for xx= 1:size(along_extrap,3)
        Pars(1,n) = along_avg_avg(n); %residual velocity (guessing average) 
        % M2 
        Pars(2,n) = 0; %M2 phase shift (guessing zero- when ran just M2 to get output)
        Pars(3,n) = (max(along_avg_re(:,n))-min(along_avg_re(:,n)))/2; %M2 amplitude (guessing range/2)
     % end
        if isnan(Pars(1,n)) 
          Pars(2,n) = NaN; % turn 0 to NaNs in phasing
        end
    end
%% ------------------------------------------------------------------------
%% copies

%make copies which contain filled matrices (non copies will be overwritten
%later when extracting vectors)
Pars_copy = Pars;
Tstar_re_copy = Tstar_re;
along_avg_re_copy = along_avg_re;
%% ------------------------------------------------------------------------
%% Run minimization with fminsearch: t fit

T_fit_result_mtx = NaN(size(Pars_copy,1), length(X));
T_fit_rms = NaN(size(X));
 
for n=1:length(X) % 18 ensembles
    % for n=3;
    Pars = Pars_copy(:,n); %isolates single column
    Tstar_re = Tstar_re_copy(:,n);
    along_avg_re= along_avg_re_copy(:,n);
        r= find(isnan(along_avg_re));
            along_avg_re(r) = [];
            Tstar_re(r) = [];
        rr= find(isnan(Pars));
            Pars(rr) = [];
            if isempty(Pars)
                continue
            else
    T_fit_result = fminsearch(@adcp_tide_fit,Pars,[],Mx_freq,Tstar_re,along_avg_re);

a = find(T_fit_result(3,:) < 0); %negative amplitudes
T_fit_result(3,a) = T_fit_result(3,a) * -1; %make +, and
T_fit_result(2,a) = T_fit_result(2,a) + pi; %shift phase
    
    T_fit_result_mtx(:,n) = T_fit_result; %fill in final matrix

T_fit_rms(:,n)= adcp_tide_fit(T_fit_result,Mx_freq,Tstar_re,along_avg_re);
            end
end
%% ------------------------------------------------------------------------
%% export m2
% table = array2table(T_fit_result_mtx);
% table.Properties.VariableNames = cellstr(string(1:width(table)));
% sheetName = date + "M2"; 
% writetable(table, 'tidefitresult.xlsx', 'Sheet', sheetName);

%% ------------------------------------------------------------------------
 %% sum of squares

ss_mtx = NaN(1, length(X));

for n=1:length(X) % 18 ensembles
    T_fit_result = T_fit_result_mtx(:,n);
    Tstar_re = Tstar_re_copy(:,n);
    along_avg_re = along_avg_re_copy(:,n);

    N_tc=length(Mx_freq);
    tide_time=0:0.01:1;

        for m=1:N_tc
        total(m,:)=T_fit_result(2*m+1,1).*sin(2*Mx_freq(m)*pi.*Tstar_re-T_fit_result(2*m,1));
        end %creates matrix where each row is another constituent at Tstar

    total_sum= sum(total,1)+T_fit_result(1); %sum of predicted tide at each Tstar (adding the avg of sine curve)
    along_avg_re= along_avg_re'; %make row vector
    d = along_avg_re - total_sum;
    ss = nansum(d.^2);
    along_avg_re= along_avg_re';
    
    ss_mtx(:,n) = ss;
end
%% ------------------------------------------------------------------------
%% export m2 ss
% table = array2table(ss_mtx);
% table.Properties.VariableNames = cellstr(string(1:width(table)));
% sheetName = date + "ss" + "M2"; 
% writetable(table, 'ss.xlsx', 'Sheet', sheetName);
%% ------------------------------------------------------------------------
%% Tc plot single ens/bin M2
ee=14;
Pars = Pars_copy(:,ee);
T_fit_result = T_fit_result_mtx(:,ee);
Tstar_re = Tstar_re_copy(:,ee);
along_avg_re = along_avg_re_copy(:,ee);

N_tc=length(Mx_freq);
tide_time=0:0.01:1;

figure('color','w')
for m=1:N_tc
    subplot(N_tc+1,1,m)
    qp_in(m,:)=Pars(2*m+1,1).*sin(2*Mx_freq(m)*pi.*tide_time-Pars(2*m,1));
    qp(m,:)=T_fit_result(2*m+1,1).*sin(2*Mx_freq(m)*pi.*tide_time-T_fit_result(2*m,1));
    plot(tide_time,qp_in(m,:),'color',[0.5,0.5,0.5]); hold on
    plot(tide_time,qp(m,:),'b')
    ylabel('Amplitude (m)'); xlabel('T*'); title(Tc_names{m})
end
subplot(N_tc+1,1,m+1)
plot(tide_time,sum(qp_in,1)+Pars(1),'color',[0.5,0.5,0.5]); hold on
plot(tide_time,sum(qp,1)+T_fit_result(1),'b')
hold on
plot(Tstar_re,along_avg_re,'r*')
ylabel('Amplitude (m)'); xlabel('T*'); title('Sum of tidal constituents')
legend('Input parameters','fminsearch result','observed data')
sgtitle(strcat('Ensemble=', num2str(ee)));
%% ------------------------------------------------------------------------
%% 2nd run
Tc_names = [{'M2'} {'M4'}]; 
Tc_periods = [12.4208, 6.2103];
Mx_freq = 24*(LT_time(2)-LT_time(1))./Tc_periods; %per tidal cycle (calculated by elapsed time between low tides divided by the tidal period of each harmonic (M2=12.4208 hours, M4=6.2103 hours, M6=4.1402 hours)

for n=1:length(X)
    Pars2(1,n) = T_fit_result_mtx(1,n);
    %M2
    Pars2(2,n) = T_fit_result_mtx(2,n); %output from M2 phase
    Pars2(3,n) = T_fit_result_mtx(3,n); %M2 amplitude
    % M4
    Pars2(4,n) = T_fit_result_mtx(2,n)-(-0.16888); %what is this value?
    Pars2(5,n) = T_fit_result_mtx(3,n)*0.0125; %""?
        if all(Pars2(1:3,n))==0
          Pars2(1:5,n) = NaN; % turn 0 to NaNs in phasing
        end
end
%% ------------------------------------------------------------------------
%% copies 2
Pars2_copy = Pars2;
%% ------------------------------------------------------------------------
%% Run minimization with fminsearch: 2nd run M2+M4

T_fit_result_mtx2 = NaN(size(Pars2_copy,1), length(X));
T_fit_rms2 = NaN(size(X));

for n=1:length(X) % 18 ensembles
    % for n=3;
    Pars2 = Pars2_copy(:,n); %isolates single column
    Tstar_re = Tstar_re_copy(:,n);
    along_avg_re= along_avg_re_copy(:,n);
        r= find(isnan(along_avg_re));
            along_avg_re(r) = [];
            Tstar_re(r) = [];
        rr= find(isnan(Pars2));
            Pars2(rr) = [];
            if isempty(Pars2)
                continue
            else
    T_fit_result2 = fminsearch(@adcp_tide_fit,Pars2,[],Mx_freq,Tstar_re,along_avg_re);
    
a = find(T_fit_result2(3,:) < 0); %negative amplitudes
T_fit_result2(3,a) = T_fit_result2(3,a) * -1; %make +, and
T_fit_result2(2,a) = T_fit_result2(2,a) + pi; %shift phase    
% later will add to this: 180-b this is an old note...?
b = find(T_fit_result2(5,:) < 0);
T_fit_result2(5,b) = T_fit_result2(5,b) * -1;
T_fit_result2(4,b) = T_fit_result2(4,b) + pi;

    T_fit_result_mtx2(:,n) = T_fit_result2; %fill in final matrix

    T_fit_rms2(:,n)= adcp_tide_fit(T_fit_result2,Mx_freq,Tstar_re,along_avg_re);
            end
end
%% ------------------------------------------------------------------------
%% export m4
% table2 = array2table(T_fit_result_mtx2);
% table2.Properties.VariableNames = cellstr(string(1:width(table2)));
% sheetName = date + "M2" + "M4"; 
% writetable(table2, 'tidefitresult.xlsx', 'Sheet', sheetName);

%% sum of squares M2+M4

ss2_mtx = NaN(1, length(X));

for n=1:length(X) % 18 ensembles
    T_fit_result2 = T_fit_result_mtx2(:,n);
    Tstar_re = Tstar_re_copy(:,n);
    along_avg_re = along_avg_re_copy(:,n);

    N_tc=length(Mx_freq);
    tide_time=0:0.01:1;

        for m=1:N_tc
            total(m,:)=T_fit_result2(2*m+1,1).*sin(2*Mx_freq(m)*pi.*Tstar_re-T_fit_result2(2*m,1));
        end %creates matrix where each row is another constituent at Tstar

    total_sum2= sum(total,1)+T_fit_result2(1); %sum of predicted tide at each Tstar (adding the avg of sine curve)
    along_avg_re= along_avg_re'; %make row vector
    d = along_avg_re - total_sum2;
    ss2 = nansum(d.^2);
    along_avg_re= along_avg_re';
    
    ss2_mtx(:,n) = ss2;
end
%% ------------------------------------------------------------------------
%% export m2m4 ss
% table = array2table(ss2_mtx);
% table.Properties.VariableNames = cellstr(string(1:width(table)));
% sheetName = date + "ss2" + "M2" + "M4"; 
% writetable(table, 'ss.xlsx', 'Sheet', sheetName);
%% ------------------------------------------------------------------------
%% plot
%avg vel vs tstar

% for xx = 1:size(along_extrap,2)
% figure('color', 'white')
% plot(Tstar_re_copy(:,xx),along_avg_re_copy(:,xx),'k*')
% xlabel('tstar')
% ylabel('avg velocity')
% title(strcat('Ensemble =', num2str(xx)));
% end
%% ------------------------------------------------------------------------
%% Tc plot single ens/bin M2+M4
ee=7;
Pars = Pars2_copy(:,ee);
T_fit_result = T_fit_result_mtx2(:,ee);
Tstar_re = Tstar_re_copy(:,ee);
along_avg_re = along_avg_re_copy(:,ee);

N_tc=length(Mx_freq);
tide_time=0:0.01:1;

figure('color','w')
for m=1:N_tc
    subplot(N_tc+1,1,m)
    qp_in(m,:)=Pars(2*m+1,1).*sin(2*Mx_freq(m)*pi.*tide_time-Pars(2*m,1));
    qp(m,:)=T_fit_result(2*m+1,1).*sin(2*Mx_freq(m)*pi.*tide_time-T_fit_result(2*m,1));
    plot(tide_time,qp_in(m,:),'color',[0.5,0.5,0.5]); hold on
    plot(tide_time,qp(m,:),'b')
    ylabel('Amplitude (m)'); xlabel('T*'); title(Tc_names{m})
end
subplot(N_tc+1,1,m+1)
plot(tide_time,sum(qp_in,1)+Pars(1),'color',[0.5,0.5,0.5]); hold on
plot(tide_time,sum(qp,1)+T_fit_result(1),'b')
hold on
plot(Tstar_re,along_avg_re,'r*')
ylabel('Amplitude (m)'); xlabel('T*'); title('Sum of tidal constituents')
legend('Input parameters','fminsearch result','observed data')
sgtitle(strcat('Ensemble=', num2str(ee)));

%% save
% save(strcat('BI_','adcp_','L',num2str(line),'_',...
%     date,'_tide_','fit','.mat'),'T_fit_result_mtx','T_fit_result_mtx2',...
%     'Tstar_re_copy','along_avg_re_copy','Pars_copy','Pars2_copy','Mx_freq',...
%     'Tc_names','Tc_periods','ss_mtx','ss2_mtx','h_binavg_avg','time_binavg_avg',...
%      'T_fit_rms','T_fit_rms2')