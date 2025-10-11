%% No6_bin_interp_extrapolation.m
% NSE
% need to incorporate cross

% Creates extrapolated matrices from surface to bottom in 0.25m increments.
% Removes vel data exceeding water depth.
% Fills in measured data to extrapolated matrices and performs surface
% extrapolations with the following options:
% 1. constant, or
% 2. 0slope at surface where constant extrap applied to short profiles.

% Bottom extrapolation options:
% 1. linear using last good bin of data (some may require using second to last).
% 2. logarithmic using the deepest bin of data*
    % *there are a handful of lines that cannot be logarithmically extrapolated
    % using the last bin of data. In these instances, there are two options:
        % a. use second to last bin for logarithmic extrapolation
        % b. use the last bin of data for linear extrapolation

% Calculates instantaneous transport.
% Conversion to sigma coordinates.

% Note:  
% along_extrap = polynomial fit to 0slope at surface extrap and log bottom extrap.
% along_sigma = along_extrap mapped to sigma coordinates.
% transport = transport for along_extrap for each ens/bin/transect.
% transport_sum = sum of transport for a single value for each transect.
% 
% along_extrap_constant= constant surface extrap and bottom linear extrap.
% along_sigma_cl = along_extrap_constant mapped to sigma coordinates.
% Transport_cl = transport for along_extrap_constant for each ens/bin/transect.
% Transport_cl_sum = sum of transport_cl for a single value for each transect.

%% manually choose one survey:

line=3; %specify line #

date= "042023";
load BI_adcp_L3_042023_cropped_grid_rotate_bininterp.mat

% date= "062823";
% load BI_adcp_L3_062823_cropped_grid_rotate_bininterp.mat

% date= "071323";
% load BI_adcp_L3_071323_cropped_grid_rotate_bininterp.mat

% date= "031524";
% load BI_adcp_L3_031524_cropped_grid_rotate_bininterp.mat

% date= "031824";
% load BI_adcp_L3_031824_cropped_grid_rotate_bininterp.mat

%% set up matrices 
dz= 0.25;
z = 0:dz:33; %32.5 deepest
along_extrap = NaN(length(z),length(X),size(along_bininterp,3)); 
cross_extrap = NaN(length(z),length(X),size(along_bininterp,3)); 

% find indices of first and last bins with data
for xx = 1:size(along_bininterp,3)
    idx_start_interp(:,:,xx) = find(~isnan(h_bininterp(:,:,xx)),1,'first');
    idx_end_interp(:,:,xx) = find(~isnan(h_bininterp(:,:,xx)),1,'last');
end

%% remove data (velocity bins) that exceed the water depth 

for xx = 1:size(along_bininterp, 3)
    h_xx = h_bininterp(:,:,xx);
    along_xx = along_bininterp(:,:,xx);
    cross_xx = cross_bininterp(:,:,xx);
    for n = idx_start_interp(xx):idx_end_interp(xx)
        height_int = h_xx(n) - bins;
        height(:, n, xx) = height_int;
        neg_int = find(height_int <= 0);
        neg(:, n, xx) = height_int <= 0; 
        height_int(neg_int) = NaN;
        along_xx(neg_int, n) = NaN;
        cross_xx(neg_int, n) = NaN;
    end
    along_bininterp(:,:,xx) = along_xx;
    cross_bininterp(:,:,xx) = cross_xx;
end

%% continue setting up matrices
% cells arranged from shallowest to deepest; find first & last bins with good data

fgb_int = NaN(1,length(h_bininterp),size(along_bininterp,3));
lgb_int = NaN(1,length(h_bininterp),size(along_bininterp,3));

for xx = 1:size(along_bininterp, 3)
for n=idx_start_interp(xx):idx_end_interp(xx)
    a=find(~isnan(along_bininterp(:,n,xx)),1,'first'); 
    if isempty(a)
        fgb_int(:,n,xx)=NaN;
    else
        fgb_int(:,n,xx)=a;
    end
   
    b=find(~isnan(along_bininterp(:,n,xx)),1,'last');
    if isempty(b)
        lgb_int(:,n,xx)=NaN;
        continue 
    else
        lgb_int(:,n,xx)=b; 
    end
end
end

%% redefine indices 

for xx = 1:size(along_bininterp, 3)
idx_e= find(~isnan(lgb_int(:,:,xx)),1,'last');
idx_end_interp(:,:,xx)= idx_e; %last real bin of data after NaNs

idx_s=find(~isnan(fgb_int(:,:,xx)),1,'first');
idx_start_interp(:,:,xx)= idx_s;
end

%% Fill in middle using measured data
%incorporating 0.25m vertical increments

for xx = 1:size(along_bininterp, 3)
for n=idx_start_interp(xx):idx_end_interp(xx)
    if ~isnan(fgb_int(:,n,xx))        
    b = find(z >= bins(fgb_int(:,n,xx)) & z<= bins(lgb_int(:,n,xx))); 
    %find where depth > first good data in bins, and < last good data in bins
    along_mid = interp1(bins(fgb_int(:,n,xx):lgb_int(:,n,xx)), ...
        along_bininterp(fgb_int(:,n,xx):lgb_int(:,n,xx),n,xx),z(b));
    cross_mid = interp1(bins(fgb_int(:,n,xx):lgb_int(:,n,xx)), ...
        cross_bininterp(fgb_int(:,n,xx):lgb_int(:,n,xx),n,xx),z(b));
        along_extrap(b,n,xx) = along_mid;
        cross_extrap(b,n,xx) = cross_mid;
    end
end
end

along_mid = along_extrap;
cross_mid = cross_extrap;

%% PLOT bininterp and extrap bininterp at 25cm vertical:
% 
% for xx = 1:size(along_bininterp,3)
% figure('color', 'white')
% subplot(2,1,1)
% pcolorjw(X,bins,along_bininterp(:,:,xx))
% hold on
% plot(X,h_bininterp(:,:,xx),'k','LineWidth',2)
% shading flat; 
% cb=colorbar; ylabel(cb,'along channel velocity (m/s)')
% caxis([-1,1])
% ax=gca;
% set(gca,'xlim',[0,max(X)])
% ax.YDir='reverse';
% ax.YLim=[0,16];
% colormap(redblue)
% xlabel('Projected Distance (km)')
% ylabel('Depth')
% title('along (bin-interp), orig')
% 
% subplot(2,1,2)
% pcolorjw(X,z,along_extrap(:,:,xx))
% hold on
% plot(X,h_bininterp(:,:,xx),'k','LineWidth',2)
% shading flat; 
% cb=colorbar; ylabel(cb,'along channel velocity (m/s)')
% caxis([-1,1])
% ax=gca;
% set(gca,'xlim',[0,max(X)])
% ax.YDir='reverse';
% ax.YLim=[0,16];
% colormap(redblue)
% xlabel('Projected Distance (km)')
% ylabel('Depth')
% title('along (bin-interp), extrap (25cm increment depth)')
% sgtitle(strcat('Transect =', num2str(xx)));
% end

%% surface extrapolation: constant 
% copies the top most bin of data to surface, even short profiles

along_extrap_constant= along_extrap;
cross_extrap_constant= cross_extrap;

for xx = 1:size(along_bininterp, 3) 
for n=idx_start_interp(xx):idx_end_interp(xx)

d= find(~isnan(along_extrap_constant(:,n,xx)),1,'first');  
along_extrap_constant(1:d-1,n,xx)= along_extrap_constant(d,n,xx); 

e= find(~isnan(cross_extrap_constant(:,n,xx)),1,'first');  
cross_extrap_constant(1:e-1,n,xx)= cross_extrap_constant(e,n,xx);
end
end

%% ------------------------------------------------------------------------
%% testing surface extraps
% xx=10;n=3;
%starting index:
start_idx = fgb_int(:,n, xx); 
%skip short profiles
a = find(~isnan(along_bininterp(:,n,xx)), 1, 'first');  
fgb_int(:,n, xx) = a;
b = find(~isnan(along_bininterp(:,n,xx)), 1, 'last');
lgb_int(:,n, xx) = b;
%add more bins...from +3  7**11*****
x_data = bins(fgb_int(:,n,xx):fgb_int(:,n,xx)+3); 
v_data = along_bininterp(fgb_int(:,n,xx):fgb_int(:,n,xx)+3,n,xx); 
%% ------------------------------------------------------------------------
%% Choose one of the below options ----------------------------------------
%% ------------------------------------------------------------------------
%% try testing intial guesses with polyfit

p = polyfit(x_data, v_data, 2);  % p = [a, b, c]
x_surface = min(x_data);  % Or wherever your slope constraint applies
a0 = p(1);
c0 = p(3);
initial_guess = [a0, c0];
fun = @(coeff) sum((coeff(1)*x_data.^2 - 2*coeff(1)*x_surface*x_data + coeff(2) - v_data).^2, 'all');

coefficients = fminsearch(fun, initial_guess);


%% try this; it removes the b coeff. still looks like a straight line even though be is actually v small
x_surface = min(x_data);% Define surface location

% Custom cost function with b = -2a*x_surface
fun = @(coeff) sum((coeff(1)*x_data.^2 - 2*coeff(1)*x_surface*x_data + coeff(2) - v_data).^2, 'all');

initial_guess = [1, 0];  % [a, c]
coefficients = fminsearch(fun, initial_guess);

% Recover coefficients
% a = coef_opt(1);
% b = -2*a*x_surface;
% c = coef_opt(2);
% coefficients = [a, b, c];
% disp(coefficients)

a_v(:,n,xx) = coefficients(1);
b_v(:,n,xx) = -2*a_v(:,n,xx)*x_surface;
c_v(:,n,xx) = coefficients(2);
%calculate surface points
d = find(z <= bins(fgb_int(:,n,xx)));
z_surf = z(d);
v_surf = a_v(:,n,xx)*z_surf.^2 + b_v(:,n,xx)*z_surf + c_v(:,n,xx); %same # which is the c_v coefficient3


% fun = @(coeff) sum((coeff(1) * x_data.^2 + coeff(2) * x_data + coeff(3) - v_data).^2, 'all');
% initial_guess = [1, 1, 1];
% 
% x_surface = min(x_data);  % Or 0, or wherever your "surface" is
% Aeq = [2*x_surface, 1, 0];  % Enforce slope = 0 at x_surface
% Beq = 0;
% 
% options = optimset('fmincon');
% coefficients = fmincon(fun, initial_guess, [], [], Aeq, Beq, [], [], [], options);
% 
% %%
% 
% x_fit = linspace(min(x_data), max(x_data), 100);
% v_fit = coefficients(1)*x_fit.^2 + coefficients(2)*x_fit + coefficients(3);
% 
% plot(x_data, v_data, 'bo'); hold on
% plot(x_fit, v_fit, 'r-'); hold off
% legend('Data', 'Fitted Polynomial');



%% 1. Constrained to 0 with optimization ----------------------------------
fun = @(coeff) sum((coeff(1) * x_data.^2 + coeff(2) * x_data + coeff(3) - v_data).^2,'all');
initial_guess = [1, 1, 1]; % Initial guess for coefficients [a, b, c]
% initial_guess = [0, 0, mean(v_data)]; %same result

% Constraint: b (coeff(2)) should be 0
Aeq = [0, 1, 0]; % Coefficient matrix for linear inequality constraint;enforce coeff(2) = 0
Beq = 0; % Right-hand side of the constraint; must be 0

% Options for optimization
options = optimset('fmincon');
% options.Display = 'iter' % Display optimization progress

% Perform constrained optimization
coefficients = NaN(1,3);
coefficients = fmincon(fun, initial_guess, [], [], Aeq, Beq, [], [], [], options);

%   FMINCON attempts to solve problems of the form:
%    min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%     X                     C(X) <= 0, Ceq(X) = 0   (nonlinear constraints)
%                              LB <= X <= UB  
%   X = FMINCON(FUN,X0,A,B,Aeq,Beq) minimizes FUN subject to the linear 
%   equalities Aeq*X = Beq as well as A*X <= B. (Set A=[] and B=[] if no 
%   inequalities exist.)
%% ------------------------------------------------------------------------
%% 2. Unconstrained no optimization ---------------------------------------
coefficients = polyfit(x_data(:), v_data, 2);  % no slope constraint
%% ------------------------------------------------------------------------
%% 3. testing soft > hard constraints w optimization ----------------------
initial_guess = [0, 0, mean(v_data)];
fun = @(coeff) sum((coeff(1)*x_data.^2 + coeff(2)*x_data + coeff(3) - v_data).^2, 'all') ...
               + 1e5 * coeff(2)^2;  % 1e3 same as no constraint, try 1e5, SAME
options = optimset('fmincon');
coefficients = fmincon(fun, initial_guess, [], [], [], [], [], [], [], options);
%% ------------------------------------------------------------------------
%% extract coefficients ---------------------------------------------------
a_v(:,n,xx) = coefficients(1);
b_v(:,n,xx) = coefficients(2);
c_v(:,n,xx) = coefficients(3);
%calculate surface points
d = find(z <= bins(fgb_int(:,n,xx)));
z_surf = z(d);
v_surf = a_v(:,n,xx)*z_surf.^2 + b_v(:,n,xx)*z_surf + c_v(:,n,xx); %same # which is the c_v coefficient3
%% ------------------------------------------------------------------------
%% plots to test profiles
figure('color','w')
plot(v_data,x_data,'r.-')
hold on
plot(v_surf,z_surf,'bo-')
set(gca, 'YDir', 'reverse')
%%

plot(v_surf,z_surf,'c*-')
legend('data','constrained: zero slope','unconstrained')

title('transect xx=10 profile n=3')

xlim([-0.85,-0.1])
ylim([0,4])
%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
%% surface extrapolation: Fit parabola to shallowest meter * need to account for cross vel from here on
% four data points, requiring zero slope at surface using optimization 
% see parabola_fit_example_zero_slope.m
% constant extrap for short profiles
% xx=10;n=3;

for xx = 1:size(along_bininterp, 3) 
for n=idx_start_interp(xx):idx_end_interp(xx)

    %starting index:
    start_idx = fgb_int(:, n, xx); 

    %check if at least 5 data point avail: *test 9 *test 13
        if start_idx + 4 > numel(bins) %4 %8 %12
            fprintf('Profile n=%d in transect xx=%d skipped: not enough data points.\n', n, xx);
            continue
        end
        if isnan(start_idx) 
            fprintf('Profile n=%d in transect xx=%d skipped: invalid start index (%g).\n', n, xx, start_idx);
            continue;
        end

    %skip short profiles
        a = find(~isnan(along_bininterp(:,n,xx)), 1, 'first');  
        if isempty(a)
            fgb_int(:, n, xx) = NaN;
        else
            fgb_int(:, n, xx) = a;
        end

        b = find(~isnan(along_bininterp(:,n,xx)), 1, 'last');
        if isempty(b)
            lgb_int(:, n, xx) = NaN;
        elseif b <= 5   %5 *13****
            fprintf('Profile n=%d in transect xx=%d skipped from extrapolation (not enough bins).\n', n, xx);
            continue  
        else
            lgb_int(:, n, xx) = b;
        end
% add more bins...from +3  7**11*****
    x_data = bins(fgb_int(:,n,xx):fgb_int(:,n,xx)+3); 
    v_data = along_bininterp(fgb_int(:,n,xx):fgb_int(:,n,xx)+3,n,xx); 
    % u_data = cross_bininterp(fgb_int(:,n,xx):fgb_int(:,n,xx)+3,n,xx); 

    % Define the optimization problem
    fun = @(coeff) sum((coeff(1) * x_data.^2 + coeff(2) * x_data + coeff(3) - v_data).^2,'all');
    initial_guess = [1, 1, 1]; % Initial guess for coefficients [a, b, c]
    % initial_guess = [0, 0, mean(v_data)]; %same result

    % Constraint: b (coeff(2)) should be 0
    Aeq = [0, 1, 0]; % Coefficient matrix for linear inequality constraint;enforce coeff(2) = 0
    Beq = 0; % Right-hand side of the constraint; must be 0

    % Options for optimization
    options = optimset('fmincon');
    % options.Display = 'iter' % Display optimization progress

    % Perform constrained optimization
    coefficients = NaN(1,3);
    coefficients = fmincon(fun, initial_guess, [], [], Aeq, Beq, [], [], [], options);

    % extract coefficients
    a_v(:,n,xx) = coefficients(1);
    b_v(:,n,xx) = coefficients(2);
    c_v(:,n,xx) = coefficients(3);

    %calculate surface points
    d = find(z <= bins(fgb_int(:,n,xx)));
    z_surf = z(d);
    v_surf = a_v(:,n,xx) * z_surf.^2 + b_v(:,n,xx) * z_surf + c_v(:,n,xx); %same # which is the c_v coefficient3

    %fill in to extrap matrices
    along_extrap(d,n,xx)=v_surf; 
end
end

% constant extrap for short profiles:
for xx = 1:size(along_bininterp, 3) 
for n=idx_start_interp(xx):idx_end_interp(xx)

    % starting index:
    start_idx = fgb_int(:, n, xx); 
        if isnan(start_idx) 
            fprintf('Profile n=%d in transect xx=%d skipped: invalid start index (%g).\n', n, xx, start_idx);
            continue;
        end
    % skip long profiles
        a = find(~isnan(along_bininterp(:,n,xx)), 1, 'first');  
        if isempty(a)
            fgb_int(:, n, xx) = NaN;
        else
            fgb_int(:, n, xx) = a;
        end
        b = find(~isnan(along_bininterp(:,n,xx)), 1, 'last');
        if isempty(b)
            lgb_int(:, n, xx) = NaN;
        elseif b > 5 %***********  check to see if this should be 4
            %fprintf('Profile n=%d in transect xx=%d skipped from CONSTANT extrapolation: log used\n', n, xx);
            continue  
        else
            fprintf('Profile n=%d in transect xx=%d constant extrap applied\n', n, xx);
            lgb_int(:, n, xx) = b;
        end
            d= find(~isnan(along_extrap(:,n,xx)),1,'first');  
            along_extrap(1:d-1,n,xx)= along_extrap(d,n,xx);
            % d= find(~isnan(cross_extrap(:,n,xx)),1,'first');  
            % cross_extrap(1:d-1,n,xx)= cross_extrap(d,n,xx);
end
end

%% plots to test profiles
figure('color','w')
plot(v_data,x_data,'r.-')
hold on
plot(v_surf,z_surf,'bo-')
set(gca, 'YDir', 'reverse')


plot(v_surf,z_surf,'m*-')
legend('data','constrained: zero slope','unconstrained')

title('transect xx=10 profile n=3')

xlim([-0.85,-0.1])
ylim([0,4])

%%
disp(coefficients)
    % 0.0000         0    1.2830
    % 0.1287   -0.5466   -0.1125



%% PLOT the difference between two surface extraps
% 
% dv= along_extrap_constant - along_extrap;
% test_con= along_extrap_constant(:,:,1);
% test_opt= along_extrap(:,:,1);
% dv1=dv(:,:,1);
% 
% dv_opt_ext = along_extrap - along_extrap_constant;
% dvoptext1=dv_opt_ext(:,:,1);
% 
% % % figure
% xx=11; 
% % % Select 10 index numbers for vertical profiles
% profiles = 1:10;
% % profiles = 9:18;
% % profiles = 20:29;
% % profiles = 22;
% % profiles = 30:39;
% % profiles = 40:49;
% % profiles = 50:59;
% % profiles = 60:69;
% % profiles = 70:79;
% % profiles = 80:89;
% % profiles = 90:99;
% % profiles = 100:109;
% 
% figure('color','w')
% for p = 1:10
%     subplot(2,5,p)
%     plot(dv_opt_ext(:,profiles(p),xx),z,'ko:')
%     hold on
%     % plot(along_bininterp(:,profiles(p),xx),bins,'b*-')
%     % plot(0,h_bininterp(:,profiles(p),xx),'ks','markerfacecolor','k')
%     set(gca,'YDir','reverse','XAxisLocation','top')
%     ylabel('Depth (m)')
%     % xlabel('Alongchannel velocity (m/s)')
%     title(strcat('profile = ',num2str(profiles(p))))
% end


%% PLOTS x sec comparing log and lin surface extraps

% for xx = 1:size(along_bininterp,3)
% figure('color', 'white')
% subplot(2,1,1)
% pcolor(X,z,along_extrap_constant(:,:,xx))
% hold on
% plot(X,h_bininterp(:,:,xx),'k','LineWidth',2)
% shading flat; 
% cb=colorbar; ylabel(cb,'along channel velocity (m/s)')
% caxis([-1.5,1.6])
% ax=gca;
% set(gca,'xlim',[0,max(X)])
% ax.YDir='reverse';
% ax.YLim=[0,16];
% colormap(redblue)
% xlabel('Projected Distance (km)')
% ylabel('Depth')
% title('constant extrap')
% 
% subplot(2,1,2)
% pcolor(X,z,along_extrap(:,:,xx))
% hold on
% plot(X,h_bininterp(:,:,xx),'k','LineWidth',2)
% shading flat; 
% cb=colorbar; ylabel(cb,'along channel velocity (m/s)')
% caxis([-1.5,1.5])
% ax=gca;
% set(gca,'xlim',[0,max(X)])
% ax.YDir='reverse';
% ax.YLim=[0,16];
% colormap(redblue)
% xlabel('Projected Distance (km)')
% ylabel('Depth')
% title('4 pt optimization extrap / constant')
% sgtitle(strcat('Transect =', num2str(xx)));
% end

%% bottom extrap: linear

z=z';
for xx = 1:size(along_bininterp, 3) 
for n=idx_start_interp(xx):idx_end_interp(xx)

    a= find(~isnan(along_extrap_constant(:,n,xx)),1,'last');  
    b= find(z <= h_bininterp(:,n,xx),1,'last');

    if b>a
        test= interp1([z(a,1) h_bininterp(:,n,xx)],[along_extrap_constant(a,n,xx) 0],[z(a+1:b)]);
        along_extrap_constant(a+1:b,n,xx)= test(:);
    else if b==a 
         fprintf('Profile n=%d in transect xx=%d skipped; used second to last bin for linear extrap to bottom\n', n, xx);
         a2= find(~isnan(along_extrap_constant(:,n,xx)),1,'last')-1; %minus one bin
        test2= interp1([z(a2,1) h_bininterp(:,n,xx)],[along_extrap_constant(a2,n,xx) 0],[z(a2+1:b)]);
        along_extrap_constant(a2+1:b,n,xx)= test2(:);
      end
   end
    
end
end

%% Bottom extrapolation: Fit logarithmic profile to bottom 
% requires curve-fitting toolbox
nearbottom = 0.08;
%make a copy to compare
along_extrap_copy = along_extrap;
%CHECK UNITS OF roughness length = 0.067 (meters?)
fit_params = fittype('a*log(x/0.067)', 'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'a'});

for xx = 1:size(along_bininterp, 3)
    for n = idx_start_interp(xx):idx_end_interp(xx)
        height2 = h_bininterp(:,n,xx) - bins;
        b = find(height2<=2 & height2>=0.1);%select cells btwn this threshold (m above bottom).

        if all(isnan(along_extrap(:, n, xx)), 'all')
            continue;%if true (all nans) then with skip rest of loop
        end

        % select bins for log fitting
        if b(end) >= lgb_int(:,n,xx) && lgb_int(:,n,xx) > 2
            x2 = height2(lgb_int(:,n,xx)-2 : lgb_int(:,n,xx));
            v2 = along_bininterp(lgb_int(:,n,xx)-2 : lgb_int(:,n,xx),n,xx);
            % u2 = cross_bininterp(lgb_int(:,n,xx)-2 : lgb_int(:,n,xx),n,xx);
        else
            x2 = height2(b);
            v2 = along_bininterp(b,n,xx);
            % u2 = cross_bininterp(b,n,xx);
        end

        %log fit
        try
            if all(~isnan(v2)) %&& all(x2 > 0)
                coeff_v = fit(x2(:), v2(:), fit_params, 'start', 0);
               % coeff_u = fit(x2(:),u2(:),fit_params,'start',0);
                fit_method = "log";
            else
                fit_method = "linear";
            end
        catch
            warning('Log fit failed at profile n=%d in transect xx=%d; falling back to linear.', xx, n);
            fit_method = "linear";
        end

        if strcmp(fit_method, "log")
            %first attempt: log extrapolation using last bin of data
            d2 = find(z>=bins(lgb_int(:,n,xx)) & z<=h_bininterp(:,n,xx)-nearbottom);
            height_extrap = h_bininterp(:,n,xx)-z;

            if ~isempty(d2)
                h = height_extrap(d2);
                v_bott = coeff_v.a .* log(h ./ 0.067);
                % u_bott = coeff_v.a .* log(h ./ 0.067);
                along_extrap(d2, n, xx) = v_bott; %fill in to extrap matrices
                % cross_extrap(d2, n, xx) = u_bott;
                continue; % done

            else
                %second attempt: log extrap using second to last bin of data
                fprintf('Profile n=%d in transect xx=%d used second to last bin LOG.\n', n, xx);
                d3 = find(z>=bins(lgb_int(:,n,xx)-1) & z<h_bininterp(:,n,xx)-nearbottom);
                if ~isempty(d3)
                    h = height_extrap(d3);
                    v_bott = coeff_v.a .* log(h ./ 0.067);
                    % u_bott = coeff_v.a .* log(h ./ 0.067);
                    along_extrap(d3, n, xx) = v_bott;%fill in to extrap matrices
                    % cross_extrap(d3, n, xx) = u_bott;
                    continue; % done
                else
                    %fallback to linear extrap
                    fit_method = "linear";
                end
            end
        end

        if strcmp(fit_method, "linear")
            a = find(~isnan(along_extrap(:,n,xx)), 1, 'last');
            b_idx = find(z<=h_bininterp(:,n,xx), 1, 'last');

            if ~isempty(a) && ~isempty(b_idx)
                if b_idx > a
                    test = interp1([z(a,1), h_bininterp(:,n,xx)],...
                                   [along_extrap(a,n,xx),0],...
                                   z(a+1:b_idx));
                    along_extrap(a+1:b_idx,n,xx) = test(:);
                    fprintf('Profile n=%d in transect xx=%d used LINEAR fallback extrapolation.\n', n, xx);
                elseif b_idx == a
                    fprintf('Profile n=%d in transect xx=%d b==a; used second to last bin linear extrap.\n', n, xx);
                    a2 = find(~isnan(along_extrap_constant(:,n,xx)),1,'last')-1;%minus one bin
                    % if a2 >= 1
                        test2 = interp1([z(a2,1), h_bininterp(:,n,xx)],...
                                        [along_extrap_constant(a2,n,xx),0],...
                                        z(a2+1:b_idx));
                        along_extrap(a2+1:b_idx,n,xx) = test2(:);
                    % end
                end
            end
        end
    end
end

%% PLOT select profiles
%1) constant (surface) and linear (bottom) extrapolation
%2) polynomial- 0slope (surface) and logarithmic (bottom) extrapolation
% 
xx=10;
% % Select 10 index numbers for vertical profiles
profiles = 2:11;
% profiles = 10:16;
% % profiles = 20:29;
% % profiles = 22;
% % profiles = 30:39;
% % profiles = 40:49;
% % profiles = 50:59;
% % profiles = 60:69;
% % profiles = 70:79;
% % profiles = 80:89;
% % profiles = 90:99;
% % profiles = 100:109;
% 
figure('color','w')
for p = 1:10
    subplot(2,5,p)
    plot(along_extrap(:,profiles(p),xx),z,'mo:')
    hold on
    plot(along_extrap_constant(:,profiles(p),xx),z,'g.:')
    plot(along_mid(:,profiles(p),xx),z,'b.:')
    % plot(along_bininterp(:,profiles(p),xx),bins,'b*-')
    plot(0,h_bininterp(:,profiles(p),xx),'ks','markerfacecolor','k')
    set(gca,'YDir','reverse','XAxisLocation','top')
    ylabel('Depth (m)')
    xlabel('Alongchannel velocity (m/s)')
    title(strcat('profile = ',num2str(profiles(p))))
    legend('poly/log','constant/lin','data','H','Location','best')
end
sgtitle(strcat('Transect =', num2str(xx)));

%% PLOT just along and saves to folder

for xx = 1:size(along_bininterp,3)
% xx=4;
figure('color', 'white')
pcolorjw(X,z,along_extrap(:,:,xx))
hold on
plot(X,h_bininterp(:,:,xx),'k','LineWidth',2)
shading flat; 
cb=colorbar; ylabel(cb,'along channel velocity (m/s)')
caxis([-1.5,1.5])
set(gca,'ydir','reverse','fontsize',12)
ax=gca;
set(gca,'xlim',[0,max(X)]);%0.45]);%0.36])%0.59]) %l2
ax.YDir='reverse';
ax.YLim=[0,16];
colormap(redblue)
xlabel('Projected Distance (km)')
ylabel('Depth (m)')
title(strcat('Date= ', date, ' Transect=', num2str(xx)));

filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No6_bin_interp_extrap', ...
    sprintf('Along_%s_line%d_transect%d', date, line, xx));
export_fig([filename, '.png'], '-m2');
savefig([filename, '.fig']);

end

%% add instantaneous transports
% along extrap (poly/log)
% X is l (25m)
l = (X(2) - X(1)) * 1000; %km to m
w = dz/100; %25cm to m
area = l*w;

transport = NaN(size(along_extrap));
for xx = 1:size(along_extrap, 3) 
    for nn = 1:length(X)
        for zz = 1:length(along_extrap)
            transport(zz,nn,xx) = along_extrap(zz,nn,xx) * area;
        end
    end
end

for xx = 1:size(along_extrap, 3) 
transport_sum(xx) = nansum(nansum(transport(:,:,xx)));
%sum cols then rows
end

% along extrap constant linear
transport_cl = NaN(size(along_extrap_constant));
for xx = 1:size(along_extrap_constant, 3) 
    for nn = 1:length(X)
        for zz = 1:length(along_extrap_constant)
            transport_cl(zz,nn,xx) = along_extrap_constant(zz,nn,xx) * area;
        end
    end
end

for xx = 1:size(along_extrap_constant, 3) 
transport_sum_cl(xx) = nansum(nansum(transport_cl(:,:,xx)));
%sum cols then rows
end

%% Convert depth to sigma coordinates
N_levels = 10; %number of vertical levels
sigma_levels = linspace(1,0,N_levels);
along_sigma = NaN(N_levels,width(along_extrap),size(along_bininterp,3));
cross_sigma = NaN(N_levels,width(cross_extrap),size(along_bininterp,3));

for xx= 1:size(along_bininterp,3)
for n=idx_start_interp(xx):idx_end_interp(xx)
    sig = (h_bininterp(:,n,xx)-z)/h_bininterp(:,n,xx);
    if all(isnan(sig))
        continue; %if true (all nans) then with skip rest of loop
    else
    d= find(sig>nearbottom); 
    vel_v = interp1([sig(d);0],[along_extrap(d,n,xx);0],sigma_levels,'linear');%,'extrap');
    vel_v(end)=0; %seafloor
    % vel_u = interp1([sig(d);0],[cross_extrap(d,n,xx);0],sigma_levels,'linear');
    % vel_u(end)=0; %seafloor
    % vel_v(1)=along_extrap(1,n,xx); %surface
    % vel_u(1)=cross_extrap(1,n,xx); 

    along_sigma(:,n,xx)=vel_v(:);
    % cross_sigma(:,n,xx)=vel_u(:);

% % Plots ALL profiles ( A LOT OF FIGURES )
    % figure('color','w')
    % plot(along_sigma(:,n,xx),sigma_levels,'ro-')
    % hold on
    % plot(along_extrap(d,n,xx),sig(d),'k*-')
    % set(gca,'XAxisLocation','top')
    % ylabel('Depth (m)')
    % xlabel('Alongchannel velocity (m/s)')
    %     title(strcat('profile = ',num2str(n),' transect =', num2str(xx)));
    %         legend('extrap','sigma','Location','best')

end
end
end
%% single profile fig

% n=7;
% xx=2;
% 
% figure('color','w')
% plot(along_sigma(:,n,xx),sigma_levels,'ro-')
% hold on
%     sig = (h_bininterp(:,n,xx)-z)/h_bininterp(:,n,xx); 
%     d= find(sig>0); 
% plot(along_extrap(d,n,xx),sig(d),'k*-')
% set(gca,'XAxisLocation','top')
% ylabel('Depth (m)')
% xlabel('Alongchannel velocity (m/s)')
% title(strcat('profile = ',num2str(n),' transect =', num2str(xx)));

%% PLOTS: Sigma x sec figs

for xx = 1:size(along_bininterp,3)
figure('color','w')
subplot(2,1,1)
pcolorjw(X,z,along_extrap(:,:,xx))
hold on
plot(X,h_bininterp(:,:,xx),'k','linewidth',2)
shading flat
cb = colorbar; ylabel(cb,'along-channel velocity (m/s)')
caxis([-1.5,1.5])
ax=gca;
set(gca,'YDir','reverse','YLim',[0,16])
set(gca,'xlim',[0,max(X)])
colormap(redblue)
xlabel('Distance (km)')
ylabel('Depth (m)')
title('Final Velocity')

subplot(2,1,2) 
pcolorjw(X,sigma_levels,along_sigma(:,:,xx))
cb = colorbar; ylabel(cb,'along-channel velocity (m/s)')
colormap(redblue)
caxis([-1.5,1.5])
ax=gca;
set(gca,'xlim',[0,max(X)])
xlabel('Distance (km)')
ylabel('Depth (m)')
shading flat
title('Sigma Coords')
% set(gca, 'YDir', 'reverse') 
sgtitle(strcat('Transect =', num2str(xx)));
% 
% filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No6_bin_interp_extrap\sigma', ...
%     sprintf('Along_Sigma_%s_line%d_transect%d', date, line, xx));
% export_fig([filename, '.png'], '-m2');
% savefig([filename, '.fig']);

end

%% Convert depth to sigma coordinates CONSTANT LIN EXTRAP
N_levels = 10; %number of vertical levels
sigma_levels = linspace(1,0,N_levels);
along_sigma_cl = NaN(N_levels,width(along_extrap_constant),size(along_bininterp,3));
cross_sigma_cl = NaN(N_levels,width(cross_extrap_constant),size(along_bininterp,3));

for xx= 1:size(along_bininterp,3)
for n=idx_start_interp(xx):idx_end_interp(xx)
    sig = (h_bininterp(:,n,xx)-z)/h_bininterp(:,n,xx);
    if all(isnan(sig))
        continue; %if true (all nans) then with skip rest of loop
    else
    d= find(sig>nearbottom); 
    vel_v = interp1([sig(d);0],[along_extrap(d,n,xx);0],sigma_levels,'linear');%,'extrap');
    vel_v(end)=0; %seafloor
    % vel_u = interp1([sig(d);0],[cross_extrap(d,n,xx);0],sigma_levels,'linear');
    % vel_u(end)=0; %seafloor
    % vel_v(1)=along_extrap_constant(1,n,xx); %surface
    % vel_u(1)=cross_extrap_constant(1,n,xx); 

    along_sigma_cl(:,n,xx)=vel_v(:);
    % cross_sigma_cl(:,n,xx)=vel_u(:);
end
end
end

%% PLOTS: Sigma constant linear x sec figs

for xx = 1:size(along_bininterp,3)
figure('color','w')
subplot(2,1,1)
pcolorjw(X,z,along_extrap_constant(:,:,xx))
hold on
plot(X,h_bininterp(:,:,xx),'k','linewidth',2)
shading flat
cb = colorbar; ylabel(cb,'along-channel velocity (m/s)')
caxis([-1.5,1.5])
ax=gca;
set(gca,'YDir','reverse','YLim',[0,16])
set(gca,'xlim',[0,max(X)])
colormap(redblue)
xlabel('Distance (km)')
ylabel('Depth (m)')
title('Final Velocity')

subplot(2,1,2) 
pcolorjw(X,sigma_levels,along_sigma_constant(:,:,xx))
cb = colorbar; ylabel(cb,'along-channel velocity (m/s)')
colormap(redblue)
caxis([-1.5,1.5])
ax=gca;
set(gca,'xlim',[0,max(X)])
xlabel('Distance (km)')
ylabel('Depth (m)')
shading flat
title('Sigma Coords')
% set(gca, 'YDir', 'reverse') 
sgtitle(strcat('Transect =', num2str(xx)));
% 
% filename = fullfile('C:\Users\nsert\Documents\MATLAB\CZM\2023_Surveys\concatenated data\No6_bin_interp_extrap\sigma\cl', ...
%     sprintf('Along_Sigma_%s_line%d_transect%d', date, line, xx));
% export_fig([filename, '.png'], '-m2');
% savefig([filename, '.fig']);

end

%% SAVE variables
save(strcat('BI_','adcp_','L',num2str(line),'_',...
    date,'_cropped_','grid_','rotate_','bininterp_','extrap','.mat'),'along_extrap',...
    'along_extrap_constant','z','height_extrap','idx_start_interp','idx_end_interp',...
    'along_sigma','along_sigma_cl','sigma_levels','X','time_bininterp','h_bininterp',...
    'transport','transport_sum','transport_cl','transport_sum_cl')
