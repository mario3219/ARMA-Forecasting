clc; clear; close all;

addpath(pwd+"/src", pwd+"/data")
load(pwd+"/data/projectData25.mat")

% Load data
pow   = data(:,2); % power
atemp = data(:,3); % air temperature
wtemp = data(:,4); % water temperature

% Columns:% 1: Observations
% 5: Years
% 6: Months
% 7: Days
% 8: Hours

%% Analyze raw data
clc; close all;

figure; sgtitle("Raw data")
subplot(311), plot(pow),   title("Power")
subplot(312), plot(atemp), title("Ambient air temp")
subplot(313), plot(wtemp), title("Supply water temp")

days   = data(:,7);
figure; plot(days);      title("Days");
figure; plot(data(:,6)); title("Months");
figure; sgtitle("Missing data spike")
subplot(411); plot(days);  title("Days");               xlim([3650 3830]);
subplot(412), plot(pow),   title("Power");              xlim([3650 3830]);
subplot(413), plot(atemp), title("Ambient air temp");   xlim([3650 3830]);
subplot(414), plot(wtemp), title("Supply water temp");  xlim([3650 3830]);

figure; plot(data(:,6)); title("Months at the missing data"); xlim([3650 3830]);

arr = unique(data(:,5));
disp(["Year span:" num2str(arr')])   

% For the days plot:
% At idx 3756 there is missing data
% Jumps from 3 hours into day 10, then to
% the last 12 hours in day 28
% In total:
% (28-10)*24+14+12=458 points of missing data
% But it all occured during the same month during december, so might not be
% statistically significant, raw data looks fine...

% Data has alot of outliers. Remove the worst using a median filter

%% Median filter
clc; close all;

windowSize = 5;    % Choose median filter window size
wtemp  = myMedianFilter(wtemp, windowSize, 1, "Water temperature");
atemp  = myMedianFilter(atemp, windowSize, 1, "Air temperature");
pow    = myMedianFilter(pow,   windowSize, 1, "Power");

%% Create datasets
clc; close all;

model_size  = 7*24*7;
val_size    = 3*24*7;             % validation data size (in weeks)
test_size   = 3*24*7;             % 1st test data size (in weeks)

model_start = 1000;
model_stop  = model_start+model_size;
val_start   = model_stop+1;
val_stop    = val_start+val_size;
test_start  = val_stop+1+24*7;
test_stop   = test_start+test_size;
test2_start = 5000-test_size;
test2_stop  = 5000;

modelInt    = model_start:model_stop;
valInt      = val_start:val_stop;
testInt     = test_start:test_stop;
testInt2    = test2_start:test2_stop;

% Modelling data
pow_m   = pow(modelInt);    % power
atemp_m = atemp(modelInt);  % air temperature
wtemp_m = wtemp(modelInt);  % water temperature

% Validation data
pow_v   = pow(valInt);      % power
atemp_v = atemp(valInt);    % air temperature
wtemp_v = wtemp(valInt);    % water temperature

% 1st Test data
pow_t   = pow(testInt);     % power
atemp_t = atemp(testInt);   % air temperature
wtemp_t = wtemp(testInt);   % water temperature

% 2nd Test data
pow_t2   = pow(testInt2);   % power
atemp_t2 = atemp(testInt2); % air temperature
wtemp_t2 = wtemp(testInt2); % water temperature

fprintf("Modelling set:   "+model_start+":"+model_stop+"\n")
fprintf("Validation set:  "+val_start+":"+val_stop+"\n")
fprintf("1st test set:    "+test_start+":"+test_stop+"\n")
fprintf("2nd test set:    "+test2_start+":"+test2_stop+"\n")

%% Plot datasets
clc; close all;

noLags = 30;

% Plot modelling data
figure('Position', [100 100 800 800])
sgtitle("Modelling data")
subplot(311), plot(pow_m),   title("Power")
subplot(312), plot(atemp_m), title("Ambient air temp")
subplot(313), plot(wtemp_m), title("Supply water temp")

% Plot validation data
figure('Position', [100 100 800 800])
sgtitle("Validation data")
subplot(311), plot(pow_v),   title("Power")
subplot(312), plot(atemp_v), title("Ambient air temp")
subplot(313), plot(wtemp_v), title("Supply water temp")

% Plot test data
figure('Position', [100 100 800 800])
sgtitle("1st Test data")
subplot(311), plot(pow_t),   title("Power")
subplot(312), plot(atemp_t), title("Ambient air temp")
subplot(313), plot(wtemp_t), title("Supply water temp")

figure('Position', [100 100 800 800])
sgtitle("2nd Test data")
subplot(311), plot(pow_t2),   title("Power")
subplot(312), plot(atemp_t2), title("Ambient air temp")
subplot(313), plot(wtemp_t2), title("Supply water temp")

% ACF of data set outputs
figure; acf(pow_v,30,0.05,1); set(gcf,'Position',[200 200 800 400]);
title('ACF (Power, validation)')
figure; acf(pow_t,30,0.05,1); set(gcf,'Position',[200 200 800 400]);
title('ACF (Power, test 1)')
figure; acf(pow_t2,30,0.05,1); set(gcf,'Position',[200 200 800 400]);
title('ACF (Power, test 2)')

%% Modelling set ACF's and PACF's
% See if there is season or trends in the modelling data
clc; close all;

noLags = 60;
figure('Position', [100 100 800 800]);
subplot(313); acf(wtemp_m, noLags, 0.05, 1); title("Water temperature");
subplot(311); acf(atemp_m, noLags, 0.05, 1); title("Air temperature");
subplot(312); acf(pow_m, noLags, 0.05, 1); title("Power");

% Season of 12 in power
% Season of 24 in air temperature
% No season in water temperature
% No trends

%% Remove outliers in water temperature data in the modelling set
clc; close all;

% Copy the array
wtemp_out = wtemp_m;

% Choose outlier intervals
%outliers = [485:563 774:787 796:808 980:1040 1060:1129];      % start=900
outliers = [384:460 881:937];                     % start=1000

% Set outlier intervals to NaN
wtemp_out(outliers) = NaN;

% Interpolate Nan values
wtemp_out = fillmissing(wtemp_out, 'linear');

% Plot new water temperature data in the modelling set
figure; sgtitle("Water Temperature with Outliers Removed");
subplot(211); plot(wtemp_m), title("Supply Water Temp (Before)"); ylim([min(wtemp_m) max(wtemp_m)]);
subplot(212); plot(wtemp_out), title("Supply Water Temp (After)");ylim([min(wtemp_m) max(wtemp_m)]);

wtemp_m = wtemp_out;

%% Does the data need transformation?
clc; close all;

figure; sgtitle("Box-cox")
subplot(311); l = bcNormPlot(pow_m,1); title("Power")
xline(max(l), 'r', sprintf('\\lambda = %.2f', max(l)), ...
      'LabelVerticalAlignment', 'bottom', ...
      'LabelHorizontalAlignment', 'center');

subplot(312); l = bcNormPlot(wtemp_m,1); title("Water temperature")
xline(max(l), 'r', sprintf('\\lambda = %.2f', max(l)), ...
      'LabelVerticalAlignment', 'bottom', ...
      'LabelHorizontalAlignment', 'center');
subplot(313); l = bcNormPlot(atemp_m,1); title("Air temperature")
xline(max(l), 'r', sprintf('\\lambda = %.2f', max(l)), ...
      'LabelVerticalAlignment', 'bottom', ...
      'LabelHorizontalAlignment', 'center');

% Plots suggest transformation, but data looks visually fine, gonna
% try not transforming the data for now

%% Save datasets to /data directory
clc; close all;

dataLocation = pwd+"/data";
save(fullfile(dataLocation, 'modellingData.mat'),'pow_m','atemp_m','wtemp_m');
save(fullfile(dataLocation, 'validationData.mat'),'pow_v','atemp_v','wtemp_v');
save(fullfile(dataLocation, 'testData.mat'),'pow_t','atemp_t','wtemp_t');
save(fullfile(dataLocation, 'testData2.mat'),'pow_t2','atemp_t2','wtemp_t2');
