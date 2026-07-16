clc; clear; close all;

addpath(pwd+"/src", pwd+"/data")
load(pwd+"/data/modellingData.mat")
noLags = 30;

%% Plot modelling data
clc; close all;

figure; sgtitle("Modelling data")
subplot(311), plot(pow_m),   title("Power")
subplot(312), plot(atemp_m), title("Ambient air temp")
subplot(313), plot(wtemp_m), title("Supply water temp")

%%
plotACFnPACF(atemp_m, noLags, "atemp"); set(gcf,'Position',[200 200 800 400])
plotACFnPACF(wtemp_m, noLags, "wtemp"); set(gcf,'Position',[200 200 800 400])
plotACFnPACF(pow_m, noLags, "power");   set(gcf,'Position',[200 200 800 400])
 
%% ARMA Air temperature
clc; close all;

% Select which coefficients to estimate
Am_params = [1 2 23];
Cm_params = [1];

% Create empty parameter arrays
Am = [1 zeros(1,30)];
Cm = [1 zeros(1,30)];
% Set selected coefficients at their indexes to 1
Am(Am_params+1) = 1;
Cm(Cm_params+1) = 1;

% Add Differential operator
s = 24;
s_poly = [1 zeros(1,s-1) -1];
Am = conv(Am,s_poly);

[modellatemp,ey,~,~] = estimateARMA(atemp_m, Am, Cm, "",noLags);
set(gcf,'Position',[200 200 800 400])
figure; whitenessTest(ey); close;

%% ARMA Power
clc; close all; 

Am_params = [1 2];
Cm_params = [24];

Am = [1 zeros(1,30)];
Cm = [1 zeros(1,30)];
Am(Am_params+1) = 1;
Cm(Cm_params+1) = 1;

s = 12;
s_poly = [1 zeros(1,s-1) -1];
Am = conv(Am,s_poly);

[modellpow,ey,~,~] = estimateARMA(atemp_m, Am, Cm, "",noLags);
set(gcf,'Position',[200 200 800 400])
figure; whitenessTest(ey); close;

%% ARMA Water temperature
clc; close all; 

Am_params = 1;
Cm_params = 1;

Am = [1 zeros(1,30)];
Cm = [1 zeros(1,30)];
Am(Am_params+1) = 1;
Cm(Cm_params+1) = 1;

[modellwtemp,ey,~,~] = estimateARMA(wtemp_m, Am, Cm, "",noLags);
set(gcf,'Position',[200 200 800 400])
figure; whitenessTest(ey); close;

%% Save models to /data directory
clc; close all;

dataLocation = pwd+"/data";
save(fullfile(dataLocation, 'modellatemp.mat'),'modellatemp');
save(fullfile(dataLocation, 'modellpow.mat'),'modellpow');
save(fullfile(dataLocation, 'modellwtemp.mat'),'modellwtemp');
