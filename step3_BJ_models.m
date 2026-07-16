clc; clear; close all;

addpath(pwd+"/src", pwd+"/data")

load('modellingData.mat');
load('validationData.mat');
load('testData.mat');

modelpow   = load("modellpow.mat").modellpow;
modelatemp = load("modellatemp.mat").modellatemp;
modelwtemp = load("modellwtemp.mat").modellwtemp;

present(modelpow), present(modelatemp), present(modelwtemp)

noLags=30;

%% Plot modelling data
clc; close all;

figure; sgtitle("Modelling data")
subplot(311), plot(pow_m),   title("Power")
subplot(312), plot(atemp_m), title("Ambient air temp")
subplot(313), plot(wtemp_m), title("Supply water temp")

%% Estimate B and A2 polynomial orders using air temperature as input
clc; close all;

% Cross-correlate input and output residuals with the found model
inputModel = modelatemp;
xM = atemp_m;
yM = pow_m;

ex = filter( inputModel.A, inputModel.C, xM ); ex = ex(length(inputModel.A)+30:end );
ey = filter( inputModel.A, inputModel.C, yM ); ey = ey(length(inputModel.A)+30:end );

figure('Position', [100 100 800 400]);
[Cxy,lags] = xcorr( ey, ex, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(ey) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between filtered in- and output')

%% Initial BJ model, estimating C1 and A1
clc; close all;

d=0;
s=2;
r=0;

A2      = [ones(1,r+1)];
B       = [zeros(1,d) ones(1,s+1)];

[model,ey,~,~] = estimateBJ( yM, xM, [], [], B, A2, '', noLags );
set(gcf,'Position',[200 200 800 400])
remove = length(model.B);
tilde_xt = xM(remove:end );

figure('Position', [100 100 800 400]);
[Cxy,lags] = xcorr( tilde_xt, ey, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(yM) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between input and residual without the influence of the input')
xlim([0 30]), ylim([-1 1])

%% Final BJ model using air temperature as input
clc; close all;

A1_params = [1 2];
C1_params = [1 24];

A1      = [1 zeros(1,30)]; A1(A1_params+1) = 1;
C1      = [1 zeros(1,30)]; C1(C1_params+1) = 1;

s = 12;
s_poly = [1 zeros(1,s-1) -1];
A1 = conv(A1,s_poly);

[foundModel_airinput,ey,~,~] = estimateBJ( yM, xM, C1, A1, B, A2, '', noLags);
set(gcf,'Position',[200 200 800 400])
remove = length(foundModel_airinput.B);

figure('Position', [100 100 800 400]);
subplot(211)
[Cxy,lags] = xcorr( ey, xM(remove:end), noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(ey) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between input and ey residuals')
ylim([-1 1]); xlim([0 30])

subplot(212)
[Cxy,lags] = xcorr( ey, yM(remove:end), noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(ey) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between output and ey residuals')
ylim([-1 1]); xlim([0 30])

%% Estimate B and A2 polynomial orders using water temperature as input
clc; close all;

% Cross-correlate input and output residuals with the found model
inputModel = modelwtemp;
xM = wtemp_m;
yM = pow_m;

ex = filter( inputModel.A, inputModel.C, xM ); ex = ex(length(inputModel.A)+30:end );
ey = filter( inputModel.A, inputModel.C, yM ); ey = ey(length(inputModel.A)+30:end );

figure;
[Cxy,lags] = xcorr( ey, ex, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(ey) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between filtered in- and output')

%% Initial BJ model, estimating C1 and A1
clc; close all;

d=0;
s=0;
r=0;

A2      = [ones(1,r+1)];
B       = [zeros(1,d) ones(1,s+1)];

[model,ey,~,~] = estimateBJ( yM, xM, [], [], B, A2, '', noLags );
remove = length(model.B);
tilde_xt = xM(remove:end );

figure
[Cxy,lags] = xcorr( tilde_xt, ey, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(yM) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between input and residual without the influence of the input')
xlim([0 30]), ylim([-1 1])

%% Final BJ model using water temperature as input
clc; close all;

A1_params = [1 2];
C1_params = [1 24];

A1      = [1 zeros(1,30)]; A1(A1_params+1) = 1;
C1      = [1 zeros(1,30)]; C1(C1_params+1) = 1;

s = 12;
s_poly = [1 zeros(1,s-1) -1];
A1 = conv(A1,s_poly);

[foundModel_waterinput,ey,~,~] = estimateBJ( yM, xM, C1, A1, B, A2, '', noLags);
remove = length(model.B);
tilde_xt = xM(remove:end );
tilde_yt = yM(remove:end );

figure;
[Cxy,lags] = xcorr( ey, tilde_xt, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(ey) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between input and xt residuals')
ylim([-1 1]); xlim([0 30])
fprintf("Whiteness for CCF (input and ehat)\n")
figure; whitenessTest(Cxy); close;

figure;
[Cxy,lags] = xcorr( ey, tilde_yt, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(ey) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between output and xt residuals')
ylim([-1 1]); xlim([0 30])
fprintf("Whiteness for CCF (Output and ehat)\n")
figure; whitenessTest(Cxy); close;

%% Save BJ model
clc;

dataLocation = pwd+"/data";
save(fullfile(dataLocation, 'BJmodel_airinput.mat'),'foundModel_airinput');
save(fullfile(dataLocation, 'BJmodel_waterinput.mat'),'foundModel_waterinput');