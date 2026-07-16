clc; clear; close all;

addpath(pwd+"/src", pwd+"/data")

% Load data
load('modellingData.mat');
load('validationData.mat');
load('testData.mat');
load('testData2.mat');

% Load trained models
inputModel_a = load("modellatemp.mat").modellatemp;
inputModel_w = load("modellwtemp.mat").modellwtemp;
outputModel = load("modellpow.mat").modellpow;
BJmodel_air = load("BJmodel_airinput.mat").foundModel_airinput;

%% Plot modelling data
clc; close all;

figure; sgtitle("Modelling data")
subplot(311), plot(pow_m),   title("Power")
subplot(312), plot(atemp_m), title("Ambient air temp")
subplot(313), plot(wtemp_m), title("Supply water temp")

%% k-step prediction
clc; close all;

foundModel = BJmodel_air;
inputModel = inputModel_a;

% -------- ---- Choose k ------- ----%

k = 7;

% ---- Set input and output data ----%

data_set = 1;

% ---- ------------------------- ----%

% Pre-differentiate trend?
% y = filter([1 -1], 1, y);
% y = y(2:end);
% x = x(2:end);

% ---- ------------------------- ----%

switch data_set
	case 1
		x = atemp_m;
		y = pow_m;
	case 2
		x = atemp_v;
		y = pow_v;
	case 3
		x = atemp_t;
		y = pow_t;
	case 4
		x = atemp_t2;
		y = pow_t2;
end

% Form polynomials
KA = conv( foundModel.D, foundModel.F );
KB = conv( foundModel.D, foundModel.B );
KC = conv( foundModel.F, foundModel.C );

[Fx, Gx] = polydivision( inputModel.C, inputModel.A, k );
[Fy, Gy] = polydivision( foundModel.C, foundModel.D, k );
[Fhh, Ghh] = polydivision( conv(Fy, KB), KC, k );

% Form input predictions
xhatk = filter(Gx, inputModel.C, x);

% Form output predictions
yhatk  = filter(Fhh, 1, xhatk) + filter(Ghh, KC, x) + filter(Gy, KC, y);

% Residuals
remove = max([length(Fhh), length(Ghh), length(Gy)]);
indV = remove:length(y);
ey = y(indV)-yhatk(indV);

% Naive predictor
[yNaive, var_naive,~] = naivePred(y, indV, k);

% Visuals

% Show only the last 200 samples
xindV = [length(y(indV))-200 length(y(indV))];

figure('Position', [100 100 800 400])
plot([y(indV) yhatk(indV) yNaive] )
legend('Output signal', 'Predicted output', 'Naive prediction')
title( sprintf('Predicted output signal, y_{t+%i|t}', k) )
xlim(xindV);

figure('Position', [100 100 800 400])
plot([x(indV) xhatk(indV)] )
legend('Input signal', 'Predicted input')
title( sprintf('Predicted input signal, x_{t+%i|t}', k) )
xlim(xindV);

figure('Position', [100 100 800 400])
str = append('y_{t+',int2str(k),'|t}');
acf( yhatk(indV), 30, 0.05, 1 );
title( sprintf('ACF (%s)',str))

% Print variances
fprintf("Residual variance k="+k+": " + var(ey) + ", normalized variance="+var(ey)/var(y(indV))+"\n");
fprintf("Residual variance naive : " + var_naive + ", normalized variance="+var_naive/var(y(indV))+"\n");
score = var(ey)/var_naive;
fprintf('sigma2e/sigma2naive : ')
score

% For discussion
% Was it necessary to add water temperature?
switch data_set
	case 1
		x_test = wtemp_m;
	case 2
		x_test = wtemp_v;
	case 3
		x_test = wtemp_t;
	case 4
		x_test = wtemp_t2;
end
myxcorr(ey,x_test(indV),30);
set(gcf,'Position',[200 200 800 400]);
title('X-corr of residuals and water temperature')
ylim([-0.15 0.15]);
