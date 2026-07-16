clc; clear; close all;

addpath(pwd+"/src", pwd+"/data")

% Load data
load('modellingData.mat');
load('validationData.mat');
load('testData.mat');
load('testData2.mat');

% Load trained models
inputModel1 = load("modellatemp.mat").modellatemp;
inputModel2 = load("modellwtemp.mat").modellwtemp;
outputModel = load("modellpow.mat").modellpow;
BJmodel_air = load("BJmodel_airinput.mat").foundModel_airinput;
BJmodel_water = load("BJmodel_waterinput.mat").foundModel_waterinput;

present(BJmodel_air); present(BJmodel_water);
noLags = 30;

% Use modelling data during estimation
y = pow_m;
x1 = atemp_m;
x2 = wtemp_m;

% Note that the function idpoly use the notation: 
%
%       A(z) y(t) = [B(z)/F(z)] u(t) + [C(z)/D(z)] e(t)
%
% This means that:
%
%   A(z) = 1,       B(z) = B(z),    F(z) = A2(z)
%   C(z) = C1(z),   D(z) = A1(z)

%% Re-estimate parameters and noise order
clc; close all;

% Extract polynomials from models, and set non-zero values to 1

B_1     = BJmodel_air.B;    B_1 = B_1 ~= 0;
A2_1    = BJmodel_air.F;    A2_1 = A2_1 ~= 0;

B_2     = BJmodel_water.B;  B_2 = B_2 ~= 0; 
A2_2    = BJmodel_water.F;  A2_2 = A2_2 ~= 0;

B{1,1}  = B_1;  % B_1;
A2{1,1} = A2_1; % A2_1;

B{1,2}  = B_2;  % B_2;
A2{1,2} = A2_2; % A2_2;

polyContainer = idpoly( 1, B, [], [], A2 );
dataContainer = iddata( y, [x1 x2] );

% Determine the parameters to estimate.
polyContainer.Structure.B(1,1).Free = B_1;
polyContainer.Structure.B(1,2).Free = B_2;

% Estimate the polynomials from data.
model = pem( dataContainer, polyContainer );
present(model);

% Residuals
remove = max([length(model.B(1)) length(model.B(2))]);
ehat = resid(model, dataContainer); ehat = ehat.OutputData; ehat = ehat(remove:end);

figure('Position', [100 100 800 400])
subplot(211)
acf( ehat, 30, 0.05, 1 );
title( sprintf('ACF'))
subplot(212)
pacf( ehat, 30, 0.05, 1 );
title( sprintf('PACF'))

%% Estimate C1 and A1
clc; close all;

A1_params = [1 2 5 6 7 23 24 25 26];
C1_params = [];

A1_0    = [1 zeros(1,30)]; A1_0(A1_params+1) = 1;
C1_0    = [1 zeros(1,30)]; C1_0(C1_params+1) = 1;

% Add differential operator?
% s = 24;
% s_poly = [1 zeros(1,s-1) -1];
% A1_0 = conv(A1_0,s_poly);

polyContainer = idpoly( 1, B, C1_0, A1_0, A2 );
dataContainer = iddata( y, [x1 x2] );

% Determine the parameters to estimate.
polyContainer.Structure.B(1,1).Free = B_1;
polyContainer.Structure.B(1,2).Free = B_2;
polyContainer.Structure.C.Free = C1_0;
polyContainer.Structure.D.Free = A1_0;
% Estimate the polynomials from data.
foundModel = pem( dataContainer, polyContainer );
present(foundModel);

% Residuals
remove = max([length(model.B(1)) length(model.B(2))]);
ehat = resid(foundModel, dataContainer); ehat = ehat.OutputData; ehat = ehat(remove:end);

figure('Position', [100 100 800 400])
subplot(211)
acf( ehat, 30, 0.05, 1 );
title( sprintf('ACF'))
subplot(212)
pacf( ehat, 30, 0.05, 1 );
title( sprintf('PACF'))
figure; whitenessTest(ehat); close;

%% k-step prediction
clc; close all;

% -------- ---- Choose k ------- ----%

k = 7;

% ---- Set input and output data ----%

data_set = 2;

% ---- ------------------------- ----%

% Pre-differentiate trend?
% y = filter([1 -1], 1, y);
% y = y(2:end);
% x = x(2:end);

% ---- ------------------------- ----%

switch data_set
	case 1
		x1 = atemp_m;
		x2 = wtemp_m;
		y = pow_m;
	case 2
		x1 = atemp_v;
		x2 = wtemp_v;
		y = pow_v;
	case 3
		x1 = atemp_t;
		x2 = wtemp_t;
		y = pow_t;
	case 4
		x1 = atemp_t2;
		x2 = wtemp_t2;
		y = pow_t2;
end

% Form polynomials
KA = conv( conv( foundModel.D, foundModel.F{1}), foundModel.F{2} );
KB = conv( conv( foundModel.D, foundModel.B{1}), foundModel.F{2} );
KC = conv( conv( foundModel.F{1}, foundModel.F{2}), foundModel.C );
KD = conv( conv( foundModel.D, foundModel.B{2}), foundModel.F{1} );

[Fy, Gy]   = polydivision( foundModel.C, foundModel.D, k );
[Fh1, Gh1] = polydivision( conv(Fy, KB), KC, k );
[Fh2, Gh2] = polydivision( conv(Fy, KD), KC, k );

% Predict the input signals.
[Fx1, Gx1] = polydivision( inputModel1.C, inputModel1.A, k );
xhatk1 = filter(Gx1, inputModel1.C, x1);
[Fx2, Gx2] = polydivision( inputModel2.C, inputModel2.A, k );
xhatk2 = filter(Gx2, inputModel2.C, x2);

% Form the predicted output signal using the predicted input signals.
yhatk  = filter(Fh1, 1, xhatk1) + filter(Gh1, KC, x1) + ...
         filter(Fh2, 1, xhatk2) + filter(Gh2, KC, x2) + ...
         filter(Gy, KC, y);

% Residuals
remove = max([length(Fh1), length(Fh2), length(Gh1), length(Gh2), length(Gy)]);
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
plot([x1(indV) xhatk1(indV)] )
legend('Input signal', 'Predicted input')
title( sprintf('Predicted input signal, x_{1,t+%i|t}', k) )
xlim(xindV);

figure('Position', [100 100 800 400])
plot([x2(indV) xhatk2(indV)] )
legend('Input signal', 'Predicted input')
title( sprintf('Predicted input signal, x_{2,t+%i|t}', k) )
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

%% Save model

dataLocation = pwd+"/data";
save(fullfile(dataLocation, 'BJdualmodel.mat'),'foundModel');
