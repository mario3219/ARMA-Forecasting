clc; clear; close all;

addpath(pwd+"/src", pwd+"/data")

% Load data
load('modellingData.mat');
load('validationData.mat');
load('testData.mat');
load('testData2.mat');

% Load trained model
foundModel = load("BJdualmodel.mat").foundModel;
inputModel1 = load("modellatemp.mat").modellatemp;
inputModel2 = load("modellwtemp.mat").modellwtemp;

present(foundModel);
present(inputModel1);
present(inputModel2);

%% k-step predictions
clc; close all;

% -------- ---- Choose k ------- ----%

% Works up until k=8
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
[~, Gx1] = polydivision( inputModel1.C, inputModel1.A, k );
[~, Gx2] = polydivision( inputModel2.C, inputModel2.A, k );

% Form input predictions
xhatk1 = filter(Gx1, inputModel1.C, x1);
xhatk2 = filter(Gx2, inputModel2.C, x2);

% Initial state values
init = [-0.956711880765908 0.212360309019174 0.137636395772135 ...
-0.003667042457479 -0.028144662871807 -0.117420928603422 ...
-0.164551365292900 0.087661594209967 0.078194102173014 ...
-0.213780201883578 -0.118368199424600 -0.002231786731556 ...
0.044725970850615 0.009103869946978 -0.037318505202701 ...
-0.033429437071952 -0.020378921667459 -0.046785348623816 ...
-0.113033665523572 0.129351158516263 0.086311195515968 ...
0.143285836615433 0.139816217926296 0.097807860320679 ...
-0.029659993673178 0.480911375240776 -0.133835324727756 ...
-0.124407274810065 0.027773398260632 0.094330999893900 ...
0.044592929106452 -0.096897592768211 -0.069692620335402 ...
0.011496398109199 0.080808219408271];

% Tracked states to remove. These indexes correspond
% to states that were estimated as zero, and are removed
% consecutively, ie index 34, array shifts, then index 14,
% array shifts, and so on.

remove = [34 14 12 27 4];

for i = 1:length(remove)
	init(remove(i)) = [];
end

N = length(y);
noPar = length(init);
startInd = 29;                                   % We use t-28, so start at t=29.
xt      = zeros(noPar, N-k);
xt(:,startInd-1) = init;

A     = eye(noPar);
Rw    = 800;                                % Measurement noise covariance matrix, R_w.
Re    = 1e-5*eye(noPar);                        % System noise covariance matrix, R_e.
Rx_t1 = 1e-4*eye(noPar);                        % Initial covariance matrix, R_{1|0}^{x,x}
Rx_k  = Rx_t1;
h_et  = zeros(N,1);                             % Estimated one-step prediction error.
yhat  = zeros(N-k,1);                           % One-step prediction \hat{y}_{t|t-1}
yhatk = zeros(N-k,1);                           % One-step prediction \hat{y}_{t+1|t}
xStd    = zeros(noPar,N-k);                       % Stores one std for the one-step prediction.
for t=startInd:N-k

    x_t1 = A*xt(:,t-1);                         % x_{t|t-1} = A x_{t-1|t-1}
    D   = [-y(t-1) -y(t-2) -y(t-5) -y(t-6) -y(t-7) -y(t-23) -y(t-24) -y(t-25) -y(t-26)];
    B1  = [x1(t) x1(t-1) x1(t-2) x1(t-3) x1(t-4) x1(t-5) x1(t-6) x1(t-7) x1(t-8) x1(t-9)...
	    x1(t-23) x1(t-24) x1(t-25) x1(t-26) x1(t-27) x1(t-28)];
    B2  = [x2(t) x2(t-1) x2(t-2) x2(t-5) x2(t-6) x2(t-7) x2(t-23) x2(t-24) x2(t-25) x2(t-26)];
    C = [D B1 B2];
        % Remove states
	for i = 1:length(remove)
		C(remove(i)) = [];
	end
    % Update the parameter estimates.
    Ry = C*Rx_t1*C' + Rw;                       % R_{t|t-1}^{y,y} = C R_{t|t-1}^{x,x} + Rw
    Kt = Rx_t1*C'/Ry;                           % K_t = R^{x,x}_{t|t-1} C^T inv( R_{t|t-1}^{y,y} )
    yhat(t) = C*x_t1;
    h_et(t) = y(t)-yhat(t);                    % One-step prediction error, \hat{e}_t = y_t - \hat{y}_{t|t-1}
    xt(:,t) = x_t1 + Kt*( h_et(t) );            % x_{t|t}= x_{t|t-1} + K_t ( y_t - Cx_{t|t-1} ) 

    % Update the covariance matrix estimates.
    Rx_t  = Rx_t1 - Kt*Ry*Kt';                  % R^{x,x}_{t|t} = R^{x,x}_{t|t-1} - K_t R_{t|t-1}^{y,y} K_t^T
    Rx_t1 = A*Rx_t*A' + Re;                     % R^{x,x}_{t+1|t} = A R^{x,x}_{t|t} A^T + Re

    % Form the k-step prediction by first constructing the future C vector
    % and the one-step prediction.
    D   = [-y(t) -y(t-1) -y(t-4) -y(t-5) -y(t-6) -y(t-22) -y(t-23) -y(t-24) -y(t-25)];
    B1  = [xhatk1(t+1) x1(t) x1(t-1) x1(t-2) x1(t-3) x1(t-4) x1(t-5) x1(t-6) x1(t-7) x1(t-8)...
	    x1(t-22) x1(t-23) x1(t-24) x1(t-25) x1(t-26) x1(t-27)];
    B2  = [xhatk2(t+1) x2(t) x2(t-1) x2(t-4) x2(t-5) x2(t-6) x2(t-22) x2(t-23) x2(t-24) x2(t-25)];
    Ck  = [D B1 B2];

    for i = 1:length(remove)
	Ck(remove(i)) = [];
    end

    yk = Ck*xt(:,t);                            % \hat{y}_{t+1|t} = C_{t+1|t} A x_{t|t}
    % the loop will require past predictions of y
    yhat_tmp = zeros(k,1);
    Rx_k = Rx_t1;
    for k0=2:k
        D   = [-yk -y(t-2+k0) -y(t-5+k0) -y(t-6+k0) -y(t-7+k0) -y(t-23+k0) -y(t-24+k0) -y(t-25+k0) -y(t-26+k0)];
        B1  = [xhatk1(t+k0) xhatk1(t-1+k0) x1(t-2+k0) x1(t-3+k0) x1(t-4+k0) x1(t-5+k0) x1(t-6+k0) x1(t-7+k0)...
		x1(t-8+k0) x1(t-9+k0) x1(t-23+k0) x1(t-24+k0) x1(t-25+k0) x1(t-26+k0) x1(t-27+k0) x1(t-28+k0)];
        B2  = [xhatk2(t+k0) xhatk2(t-1+k0) x2(t-2+k0) x2(t-5+k0) x2(t-6+k0) x2(t-7+k0)...
		x2(t-23+k0) x2(t-24+k0) x2(t-25+k0) x2(t-1+k0)];
	% Update C array to use future values, depending on current k
        if k0>2 D(2)=-yhat_tmp(k0-1); B1(3)=xhatk1(t-2+k0); B2(3)=xhatk2(t-2+k0); end
	    if k0>3 B1(4)=xhatk1(t-3+k0); end
	    if k0>4 B1(5)=xhatk1(t-4+k0); end
        if k0>5 D(3)=-yhat_tmp(k0-4); B1(6)=xhatk1(t-5+k0); B2(4)=xhatk2(t-5+k0); end
        if k0>6 D(4)=-yhat_tmp(k0-5); B1(7)=xhatk1(t-6+k0); B2(5)=xhatk2(t-6+k0); end
        if k0>7 D(5)=-yhat_tmp(k0-6); B1(8)=xhatk1(t-7+k0); B2(6)=xhatk2(t-7+k0); end
        Ck = [D B1 B2];                           % Store the old prediction (because Ck needs y(t+1) and y(t+2))

	for i = 1:length(remove)
		Ck(remove(i)) = [];
	end

	yk = Ck*A^k*xt(:,t);                    % \hat{y}_{t+k|t} = C_{t+k|t} A^k x_{t|t}
        yhat_tmp(k0) = yk;
        Rx_k = A*Rx_k*A' + Re;                  % R_{t+k+1|t}^{x,x} = A R_{t+k|t}^{x,x} A^T + Re  
    end
    yhatk(t+k) = yk;

    % Estimate a one std confidence interval of the estimated parameters.
    xStd(:,t) = sqrt( diag(Rx_t) );             % This is one std for each of the parameters for the one-step prediction.

end

% Limits
remove = 100; % Remove a large amount, just to be safe
indV = remove:length(y);

% Naive predictors
[yNaive, var_naive,~] = naivePred(y, indV, k);
[~, var_naive_air,~] = naivePred(x1, indV, k);
[~, var_naive_water,~] = naivePred(x2, indV, k);

% Residuals
ey = y(indV)-yhatk(indV);

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

figure('Position', [100 100 800 400])
plotWithConf( 1:N-k, xt', xStd', [] );
axis([startInd-1 N -0.3 0.75])
title(sprintf('Estimated parameters, with Re = %7.6f and Rw = %4.3f', Re(1,1), Rw(1,1)))
xlabel('Time')

% Code used to find the zero states when forming the initial state values

initial = xt(:,end);
initial = abs(initial);
idx = 1:length(initial);
toremove = [idx' initial];
toremove = sortrows(toremove, 2);
