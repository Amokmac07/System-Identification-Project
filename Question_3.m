%% Question 3
clc;clear;
rng(1)
load dryer_data.mat;
%% Load data
uk_1 = IO_data(:,1);
uk_2 = IO_data(:,2);
uk_3 = IO_data(:,3);
yk_1 = IO_data(:,4);
yk_2 = IO_data(:,5);
yk_3 = IO_data(:,6);

%% Estimating the delays
% Define a function to calculate delay using cross-correlation
% Sampling time
Ts = 0.1;

% Calculate delay for each input-output pair
delay_uk1_yk1 = calculate_delay(uk_1, yk_1, Ts);
delay_uk2_yk1 = calculate_delay(uk_2, yk_1, Ts);
delay_uk3_yk1 = calculate_delay(uk_3, yk_1, Ts);

delay_uk1_yk2 = calculate_delay(uk_1, yk_2, Ts);
delay_uk2_yk2 = calculate_delay(uk_2, yk_2, Ts);
delay_uk3_yk2 = calculate_delay(uk_3, yk_2, Ts);

delay_uk1_yk3 = calculate_delay(uk_1, yk_3, Ts);
delay_uk2_yk3 = calculate_delay(uk_2, yk_3, Ts);
delay_uk3_yk3 = calculate_delay(uk_3, yk_3, Ts);
%% Define parameters for each output

% We have to define the Laguerre function properties for each output
% differently as the delay and poles are  different 
%
p = 0.6; % Pole
n1 = 3; % 1st input
n2 = 2; % 2nd input
n3 = 2; % 3rd input
d = 10*0.1; % Process delay
Ts = 0.1; % Sampling time

% Generate Laguerre basis functions
F1 = laguerre(p, n1, d);
F2 = laguerre(p, n2, d);
F3 = laguerre(p, n3, d);

Phi1 = construct_regression_matrix(F1, uk_1);
Phi2 = construct_regression_matrix(F2, uk_2);
Phi3 = construct_regression_matrix(F3, uk_3);
Phi_t = [Phi1 Phi2 Phi3];

% Process each output individually

%% For yk_1

% Estimating coefficients for yk_1
c1 = Phi_t \ yk_1;
c1_1 = c1(1:n1);
c1_2 = c1(n1+1:n1+n2);
c1_3 = c1(n1+n2+1:end);

% Construct the MISO OBF model for yk_1
obf_model1 = construct_obf_model(c1_1, F1, c1_2, F2, c1_3, F3);
obf_model1.Ts = Ts;


% Residual analysis of MISO OBF model for yk_1
dataset1 = iddata(yk_1, [uk_1 uk_2 uk_3], Ts);

% Comparing model 
compare(dataset1,obf_model1);
vk1 = resid(dataset1, obf_model1);

% Plot autocorrelation of residuals for yk_1
figure;
autocorr(vk1.OutputData);

% Fit an AR model for the residuals of yk_1
mdl_ar1 = ar(vk1.OutputData, 3);
figure;
resid(vk1.OutputData, mdl_ar1);

%% For yk_2
c2 = Phi_t \ yk_2;
c2_1 = c2(1:n1);
c2_2 = c2(n1+1:n1+n2);
c2_3 = c2(n1+n2+1:end);

% Construct the MISO OBF model for yk_2
obf_model2 = construct_obf_model(c2_1, F1, c2_2, F2, c2_3, F3);
obf_model2.Ts = Ts;

% Residual analysis of MISO OBF model for yk_2
dataset2 = iddata(yk_2, [uk_1 uk_2 uk_3], Ts);

% Comparing model 
compare(dataset2,obf_model2);

vk2 = resid(dataset2, obf_model2);

% Plot autocorrelation of residuals for yk_2
figure;
autocorr(vk2.OutputData);

% Fit an AR model for the residuals of yk_2
mdl_ar2 = ar(vk2.OutputData, 2);
figure;
resid(vk2.OutputData, mdl_ar2);

%% For yk_3
c3 = Phi_t \ yk_3;
c3_1 = c3(1:n1);
c3_2 = c3(n1+1:n1+n2);
c3_3 = c3(n1+n2+1:end);

% Construct the MISO OBF model for yk_3
obf_model3 = construct_obf_model(c3_1, F1, c3_2, F2, c3_3, F3);
obf_model3.Ts = Ts;

% Residual analysis of MISO OBF model for yk_3
dataset3 = iddata(yk_3, [uk_1 uk_2 uk_3], Ts);

% Comparing model 
compare(dataset3,obf_model3);

vk3 = resid(dataset3, obf_model3);

% Plot autocorrelation of residuals for yk_3
figure;
autocorr(vk3.OutputData);

% Fit an AR model for the residuals of yk_3
mdl_ar3 = ar(vk3.OutputData, 3);
figure;
resid(vk3.OutputData, mdl_ar3);


%% Functions
function delay = calculate_delay(input, output, Ts)
    [cross_corr, lags] = xcorr(output, input, 'normalized');
    [~, max_idx] = max(cross_corr);
    delay = lags(max_idx) * Ts;
end
% Function to construct regression matrix using Laguerre basis functions
function Phi = construct_regression_matrix(F, uk)
    n = length(F);
    Phi = zeros(length(uk), n);
    for i = 1:n
        y_est = sim(idpoly(F{i}), uk); % Simulate the ith Laguerre basis
        Phi(:, i) = y_est; % Populate the regression matrix
    end
end

% Function to construct the OBF model using estimated coefficients and Laguerre basis functions
function obf_model = construct_obf_model(c1, F1, c2, F2, c3, F3)
    obf_model = 0;
    for a = 1:length(c1)
        obf_model = obf_model + c1(a) * F1{a};
    end
    for b = 1:length(c2)
        obf_model = obf_model + c2(b) * F2{b};
    end
    for c = 1:length(c3)
        obf_model = obf_model + c3(c) * F3{c};
    end
end



% To improve the predictions :-
% 1. the delay and pole relations between the outputs and inputs must be carefully implemented to get a bespoke laguerre function for each relation 
% 2. we also need to account for the hidden relations hence the truncation
% error resulting from assumption of a MIMO system as 3 diiferent MISO
% systems 
