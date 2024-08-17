%% Question 2
clc;clear;
rng(1)
%% Question 2 A
mod_dgp=tf([0.2122 0.7778 0.1781],[1 -2.6714 2.3772 -0.7047 0 0 0],1);
mod_dgp=idpoly(mod_dgp);
uk = idinput(1000,'prbs',[0 0.1],[-1 1]);
xk = sim(mod_dgp,uk);
mod_dgp.Noisevariance = var(xk)/10; %noise variance SNR 10
yk = sim(mod_dgp, uk, simOptions('AddNoise', true));
zk = iddata(yk,uk,1);
zk=detrend(zk,0);
figure
plot(zk)


% Estimate the impulse response using impulseest
mod_fir= impulseest(zk);

% Plot the estimated impulse response
figure;
impulse(mod_fir,80);
title('Estimated Impulse Response');

% from plot delay is 0
figure;
compare(zk,mod_fir)
%% Question B
% Generate Laguerre basis functions
p = 0.8;  % Pole
n = 5;    % Number of terms in Laguerre expansion
d = 0;    % Process delay

F = laguerre(p, n, d);
% Initialize the regression matrix
N = length(uk);
Phi = zeros(N, n);

% Construct the regression matrix
for i = 1:n
    y_est = sim(idpoly(F{i}), uk);  % Simulate the ith Laguerre basis filter
    Phi(:, i) = y_est;  % Populate the regression matrix
end
% Estimate the coefficients using least squares
c = Phi \ yk;

% Display the estimated coefficients
disp('Estimated coefficients:');
disp(c);

% Construct the OBF model
obf_model = 0;
for i = 1:n
    obf_model = obf_model + c(i) * F{i};
end

% Display the OBF model
% disp('OBF Model:');
% disp(obf_model);
 figure;
compare(zk,obf_model)

% Simulate the OBF model with the input data
yk_obf = sim(idpoly(obf_model), uk);

% Calculate the mean square error (MSE)
mse_obf = mean((yk - yk_obf).^2)/1000;

% Display the mean square error
fprintf('Mean Square Error of OBF Model: %.6f\n', mse_obf);

