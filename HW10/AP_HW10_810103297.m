clc, clear, close all;

%% Signal Properties

c = 3e8;
fc = 150e6;
fs = 1e6;
ts = 1/fs;

tau = 1e-6;
PRI = 250e-6;
Trec = 25e-3;
delta_f=1/Trec;

tau_num = round(tau/ts);
num_pulse = round(Trec/PRI);

t = 0:ts:Trec-ts;
f = -fs/2:delta_f:fs/2-delta_f;
sl_PRI = [ones(1, round(tau/ts)), zeros(1, round(PRI/ts)-round(tau/ts))];
sl = repmat(sl_PRI, 1, num_pulse);
slf = fftshift(fft(sl));

%% Array Properties

lambda = c/fc;
k = 2*pi/lambda;

Mx = 20;
Mz = 10;
M = Mx*Mz;

delta_d = lambda/2;
dx = 0:delta_d:(Mx-1)*delta_d;
dz = 0:delta_d:(Mz-1)*delta_d;

delta_theta = 1;
theta = -90:delta_theta:90;

delta_psi = 1;
psi = -90:delta_psi:90;

A = zeros(M, length(psi)*length(theta));

for i=1:length(theta)
    for j=1:length(psi)
        Temp = exp(-1j*k*dx.'*sin(pi*psi(j)/180))...
            *exp(-1j*k*dz*sin(pi*theta(i)/180));
        A(:, (i-1)*length(psi)+j) = Temp(:);
    end
end

%% Question 1-1

% 1st target
r_1 = 6e3;
td_1 = 2*r_1/c;
fd_1 = 400;
alpha_1 = 1;

theta_1 = 45;
phi_1 = 37;
psi_1 = asind(cosd(theta_1)*cosd(phi_1));

yr_1 = alpha_1*exp(1j*2*pi*fd_1*t).*circshift(sl, round(td_1/ts));

% 2nd target
r_2 = 18e3;
td_2 = 2*r_2/c;
fd_2 = 1200;
alpha_2 = 1;

theta_2 = 30;
phi_2 = 53;
psi_2 = asind(cosd(theta_2)*cosd(phi_2));

yr_2 = alpha_2*exp(1j*2*pi*fd_2*t).*circshift(sl, round(td_2/ts));

Temp = exp(-1j*k*dx.'*sin(pi*psi_1/180))*exp(-1j*k*dz*sin(pi*theta_1/180));
ar_1 = Temp(:);
Yr_1 = ar_1*yr_1;

Temp = exp(-1j*k*dx.'*sin(pi*psi_2/180))*exp(-1j*k*dz*sin(pi*theta_2/180));
ar_2 = Temp(:);
Yr_2 = ar_2*yr_2;

Yr = Yr_1 + Yr_2;

%% Range-Doppler Processing

COEFF_MEMORY = [];
for m = 1:M
    yr = Yr(m, :);
    [COEFF, DOPPLER, RANGE] = Range_Doppler_Process(yr, sl, ts, PRI, c);
    COEFF_MEMORY(:, :, m) = COEFF;
end

COEFF_ALL = sum(abs(COEFF_MEMORY).^2, 3);
figure
subplot(1, 2, 1)
mesh(DOPPLER, RANGE, COEFF_ALL)
ylabel('Range')
xlabel('Doppler')
subplot(1, 2, 2)
imagesc(DOPPLER, RANGE, COEFF_ALL)
ylabel('Range')
xlabel('Doppler')
sgtitle("Range-Doppler Matrix Magnitude")

[idx_R, idx_D] = find(COEFF_ALL > 1e-3);
K = length(idx_R);

%% Matching Pursuit

theta_est = zeros(K,1);
psi_est   = zeros(K,1);
phi_est   = zeros(K,1);

for kk=1:K
    y = squeeze(COEFF_MEMORY(idx_R(kk), idx_D(kk), :));
    idx_A = MP(y(:), A, 1);

    i_theta = ceil(idx_A / length(psi));
    j_psi   = idx_A - (i_theta-1)*length(psi);

    theta_est(kk) = theta(i_theta);
    psi_est(kk) = psi(j_psi);
    phi_est(kk) = acosd(sind(psi_est(kk))/cosd(theta_est(kk)));
end

disp("Matching Pursuit:")
for kk=1:K
    disp([' -Target ', num2str(kk)])
    disp(['     Range   = ', num2str(RANGE(idx_R(kk))), ' m'])
    disp(['     Doppler = ', num2str(DOPPLER(idx_D(kk))), ' Hz'])
    disp(['     theta   = ', num2str(theta_est(kk)), ' deg'])
    disp(['     psi     = ', num2str(psi_est(kk)), ' deg'])
    disp(['     phi     = ', num2str(phi_est(kk)), ' deg'])
end

%% SVD

theta_est = zeros(K,1);
psi_est   = zeros(K,1);
phi_est   = zeros(K,1);

d = 1;

for kk = 1:K
    y = squeeze(COEFF_MEMORY(idx_R(kk), idx_D(kk), :));
    Y = reshape(y, Mx, Mz);

    [U,S,V] = svd(Y, 'econ');

    Unull = U(:, d+1:end);
    Vnull = V(:, d+1:end);

    ax = exp(-1j*k*dx.'*sind(psi));
    psi_music = 1 ./ (abs(diag(ax'*(Unull*Unull')*ax)) + 1e-12);

    az = exp(+1j*k*dz.'*sind(theta));
    theta_music = 1 ./ (abs(diag(az'*(Vnull*Vnull')*az)) + 1e-12);

    [~, idx_P] = max(psi_music);
    [~, idx_T] = max(theta_music);

    psi_est(kk)   = psi(idx_P);
    theta_est(kk) = theta(idx_T);
    phi_est(kk) = acosd(sind(psi_est(kk))/cosd(theta_est(kk)));
end

disp("SVD+MUSIC:")
for kk=1:K
    disp([' -Target ', num2str(kk)])
    disp(['     Range   = ', num2str(RANGE(idx_R(kk))), ' m'])
    disp(['     Doppler = ', num2str(DOPPLER(idx_D(kk))), ' Hz'])
    disp(['     theta   = ', num2str(theta_est(kk)), ' deg'])
    disp(['     psi     = ', num2str(psi_est(kk)), ' deg'])
    disp(['     phi     = ', num2str(phi_est(kk)), ' deg'])
end

figure
subplot(1, 2, 1)
plot(psi, psi_music)
xlabel("psi (deg)")
grid on
subplot(1,2,2)
plot(theta, theta_music)
xlabel("theta (deg)")
grid on
sgtitle("Music Algorithm Performance on Target 2")

%% Uniform Array

theta_est = zeros(K,1);
psi_est   = zeros(K,1);
phi_est   = zeros(K,1);

for kk = 1:K
    y = squeeze(COEFF_MEMORY(idx_R(kk), idx_D(kk), :));
    y = y(:);

    Y = reshape(y, Mx, Mz);

    S = fftshift(fft2(Y));
    P = abs(S).^2;
    
    [~, idxMax] = max(P(:));
    [ix, iz] = ind2sub(size(P), idxMax);
    

    
    fx = (ix - (Mx/2+1))/Mx;   % for even Mx
    fz = (iz - (Mz/2+1))/Mz;   % for even Mz
    
    psi_est(kk) = asind(-(lambda/delta_d) * fx);
    theta_est(kk) = asind(-(lambda/delta_d) * fz);
    phi_est(kk) = acosd(sind(psi_est(kk))/cosd(theta_est(kk)));
end

figure
subplot(1,2,1)
mesh( (-Mz/2:Mz/2-1)/Mz, (-Mx/2:Mx/2-1)/Mx, 10*log10(P./max(P(:))) );
xlabel('f_z'); ylabel('f_x'); zlabel('dB');
grid on
subplot(1,2,2)
imagesc( (-Mz/2:Mz/2-1)/Mz, (-Mx/2:Mx/2-1)/Mx, 10*log10(P./max(P(:))) );
axis xy; colorbar;
xlabel('f_z'); ylabel('f_x');
sgtitle("2D Spatial Spectrum for Target 2 (dB, normalized)")

disp("Uniform Array:")
for kk=1:K
    disp([' -Target ', num2str(kk)])
    disp(['     Range   = ', num2str(RANGE(idx_R(kk))), ' m'])
    disp(['     Doppler = ', num2str(DOPPLER(idx_D(kk))), ' Hz'])
    disp(['     theta   = ', num2str(theta_est(kk)), ' deg'])
    disp(['     psi     = ', num2str(psi_est(kk)), ' deg'])
    disp(['     phi     = ', num2str(phi_est(kk)), ' deg'])
end