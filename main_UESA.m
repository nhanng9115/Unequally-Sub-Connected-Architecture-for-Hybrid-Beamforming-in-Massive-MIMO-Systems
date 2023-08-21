clc;
clear all;
warning off
% close all;
% figure
Nt = 1;
Nr = 64;
N = 4; K = 4;
do_find_M = 0;
%% load m space
% file_name = strcat('.\data\dat_m_',num2str(Nr),'_', num2str(N));
% if do_find_M == 1
%     [m_space_ES, m_space_RES, m_space_RES_EE] = find_m_space(Nr,N);
%     save(file_name,'m_space_ES','m_space_RES','m_space_RES_EE');
% end
% m_dat = load(file_name);
% m_dat = [];
[count_max, n_stop] = get_params(Nr);
%%

n_path = 10; % number of paths
q = 2^4;
ITER = 400;
SNR_dB_vec = 0:2:12;
fig = [1, 0];
gamma = 3;
n_update_max = 80;

%% schemes:      ESA	UESA-ES     RES     RES-ET	Fast-UESA   AS      SVD     AHP
combine_scheme = [1,	1,          0,      0,      0,          0,      0,      0];
n_scheme = length(combine_scheme);

R = zeros(n_scheme,ITER,length(SNR_dB_vec));
R_ub = zeros(n_scheme,ITER,length(SNR_dB_vec));
sigma = zeros(n_scheme,ITER,N);
trace = zeros(n_scheme,ITER,N);
n_cand = zeros(ITER,1);

for ss = 1:length(SNR_dB_vec)
    snr_dB = SNR_dB_vec(ss);
    disp(snr_dB)
    snr = 10^(snr_dB/10);
    Es = 1; % normalized transmit power
    N0 = Es / snr; % noise variance
    rho = 1/N0;
    
    M_vec = zeros(1,N);
    
    for ii = 1:ITER
        
        [H, a_rx_best] = genChannel(Nr,Nt,K,n_path,1);
        %H = 1/sqrt(2)*(randn(Nr,K) + 1i*randn(Nr,K));
        %% ESA
        if combine_scheme(1) == 1
            [R(1,ii,ss), R_ub(1,ii,ss), sigma(1,ii,:), trace(1,ii,:)] = ESA(Nr, H, N, K, q, rho);
        end
        
        %% UESA
        %H = order_channel(H,'row','descend'); % row, element
        
        %% UESA-ES
        if combine_scheme(2) == 1
            [R(2,ii,ss), R_ub(2,ii,ss), sigma(2,ii,:), trace(2,ii,:)] = ...
                ESA_test(Nr, H, N, K, q, rho);
            %UESA(Nr, H, N, K, q, rho, 'ES', 0, m_dat);
        end
        
        %% UESA-RES
        if combine_scheme(3) == 1
            [R(3,ii,ss), R_ub(3,ii,ss), sigma(3,ii,:), trace(3,ii,:), n_cand(ii), m_vec] = ...
                UESA(Nr, H, N, K, q, rho, 'RES', 0, m_dat);
        end
        
        %% UESA-RES-EE
        if combine_scheme(4) == 1
            [R(4,ii,ss), R_ub(4,ii,ss), sigma(4,ii,:), trace(4,ii,:), n_cand(ii), m_vec] = ...
                UESA(Nr, H, N, K, q, rho, 'ET', count_max, m_dat);
        end
        
        %% Fast-UESA
        if combine_scheme(5) == 1
            [R(5,ii,ss), R_ub(5,ii,ss), sigma(5,ii,:), trace(5,ii,:), trace_Q, n_cand(ii), m_vec] = ...
                Fast_UESA(Nr, H, N, K, q, rho, m_dat, gamma, n_update_max);
        end
        
        %% AS
        if combine_scheme(6) == 1
            R(6,ii,ss) = UESA_AS(H, n_bit, Nr, N, K, rho);
        end
        
        %% SVD
        if combine_scheme(7) == 1
            R(7,ii,ss) = SVD(H, q, Nr, N, K, rho, m_vec);
        end
        
        if combine_scheme(8) == 1
            R(8,ii,ss) = AHP(Nr, H, N, K, q, rho, count_max, 1, m_dat);
        end
        %m_vec(1)
    end
end

n_cand_mean = mean(n_cand,1);


sigma_mean = zeros(N,n_scheme);
trace_mean = zeros(N,n_scheme);
for n = 1:n_scheme
    for i = 1:N
        sigma_mean(i,n) = mean(sigma(n,:,i));
        trace_mean(i,n) = mean(trace(n,:,i));
    end
end

% mean rate
R_mean = zeros(length(SNR_dB_vec), n_scheme);
R_ub_mean = zeros(length(SNR_dB_vec), n_scheme);
for n = 1:n_scheme
    for s = 1:length(SNR_dB_vec)
        R_mean(s,n) = mean(R(n,:,s));
        R_ub_mean(s,n) = mean(R_ub(n,:,s));
    end
end
% figure
% plot_fig(fig, SNR_dB_vec, N, R_mean, R_ub_mean, sigma_mean, trace_mean, combine_scheme);

%% plot figures to compare Fast-UESA, AS, and SVD
hold all
line_style = {'-ko','--ko', '-gs','--gs','-g*','--g*','-b^','--b^','-r+','--m+','-gp','-gp'};
leg_str = {'','','','','Fast-UESA', '', 'SVD', 'AHP'};
line = zeros(n_scheme,1);
leg = zeros(n_scheme,1);
figure
for n = 1:n_scheme
    if combine_scheme(n) == 1
        line(n) = semilogy(SNR_dB_vec, R_mean(:,n), line_style{n}, 'LineWidth', 1.4); hold on;
        leg(n) = 1;
    end
end

legend(line(line~=0), leg_str(leg~=0));
hold all;

xlabel('SNR [dB]');
ylabel('Total achievable rate [bps/Hz');
box on
grid on
axis tight


