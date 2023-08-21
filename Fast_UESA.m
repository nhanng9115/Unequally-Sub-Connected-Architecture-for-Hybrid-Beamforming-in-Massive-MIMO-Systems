function [R1, R_ub, sigma, trace_T, trace_Q, n_cand, m_vec] = Fast_UESA(Nr, H, N, K, q, rho, m_dat, gamma, n_update_max)

m_space = m_dat.m_space_RES_EE;
m_init = m_space(1,:);
n_cand = 0;

% m1 = 1; m2 = 1; m3 = 1; m4 = 1; m5 = 1; m6 = 1; m7 = 1; m8 = Nr - (m1+m2+m3+m4+m5+m6+m7);
% m_init = [m1,m2,m3,m4,m5,m6,m7,m8];

[W_tmp, sigma, trace_T] = fatorize_AC(H, m_init, q, Nr, N, K, rho);
He = W_tmp'*H; % effective channel
I = eye(N);
tau = sum(sigma);
R1 = log2(det(I + rho * (He'*He)));
R_ub = N*log2(1 + 1/N*rho*sum(eig(He'*He)));
trace_Q = zeros(N,1);
m_vec = m_init;
m_tmp = m_init;
for i = 1:n_update_max
    delta_sigma = sigma - mean(sigma);
    for n = 1:N
        if delta_sigma(n) < gamma
            m_tmp(n) = m_tmp(n) + 1;
            n_cand = n_cand + 1;
        end
        m_tmp(end) = Nr - sum(m_tmp(1:N-1));
        if m_tmp(end) < m_tmp(N-1)
            break
        end
        
        %m_tmp
        if isequal(ismember(m_tmp, m_space), [1,1,1,1])
            if m1<=m2 && m2<=m3 && m3<=m4 && m4<=m5 && m5<=m6 && m6<=m7 && m7<=m8
                [W_tmp, sigma_tmp, trace_Ttmp, trace_Qtmp] = fatorize_AC(H, m_tmp, q, Nr, N, K, rho);
                He = W_tmp'*H; % effective channel
                Rtmp = log2(det(I + rho * (He'*He)));
                sum_sm = sum(sigma_tmp);
                if Rtmp > R1
                    m_vec = m_tmp;
                    tau = sum_sm;
                    R1 = Rtmp;
                    R_ub = N*log2(1 + 1/N*rho*sum(eig(He'*He)));
                    sigma = sigma_tmp;
                    trace_T = trace_Ttmp;
                    trace_Q = trace_Qtmp;
                end
            end
        end
    end
    
end % eof