function [W, sigma_max, trace_Tn, trace_Qn] = fatorize_AC(H, m_vec, q, Nr, N, K, rho)

W = zeros(Nr,N);
I = eye(K);
Q = I;
r_end = 0;
sigma_max = zeros(N,1);
trace_Tn = zeros(N,1);
trace_Qn = zeros(N,1);
% Tn_mean = zeros(Nr/N,Nr/N);
Qn = zeros(K,K,N);
Gn = zeros(K,K,N);
for nn = 1:N
    M = m_vec(nn);
%     m_vec
%     sum(m_vec)
    r_start = r_end + 1;
    r_end = r_start + M - 1;
    Dict = 1/sqrt(M) * exp(1i * 2*pi/q * [0:q-1]).';
    
    H_n = H(r_start:r_end,:);
    T_n = H_n*(Q^(-1))*H_n';
    [Ut,Dt,~] = svd(T_n);
    
    sigma_vec_tmp = diag(Dt);
    sigma_vec = sigma_vec_tmp(sigma_vec_tmp>1e-4);
    
    sigma_max(nn) = max(sigma_vec);
    trace_Tn(nn) = real(trace(T_n));
    trace_Qn(nn) = trace(Q^-1);
    
    
    u_n = Ut(:,1);
    w_n = quant_sub(M,Dict,N,u_n); % quantize u
    W(r_start:r_end, nn) = w_n;
    
    g_n = H_n'*w_n;
    G_n = g_n*g_n';
    Gn(:,:,nn) = G_n;
    E_n = I + rho*Q^(-1)*G_n;
    %Q^(-1)
    Q = Q*E_n;
    if nn == 1
        Qn(:,:,nn) = E_n;
    else
        Qn(:,:,nn) = Qn(:,:,nn-1)*E_n;
    end
end
end % eof