function [R, R_ub, sigma_max, trace_Tn] = ESA(Nr, H, N, K, q, rho)

I = eye(K);
m_vec = Nr/N * ones(1,N);
[W, sigma_max, trace_Tn] = fatorize_AC(H, m_vec, q, Nr, N, K, rho);
He = W'*H; % effective channel
R = log2(det(I + rho * (He'*He)));
R_ub = N*log2(1 + 1/N*rho*sum(eig(He'*He)));


end %eof

