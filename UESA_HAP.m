function R = UESA_HAP(H, q, Nr, N, K, rho, m_dat, count_max)
% m_vec
W = zeros(Nr,N);
M = Nr/N;
I = eye(K);
Dict = 1/sqrt(M) * exp(1i * 2*pi/q * [0:q-1]).';

m_space = m_dat.m_space_RES;
n_cand = size(m_space,1);

R = 0;
m_vec = zeros(1,N);
count = 0;
for n = 1:n_cand
    m_tmp = m_space(n, :);
    Htmp = H;
    %r_end = 0;
    for k = 1:K
        M = m_tmp(k);
        %r_start = randi([1,Nr/4-1]);%r_end + 1;
        %r_end = randi([Nr/2,Nr]);%r_start + M - 1;

        for m = 1:M
            hk = Htmp(:,k); % channel vector of user k
            % find the channel coefficient with max abs
            [~, j0] = max(abs(hk));
            hk_max = hk(j0);
            
            %phase of hk_max
            phase_diff = abs(hk_max/abs(hk_max) - Dict);
            [~, n_hat] = min(phase_diff);
            W(j0,k) = Dict(n_hat);
            Htmp(j0,:) = 0;
        end
    end
    %W
    He = W'*H; % effective channel
    R_tmp = log2(det(I + rho * (He'*He)));
    if R_tmp > R
        R = R_tmp;
        m_vec = m_tmp;
        count = 0;
    else
        count = count + 1;
    end
    if count == count_max
        break
    end
end

end % eof