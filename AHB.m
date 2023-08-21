function R = AHB(Nr, H, N, K, q, rho, count_max, esa, m_dat, ET)

R = 0;
W = zeros(Nr,N);
M = Nr/N;
if rem(Nr,N) ~= 0
    error('error');
end
    
I = eye(K);
Htmp = H;

if esa == 1
    % consider K users
    phi_dict = 2*pi/q * [0:q-1];
    for m = 1:M
        % channel vector corresponding to user k
        for k = 1:K
            hk = Htmp(:,k);
            % find the channel coefficient with max abs
            [~, j0] = max(abs(hk));
            hk_max = hk(j0);
            
            %% paper's method ===============================
            a = hk_max/abs(hk_max);
            phase_diff = abs(a - exp(1i * 2*pi/q * [0:q-1]));
            [~, n_hat] = min(phase_diff);
            W(j0,k) = 1/sqrt(M) * exp(1i * 2*pi/q * n_hat);
            
            %% compare phases
            phi_hmax = angle(hk_max);            
            delta_phi = abs(wrapToPi(phi_dict - phi_hmax));
            [~, iMin] = min(delta_phi);
            W(j0,k) = 1/sqrt(M) * exp(1i * 2*pi/q * iMin);
            %% ==============================================
            Htmp(j0,:) = 0;
        end
    end
    % W
    He = W'*H; % effective channel
    R = log2(det(I + rho * (He'*He)));
else
    %count_max = Nr*2;
    m_space = m_dat.m_space_RES;
    n_cand = size(m_space,1);
    count = 0;
    for n = 1:n_cand
        m_tmp = m_space(n, :);
        Htmp = H;
        for k = 1:K
            M = m_tmp(k);
            dict = 1/sqrt(M) * exp(1i * 2*pi/q * [0:q-1]).';
            for m = 1:M
                hk = Htmp(:,k); % channel vector of user k
                [~, j0] = max(abs(hk));
                hk_max = hk(j0);
                
                %phase of hk_max
                phase_diff = abs(hk_max/abs(hk_max) - dict);
                [~, n_hat] = min(phase_diff);
                W(j0,k) = dict(n_hat);
                Htmp(j0,:) = 0;
            end
        end
        %W
        He = W'*H; % effective channel
        R_tmp = log2(det(I + rho * (He'*He)));
        if R_tmp > R
            R = R_tmp;
            count = 0;
        else
            count = count + 1;
        end
        if count == count_max && ET == 1
            break
        end
    end
end

end % eof