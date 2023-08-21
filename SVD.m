function R = SVD(H, q, Nr, N, K, rho, m_vec)


R = 0;
W = zeros(Nr,N);
M = Nr/N;
I = eye(K);
Dict = 1/sqrt(1) * exp(1i * 2*pi/q * [0:q-1]).';
ESA = 0;

if ESA == 1
    % consider K users
    r_end = 0;
    for n = 1:N
        r_start = r_end + 1;
        r_end = r_start + M - 1;
        Hn = H(r_start:r_end,:);
        
        [U,D,~] = svd(Hn);
        [dmax, imax] = max(diag(D));
        w = U(:,imax);
        W(r_start:r_end,n) = w;
    end
    % W
    He = W'*H; % effective channel
    R = log2(det(I + rho * (He'*He)));
else
    
    r_end = 0;
    for n = 1:N
        M = m_vec(n);
        r_start = r_end + 1;
        r_end = r_start + M - 1;
        Hn = H(r_start:r_end,:);

        [U,D,~] = svd(Hn);
        w = U(:,1);

        W(r_start:r_end,n) = w;
    end
    % W
    He = W'*H; % effective channel
    R = log2(det(I + rho * (He'*He)));
end % eof