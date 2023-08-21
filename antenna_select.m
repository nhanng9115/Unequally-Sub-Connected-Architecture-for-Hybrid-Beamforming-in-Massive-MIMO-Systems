function best_rows = antenna_select(H, N, rho)

%% Antenna selection
M = size(H,1);
I = 1:M;
B = eye(N);
A = H*H';
alpha = diag(A);
Itmp = I;
for n = 1:N
    [alpha_max, J] = max(alpha);
    Itmp = setdiff(Itmp,J);
    if n < N
        a = 1/sqrt(N*rho + alpha_max) * B * H(J,:)';
        B = B - a*a';
        
        for t = 1:length(Itmp)
            j = Itmp(t);
            alpha(j) = alpha(j) - abs(a'*H(j,:)')^2;
        end
    end
    alpha(J) = 0;
end
best_rows = setdiff(I, Itmp);

end % eof