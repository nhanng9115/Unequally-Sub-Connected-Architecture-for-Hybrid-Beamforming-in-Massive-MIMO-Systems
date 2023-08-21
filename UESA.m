function [R, R_ub, sigma_max, trace_Tn, nExamCand, m_vec] = UESA(Nr, H, N, K, q, rho, search_scheme, count_max, m_dat)

T = 0;
count = 0;
I = eye(N);
break_flag = 0;
nExamCand = 0;
m_vec = zeros(1,N);

% load all candidates
if strcmp(search_scheme,'ET')
    m_space = m_dat.m_space_RES_EE;
elseif strcmp(search_scheme,'RES')
    m_space = m_dat.m_space_RES;
else
    m_space = m_dat.m_space_ES;
end

nCand = size(m_space,1);
for i = 1:nCand
    
    m_vec_tmp = m_space(i,:);
    
    if strcmp(search_scheme,'ET') %% UESA-RES-ET
        nExamCand = nExamCand + 1;
        if count < count_max
            [W_tmp, sigma_max_tmp, trace_Tn_tmp] = fatorize_AC(H, m_vec_tmp, q, Nr, N, K, rho);
            He = W_tmp'*H; % effective channel
            R_tmp = log2(det(I + rho * (He'*He)));
            R_ub_tmp = N*log2(1 + 1/N*rho*sum(eig(He'*He)));
            if R_tmp > T
                R = R_tmp;
                R_ub = R_ub_tmp;
                T = R_tmp;
                count = 0;
                sigma_max = sigma_max_tmp;
                trace_Tn = trace_Tn_tmp;
                m_vec = m_vec_tmp;
            else
                count = count + 1;
            end
        else
            %disp('exceed');
            break_flag = 1;
        end
        if break_flag == 1
            break;
        end
        
    else %% UESA-RES
        [W_tmp, sigma_max_tmp, trace_Tn_tmp] = fatorize_AC(H, m_vec_tmp, q, Nr, N, K, rho);
        He = W_tmp'*H; % effective channel
        R_tmp = log2(det(I + rho * (He'*He)));
        R_ub_tmp = N*log2(1 + 1/N*rho*sum(eig(He'*He)));
        if R_tmp > T
            T = R_tmp;
            R = R_tmp;
            R_ub = R_ub_tmp;
            sigma_max = sigma_max_tmp;
            trace_Tn = trace_Tn_tmp;
            m_vec = m_vec_tmp;
            %test = test + sum(log2(m_vec./(Nr/N)))/nCand;
        end
    end
end

end %eof

