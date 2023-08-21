function [H] = order_channel(H,orderType,order)
Nr = size(H,1);
K = size(H,2);

switch orderType
    case 'row'
        normH = zeros(Nr,1);
        for rr = 1:Nr
            normH(rr) = norm(H(rr,:))^2;
        end
        [~,iH] = sort(normH,order);
        H = H(iH,:);
    case 'col'
        normH = zeros(K,1);
        for cc = 1:K
            normH(cc) = norm(H(cc,:))^2;
        end
        [~,iH] = sort(normH,order);
        H = H(:,iH);
    case 'element'
        for k = 1:K
            h_k = H(:,k);
            abs_h_k = abs(h_k).^2;
            [~,sort_ind] = sort(abs_h_k,order);
            h_k_sort = h_k(sort_ind);
            H(:,k) = h_k_sort;
        end
end


end