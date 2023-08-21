function [m_space_ES, m_space_RES, m_space_RES_EE] = find_m_space(Nr,N)
m_space_ES = [];
m_space_RES = [];
m_space_RES_EE = [];

delta_vec = [];
switch N
    case 2
        for m1 = 1:Nr-1
            m2 = Nr-m1;
            m_tmp = [m1,m2];
            m_space_ES = cat(1,m_space_ES,m_tmp);
            if m1<=m2
                m_space_RES = cat(1,m_space_RES,m_tmp);
                delta_vec = cat(1,delta_vec,prod(abs(diff(m_tmp))));
            end
        end
        
    case 3
        for m1 = 1:Nr-2
            for m2 = 1:Nr-m1-1
                m3 = Nr-m1-m2;
                m_tmp = [m1,m2,m3];
                m_space_ES = cat(1,m_space_ES,m_tmp);
                if m1<=m2 && m2<=m3
                    m_space_RES = cat(1,m_space_RES,m_tmp);
                    delta_vec = cat(1,delta_vec,prod(abs(diff(m_tmp) + 1)));
                end
            end
        end
        
    case 4
        for m1 = 1:Nr-3
            for m2 = 1:Nr-m1-2
                for m3 = 1:Nr-m1-m2-1
                    m4 = Nr-m1-m2-m3;
                    m_tmp = [m1,m2,m3,m4];
                    m_space_ES = cat(1,m_space_ES,m_tmp);
                    if m1<=m2 && m2<=m3 && m3<=m4
                        m_space_RES = cat(1,m_space_RES,m_tmp);
                        delta_vec = cat(1,delta_vec,prod(abs(diff(m_tmp) + 1)));
                    end
                end
            end
        end
    case 5
        for m1 = 1:Nr-4
            for m2 = 1:Nr-m1-3
                for m3 = 1:Nr-m1-m2-2
                    for m4 = 1:Nr-m1-m2-m3-1
                        m5 = Nr-m1-m2-m3-m4;
                        m_tmp = [m1,m2,m3,m4,m5];
                        m_space_ES = cat(1,m_space_ES,m_tmp);
                        if m1<=m2 && m2<=m3 && m3<=m4 && m4<=m5
                            m_space_RES = cat(1,m_space_RES,m_tmp);
                            delta_vec = cat(1,delta_vec,prod(abs(diff(m_tmp) + 1)));
                        end
                    end
                end
            end
        end
    case 6
        for m1 = 1:4%Nr-5
            for m2 = 1:6%Nr-m1-4
                for m3 = 1:12%Nr-m1-m2-3
                    for m4 = 1:16%Nr-m1-m2-m3-2
                        for m5 = 1:Nr-m1-m2-m3-m4-1
                            m6 = Nr-m1-m2-m3-m4-m5;
                            m_tmp = [m1,m2,m3,m4,m5,m6];
                            m_space_ES = cat(1,m_space_ES,m_tmp);
                            if m1<=m2 && m2<=m3 && m3<=m4 && m4<=m5 && m5<=m6
                                m_space_RES = cat(1,m_space_RES,m_tmp);
                                delta_vec = cat(1,delta_vec,prod(abs(diff(m_tmp) + 1)));
                            end
                        end
                    end
                end
            end
        end
    case 8
        for m1 = 1:4
            for m2 = 1:6
                for m3 = 1:8
                    for m4 = 1:12
                        for m5 = 2:16
                            for m6 = 3:24
                                for m7 = 4:Nr-m1-m2-m3-m4-m5-m6-1
                                    m8 = Nr-m1-m2-m3-m4-m5-m6-m7;
                                    m_tmp = [m1,m2,m3,m4,m5,m6,m7,m8];
                                    m_space_ES = cat(1,m_space_ES,m_tmp);
                                    if m1<=m2 && m2<=m3 && m3<=m4 && m4<=m5...
                                            && m5<=m6 && m6<=m7 && m7<=m8
                                        m_space_RES = cat(1,m_space_RES,m_tmp);
                                        delta_vec = cat(1,delta_vec,prod(abs(diff(m_tmp) + 1)));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
end

[delta_vec_order, i_order] = sort(delta_vec, 'descend');
m_space_RES_EE = m_space_RES(i_order,:);
end