function [A, neighbor_idx, neighbor_count] = build_adjacency(mode)
    N = mode^2;
    A = sparse(N, N);
    neighbor_idx = zeros(N, 4);
    neighbor_count = zeros(N, 1);
    
    for i = 1:mode
        for j = 1:mode
            idx = (i-1)*mode + j;
            cnt = 0;
            if i > 1
                nb = (i-2)*mode + j;
                A(idx, nb) = 1;
                cnt = cnt + 1;
                neighbor_idx(idx, cnt) = nb;
            end
            if i < mode
                nb = i*mode + j;
                A(idx, nb) = 1;
                cnt = cnt + 1;
                neighbor_idx(idx, cnt) = nb;
            end
            if j > 1
                nb = (i-1)*mode + (j-1);
                A(idx, nb) = 1;
                cnt = cnt + 1;
                neighbor_idx(idx, cnt) = nb;
            end
            if j < mode
                nb = (i-1)*mode + (j+1);
                A(idx, nb) = 1;
                cnt = cnt + 1;
                neighbor_idx(idx, cnt) = nb;
            end
            neighbor_count(idx) = cnt;
        end
    end
end

