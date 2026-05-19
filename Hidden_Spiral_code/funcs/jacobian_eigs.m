function eigs_out = jacobian_eigs(inp_frame, c, d, k, neighbor_idx, neighbor_count)
    mode = size(inp_frame, 1);
    N = mode^2;
    theta = inp_frame(:);

    StaMa = zeros(N, N);

    for i = 1:N
        row_sum = 0;
        for n = 1:neighbor_count(i)
            j = neighbor_idx(i, n);
            diff = theta(j) - theta(i);
            val = cos(diff) - c*sin(diff) - 2*d*sin(2*diff) + 2*k*cos(2*diff);
            StaMa(i, j) = val;
            row_sum = row_sum + val;
        end
        StaMa(i, i) = -row_sum;
    end

    eigs_out = sort(real(eig(StaMa)));
end

