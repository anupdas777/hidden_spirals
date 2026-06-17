function [c, d, k, min_eig] = optimize_coupling(A, options, neighbor_idx, neighbor_count)
    theta = A(:);
    obj = @(p) spectral_objective(theta, p, neighbor_idx, neighbor_count);

    [params, min_eig] = fmincon(obj, [0 0 0], [], [], [], [], ...
        [-0.75 -0.75 -0.75], [0.75 0.75 0.75], [], options);
    c = params(1); d = params(2); k = params(3);
end

function f = spectral_objective(theta, params, neighbor_idx, neighbor_count)
    c = params(1); d = params(2); k = params(3);
    N = numel(theta);
    StaMa = zeros(N, N);

    for i = 1:N
        row_sum = 0;
        for n = 1:neighbor_count(i)
            j = neighbor_idx(i, n);
            diff = theta(j) - theta(i);
            val = cos(diff) - c * sin(diff) - 2 * d * sin(2 * diff) + 2 * k * cos(2 * diff);
            StaMa(i, j) = val;
            row_sum = row_sum + val;
        end
        StaMa(i, i) = -row_sum;
    end

    opts_eig.tol = 1e-6;
    opts_eig.maxit = 200;
    opts_eig.disp = 0;
    f = real(eigs(sparse(StaMa), 1, 'largestreal', opts_eig));

    if ~isfinite(f)
        f = 1e10;
    end
end

