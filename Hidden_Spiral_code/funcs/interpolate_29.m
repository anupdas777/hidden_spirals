function ready = interpolate_29(Phase, num_steps, dt, conv_tol, checkerboard_tol)
    target_mode = 29;
    phases = Phase;

    while size(phases, 1) < target_mode
        [phases_next, converged] = interpolate_stage(phases, num_steps, dt, conv_tol);
        if ~converged || checkerboard_score(phases_next) > checkerboard_tol
            error('Interpolation stage %d->%d did not converge cleanly.', size(phases, 1), 2*size(phases, 1) - 1);
        end

        phases = phases_next;
    end

    if size(phases, 1) ~= target_mode
        error('Interpolation size mismatch: %d vs %d.', size(phases, 1), target_mode);
    end

    ready = mod(phases - phases(1,1), 2*pi);
end

function [phases_next, converged] = interpolate_stage(phases, num_steps, dt, conv_tol)
    M = 2 * size(phases, 1) - 1;
    phases_next = zeros(M);
    phases_next(1:2:M, 1:2:M) = phases;

    fixed_mask = false(M);
    fixed_mask(1:2:M, 1:2:M) = true;
    fixed_vals = phases_next;
    converged = false;

    for s = 1:num_steps
        phases_old = phases_next;
        up = circshift(phases_next, [1, 0]); up(1, :) = phases_next(1, :);
        dn = circshift(phases_next, [-1, 0]); dn(end, :) = phases_next(end, :);
        lt = circshift(phases_next, [0, 1]); lt(:, 1) = phases_next(:, 1);
        rt = circshift(phases_next, [0, -1]); rt(:, end) = phases_next(:, end);

        phases_next = phases_next + dt * ( ...
            sin(up - phases_next) + sin(dn - phases_next) + ...
            sin(lt - phases_next) + sin(rt - phases_next));
        phases_next(fixed_mask) = fixed_vals(fixed_mask);

        if max(abs(phases_next(:) - phases_old(:))) < conv_tol
            converged = true;
            break;
        end
    end
end

function score = checkerboard_score(theta)
    [ii, jj] = ndgrid(1:size(theta, 1), 1:size(theta, 2));
    alt = ones(size(theta));
    alt(mod(ii + jj, 2) == 1) = -1;
    score = abs(mean(exp(1i * theta(:)) .* alt(:)));
end

