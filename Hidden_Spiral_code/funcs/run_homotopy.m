function [success, homotopy_record, temp_curve] = run_homotopy(...
    Target_temp, c, d, k, mode, stab_tol, pos_tol, zero_tol, expected_max, ...
    neighbor_idx, neighbor_count, num_steps, dt, conv_tol, cleanup_steps, perturb_sigma, checkerboard_tol)
    success = false;
    homotopy_record = zeros(mode, mode, 81);
    temp_curve = zeros(1, 81);

    [success, homotopy_record, temp_curve] = run_homotopy_attempt(...
        Target_temp, c, d, k, mode, stab_tol, pos_tol, zero_tol, expected_max, ...
        neighbor_idx, neighbor_count, num_steps, dt, conv_tol, cleanup_steps, perturb_sigma);

    if success && max_checkerboard_record(homotopy_record) > checkerboard_tol
        success = false;
        homotopy_record = zeros(mode, mode, 81);
        temp_curve = zeros(1, 81);
    end
end

function [success, homotopy_record, temp_curve] = run_homotopy_attempt(...
    Target_temp, c, d, k, mode, stab_tol, pos_tol, zero_tol, expected_max, ...
    neighbor_idx, neighbor_count, num_steps, dt, conv_tol, cleanup_steps, perturb_sigma)
    success = false;
    homotopy_record = zeros(mode, mode, 81);
    temp_curve = zeros(1, 81);

    dhp = @(x) sin(x) + c * (cos(x) - 1) + d * (cos(2 * x) - 1) + k * sin(2 * x);
    w = natural_frequency(Target_temp, dhp);

    homotopy_record(:, :, 1) = Target_temp;
    temp_curve(1) = 1.0;

    for p = 2:81
        phases = homotopy_record(:, :, p - 1) + perturb_sigma * rand(mode);
        scale = 1.0125 - 0.0125 * p;

        phases = evolve_phases(phases, w * scale, dhp, num_steps, dt, conv_tol);
        if scale == 0 && cleanup_steps > 0
            phases = evolve_phases(phases, zeros(mode), dhp, cleanup_steps, dt, conv_tol);
        end
        homotopy_record(:, :, p) = mod(phases - phases(1, 1), 2 * pi);

        homo_eig = jacobian_eigs(homotopy_record(:, :, p), c, d, k, neighbor_idx, neighbor_count);
        [is_stable, ~, ~, ~] = check_spectrum_stability(homo_eig, stab_tol, pos_tol, zero_tol, expected_max);
        if ~is_stable
            return;
        end

        diff_phase = mod(homotopy_record(:, :, p) - homotopy_record(:, :, 1), 2 * pi);
        temp_curve(p) = abs(mean(exp(1i * diff_phase(:))));
    end

    success = true;
end

function phases = evolve_phases(phases_init, w, dhp, num_steps, dt, conv_tol)
    phases = phases_init;

    for s = 1:num_steps
        phases_old = phases;

        up = circshift(phases, [1, 0]); up(1, :) = phases(1, :);
        dn = circshift(phases, [-1, 0]); dn(end, :) = phases(end, :);
        lt = circshift(phases, [0, 1]); lt(:, 1) = phases(:, 1);
        rt = circshift(phases, [0, -1]); rt(:, end) = phases(:, end);

        force = dhp(up - phases) + dhp(dn - phases) + dhp(lt - phases) + dhp(rt - phases);
        phases = phases + dt * (w + force);

        if max(abs(phases(:) - phases_old(:))) < conv_tol
            break;
        end
    end

    phases = mod(phases - phases(1, 1), 2 * pi);
end

function w = natural_frequency(Target, dhp)
    up = circshift(Target, [1, 0]); up(1, :) = Target(1, :);
    dn = circshift(Target, [-1, 0]); dn(end, :) = Target(end, :);
    lt = circshift(Target, [0, 1]); lt(:, 1) = Target(:, 1);
    rt = circshift(Target, [0, -1]); rt(:, end) = Target(:, end);
    w = -dhp(lt - Target) - dhp(rt - Target) - dhp(up - Target) - dhp(dn - Target);
end

