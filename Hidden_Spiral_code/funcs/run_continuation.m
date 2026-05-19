function [success, branch_jump, homotopy_record, temp_curve] = run_continuation(...
    Target_temp, c, d, k, mode, A, lambda_start, lambda_end, ds_init, ds_min, ds_max, ...
    max_steps, newton_tol, newton_max_iter, rcond_thresh, ...
    seed_count, seed_dlambda, seed_tol, seed_max_iter, ...
    jump_dstep_factor, ...
    neighbor_idx, neighbor_count, stab_tol, pos_tol, zero_tol, expected_neutral_max, ...
    continuation_full_zero_tol, checkerboard_tol)
    success = false; branch_jump = false;
    homotopy_record = zeros(mode, mode, 81);
    temp_curve = zeros(1, 81);
    N = mode^2;
    theta0 = Target_temp(:);
    theta_ref = theta0(1);

    coupling0 = coupling_force(theta0, A, c, d, k);
    omega = -coupling0;

    branch_theta = zeros(N, max_steps);
    branch_lambda = zeros(max_steps, 1);
    branch_theta(:,1) = theta0;
    branch_lambda(1) = lambda_start;
    homotopy_record(:,:,1) = reshape(theta0, mode, mode);
    temp_curve(1) = abs(mean(exp(1i*theta0)));

    actual_seeds = 1;
    for s = 2:seed_count
        lambda_t = max(lambda_end, lambda_start - (s-1)*seed_dlambda);
        theta_g = branch_theta(:, s-1);
        [theta_s, ok] = newton_solve(theta_g, omega, lambda_t, A, theta_ref, c, d, k, seed_tol, seed_max_iter, rcond_thresh);
        if ~ok, break; end
        actual_seeds = s;
        branch_theta(:,s) = theta_s;
        branch_lambda(s) = lambda_t;
    end
    if actual_seeds < 2, return; end
    
    step = actual_seeds;
    ds = ds_init;

    while step < max_steps && branch_lambda(step) > lambda_end
        t_theta = branch_theta(:,step) - branch_theta(:,step-1);
        t_lambda = branch_lambda(step) - branch_lambda(step-1);
        tnorm = sqrt(norm(t_theta)^2 + t_lambda^2);
        t_theta = t_theta/tnorm; t_lambda = t_lambda/tnorm;

        theta_pred = branch_theta(:,step) + ds*t_theta;
        lambda_pred = max(lambda_end, min(lambda_start, branch_lambda(step) + ds*t_lambda));
        theta_pred(1) = theta_ref;

        theta_c = theta_pred; lambda_c = lambda_pred;
        converged = false;

        for nit = 1:newton_max_iter
            F = system_residual(theta_c, omega, lambda_c, A, theta_ref, c, d, k);
            g = sum((theta_c - branch_theta(:,step)).*t_theta) + (lambda_c - branch_lambda(step))*t_lambda - ds;
            res = [F; g];
            if norm(res) < newton_tol, converged = true; break; end
            
            J = jacobian_matrix(theta_c, A, c, d, k);
            dF_dl = omega; dF_dl(1) = 0;
            J_ext = [J, dF_dl; t_theta', t_lambda];
            
            if rcond(full(J_ext)) < rcond_thresh
                delta = -pinv(full(J_ext))*res;
            else
                delta = -J_ext\res;
            end
            
            theta_c = theta_c + delta(1:N);
            lambda_c = max(lambda_end, min(lambda_start, lambda_c + delta(N+1)));
            theta_c(1) = theta_ref;
        end
        
        if ~converged
            ds = ds/2;
            if ds < ds_min, return; end
            continue;
        end

        d_step = norm(wrap_to_pi(theta_c - theta_c(1) - (branch_theta(:,step) - branch_theta(1,step))))/sqrt(N);
        if d_step > jump_dstep_factor*ds, branch_jump = true; return; end

        step = step + 1;
        branch_theta(:,step) = theta_c;
        branch_lambda(step) = lambda_c;

        lambda_prev = branch_lambda(step-1);
        check_thresholds = [0.75, 0.5, 0.25, 0.1];
        for thresh = check_thresholds
            if lambda_prev > thresh && lambda_c <= thresh
                step_frame = reshape(theta_c, mode, mode);
                step_eig = jacobian_eigs(step_frame, c, d, k, neighbor_idx, neighbor_count);
                [step_stable, ~, ~, ~] = check_spectrum_stability(step_eig, stab_tol, pos_tol, zero_tol, expected_neutral_max);
                if ~step_stable
                    return;
                end
                break;
            end
        end

        if nit <= 3, ds = min(ds*1.3, ds_max);
        elseif nit >= 7, ds = max(ds/1.5, ds_min); end

        if lambda_c <= lambda_end + 1e-12, break; end
    end

    if branch_lambda(step) > lambda_end + continuation_full_zero_tol
        success = false;
        return;
    end

    lambda_targets = linspace(lambda_start, lambda_end, 81);
    lambda_hist = branch_lambda(1:step);
    theta_hist = branch_theta(:, 1:step);
    for p = 1:81
        idx = find(lambda_hist <= lambda_targets(p), 1, 'first');
        if isempty(idx), idx = length(lambda_hist); end
        if idx == 1
            th = theta_hist(:,1);
        else
            a = (lambda_targets(p) - lambda_hist(idx-1))/(lambda_hist(idx) - lambda_hist(idx-1) + 1e-12);
            th = theta_hist(:,idx-1) + a*(theta_hist(:,idx) - theta_hist(:,idx-1));
        end
        homotopy_record(:,:,p) = mod(reshape(th, mode, mode), 2*pi);
        diff_phase = mod(homotopy_record(:,:,p) - homotopy_record(:,:,1), 2*pi);
        temp_curve(p) = abs(mean(exp(1i*diff_phase(:))));
    end

    check_points = [21, 41, 61, 81];
    for cp = check_points
        interp_eig = jacobian_eigs(homotopy_record(:,:,cp), c, d, k, neighbor_idx, neighbor_count);
        [interp_stable, ~, ~, ~] = check_spectrum_stability(interp_eig, stab_tol, pos_tol, zero_tol, expected_neutral_max);
        if ~interp_stable
            success = false;
            return;
        end
    end

    if max_checkerboard_record(homotopy_record) > checkerboard_tol
        success = false;
        return;
    end

    success = true;
end

function [theta_sol, ok] = newton_solve(theta, omega, lambda, A, theta_ref, c, d, k, tol, max_iter, rcond_thresh)
    ok = false;
    N = length(theta);

    for it = 1:max_iter
        F = system_residual(theta, omega, lambda, A, theta_ref, c, d, k);
        if norm(F) < tol
            ok = true;
            break;
        end

        J = jacobian_matrix(theta, A, c, d, k);
        if rcond(full(J)) < rcond_thresh
            delta = -pinv(full(J)) * F;
        else
            delta = -J \ F;
        end

        theta = theta + delta(1:N);
        theta(1) = theta_ref;
    end

    theta_sol = theta;
end

function F = system_residual(theta, omega, lambda, A, theta_ref, c, d, k)
    coupling = coupling_force(theta, A, c, d, k);
    F = lambda * omega + coupling;
    F(1) = theta(1) - theta_ref;
end

function J = jacobian_matrix(theta, A, c, d, k)
    N = length(theta);
    J = sparse(N, N);
    hp = @(x) cos(x) - c * sin(x) - 2 * d * sin(2 * x) + 2 * k * cos(2 * x);
    [rows, cols] = find(A);

    for e = 1:length(rows)
        i = rows(e);
        j = cols(e);
        val = hp(theta(j) - theta(i));
        J(i, j) = J(i, j) + val;
        J(i, i) = J(i, i) - val;
    end

    J(1, :) = 0;
    J(1, 1) = 1;
end

function force = coupling_force(theta, A, c, d, k)
    N = length(theta);
    force = zeros(N, 1);
    h = @(x) sin(x) + c * (cos(x) - 1) + d * (cos(2 * x) - 1) + k * sin(2 * x);
    [rows, cols] = find(A);

    for e = 1:length(rows)
        i = rows(e);
        j = cols(e);
        force(i) = force(i) + h(theta(j) - theta(i));
    end
end

function y = wrap_to_pi(x)
    y = mod(x + pi, 2 * pi) - pi;
end

