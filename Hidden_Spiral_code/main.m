this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'funcs'));

%% Configuration
% directory = ;                  % output folder
% current_mode = ;               % 15 to 29
% interp_steps = ;               % 1e5 to 1.2e6
% evolve_steps = ;               % 6e3 to 9e4
% interp_dt = ;                  % 0.005 to 0.05
% evolve_dt = ;                  % 0.01 to 0.10
% convergence_tol = ;            % 1e-10 to 1e-6
% lambda_zero_cleanup_steps = ;  % 6e3 to 9e4
% homotopy_perturb_sigma = ;     % 1e-5 to 1e-2
% checkerboard_tol = ;           % 0.02 to 0.10
% continuation_full_zero_tol = ; % 1e-8 to 1e-4
% lambda_start = ;               % 1.0
% lambda_end = ;                 % 0.0
% ds_init = ;                    % 0.005 to 0.05
% ds_min = ;                     % 1e-5 to 1e-3
% ds_max = ;                     % 0.02 to 0.20
% max_steps = ;                  % 200 to 4000
% newton_tol = ;                 % 1e-8 to 1e-4
% newton_max_iter = ;            % 5 to 30
% rcond_thresh = ;               % 1e-14 to 1e-8
% seed_count = ;                 % 2 to 6
% seed_dlambda = ;               % 0.005 to 0.05
% seed_tol = ;                   % 1e-9 to 1e-5
% seed_max_iter = ;              % 5 to 50
% jump_dstep_factor = ;          % 3 to 20
% stab_tol = ;                   % 1e-8 to 1e-4
% pos_tol = ;                    % 1e-8 to 1e-4
% zero_tol = ;                   % 1e-6 to 1e-2
% expected_neutral_max = ;       % 1 to 2

if ~(exist(directory, 'dir') == 7), mkdir(directory); end

[A_sparse, neighbor_idx, neighbor_count] = build_adjacency(current_mode);

%% Trial list
num_trials = size(DataCoordinate, 2);
method_record = zeros(1, num_trials);

trial_selection = unique(trial_selection(:)');
trial_indices = trial_selection(trial_selection >= 1 & trial_selection <= num_trials);

optim_options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 150, 'OptimalityTolerance', 5e-4, 'StepTolerance', 1e-5);

%% Trial data
num_to_process = length(trial_indices);

fprintf('Extracting trial data...\n');
trial_data_cell = cell(num_to_process, 1);
for i = 1:num_to_process
    tri_num = trial_indices(i);
    [trial, start_time, end_time] = decode_trial(tri_num, DataCoordinate);
    frm = floor((start_time + end_time) / 2);
    Phase8 = unwrap_phases(trial, Phases);
    trial_data_cell{i} = struct('tri_num', tri_num, 'Phase8_frame', Phase8(:,:,frm));
end

%% Main loop
fprintf('Processing %d trials serially...\n', num_to_process);
serial_timer = tic;

for idx = 1:num_to_process
    td = trial_data_cell{idx};
    tri_num = td.tri_num;

    fprintf('\n--- Trial %d (%d/%d) ---\n', tri_num, idx, num_to_process);

    result = struct('tri_num', tri_num, 'method', 0, ...
        'eigenvalues', [], 'c', 0, 'd', 0, 'k', 0, 'opt_temp', 0, ...
        'Target_temp', [], 'homotopy_record', [], 'temp_curve', []);

    Phase8_frame = td.Phase8_frame;
    Target_temp = interpolate_29(...
        Phase8_frame, interp_steps, interp_dt, convergence_tol, checkerboard_tol);

    result.Target_temp = Target_temp;

    [c_val, d_val, k_val, opt_temp] = optimize_coupling(...
        Target_temp, optim_options, neighbor_idx, neighbor_count);
    result.c = c_val; result.d = d_val; result.k = k_val; result.opt_temp = opt_temp;

    if opt_temp < stab_tol
        all_eigenvalues = jacobian_eigs(Target_temp, c_val, d_val, k_val, ...
            neighbor_idx, neighbor_count);
        result.eigenvalues = all_eigenvalues;

        [is_stable, ~, ~, ~] = check_spectrum_stability(all_eigenvalues, stab_tol, pos_tol, zero_tol, expected_neutral_max);
        if is_stable
            [cont_success, branch_jump, homotopy_record, temp_curve] = run_continuation(...
                Target_temp, c_val, d_val, k_val, current_mode, A_sparse, ...
                lambda_start, lambda_end, ds_init, ds_min, ds_max, max_steps, ...
                newton_tol, newton_max_iter, rcond_thresh, ...
                seed_count, seed_dlambda, seed_tol, seed_max_iter, ...
                jump_dstep_factor, ...
                neighbor_idx, neighbor_count, stab_tol, pos_tol, zero_tol, expected_neutral_max, ...
                continuation_full_zero_tol, checkerboard_tol);

            if cont_success && ~branch_jump
                final_eig_cont = jacobian_eigs(homotopy_record(:,:,81), c_val, d_val, k_val, ...
                    neighbor_idx, neighbor_count);
                [final_stable_cont, ~, ~, ~] = check_spectrum_stability(final_eig_cont, stab_tol, pos_tol, zero_tol, expected_neutral_max);

                if final_stable_cont
                    result.method = 1;
                    result.homotopy_record = homotopy_record;
                    result.temp_curve = temp_curve;
                end
            end

            if result.method == 0
                [homo_success, homotopy_record, temp_curve] = run_homotopy(...
                    Target_temp, c_val, d_val, k_val, current_mode, ...
                    stab_tol, pos_tol, zero_tol, expected_neutral_max, ...
                    neighbor_idx, neighbor_count, evolve_steps, evolve_dt, convergence_tol, ...
                    lambda_zero_cleanup_steps, homotopy_perturb_sigma, checkerboard_tol);

                if homo_success
                    final_eig_homo = jacobian_eigs(homotopy_record(:,:,81), c_val, d_val, k_val, ...
                        neighbor_idx, neighbor_count);
                    [final_stable_homo, ~, ~, ~] = check_spectrum_stability(final_eig_homo, stab_tol, pos_tol, zero_tol, expected_neutral_max);

                    if final_stable_homo
                        result.method = 2;
                        result.homotopy_record = homotopy_record;
                        result.temp_curve = temp_curve;
                    end
                end
            end
        end
    end

    method_record(tri_num) = result.method;

    method_val = result.method;
    save(fullfile(directory, sprintf('trial_%d_method.mat', tri_num)), 'method_val');

    if ~isempty(result.eigenvalues)
        all_eigenvalues_29 = result.eigenvalues;
        c_29 = result.c; d_29 = result.d; k_29 = result.k;
        opt_temp_29 = result.opt_temp; Target_temp_29 = result.Target_temp;
        save(fullfile(directory, sprintf('trial_%d_eigenvalues_29.mat', tri_num)), ...
            'all_eigenvalues_29', 'c_29', 'd_29', 'k_29', 'opt_temp_29', 'Target_temp_29');
    end

    if result.method > 0
        homotopy_record = result.homotopy_record;
        temp_curve = result.temp_curve;
        mode = current_mode;
        save(fullfile(directory, sprintf('trial_%d_homotopy.mat', tri_num)), ...
            'homotopy_record', 'temp_curve', 'mode');

        method_str = {'FAILED', 'Cont@29', 'Homo@29'};

        f1 = figure('Visible', 'off');
        imagesc(homotopy_record(:,:,81)); colormap('hsv'); clim([0 2*pi]); colorbar;
        title(sprintf('Trial %d, %s', tri_num, method_str{result.method+1}));
        print(f1, fullfile(directory, sprintf('homotopy_image--%d.jpg', tri_num)), '-djpeg', '-r150');
        savefig(f1, fullfile(directory, sprintf('homotopy_image--%d.fig', tri_num)));
        close(f1);

        f2 = figure('Visible', 'off');
        scatter(1:81, temp_curve, 5);
        xlabel('Step'); ylabel('Order parameter');
        title(sprintf('Trial %d, %s', tri_num, method_str{result.method+1}));
        print(f2, fullfile(directory, sprintf('homotopy_scatter--%d.jpg', tri_num)), '-djpeg', '-r150');
        savefig(f2, fullfile(directory, sprintf('homotopy_scatter--%d.fig', tri_num)));
        close(f2);
    end

    save(fullfile(directory, 'method_record.mat'), 'method_record');

    fprintf('Trial %d complete. Running totals - Cont: %d | Homo: %d | Failed: %d\n', ...
        tri_num, sum(method_record==1), sum(method_record==2), sum(method_record==0));
end

elapsed_time = toc(serial_timer);
fprintf('\n=== Final Summary ===\n');
fprintf('Continuation: %d | Homotopy: %d | Failed: %d\n', ...
    sum(method_record==1), sum(method_record==2), sum(method_record==0));
fprintf('Total time: %.1f sec (%.2f sec/trial)\n', elapsed_time, elapsed_time/num_to_process);

