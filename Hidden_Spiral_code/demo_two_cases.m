this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'funcs'));

out_dir = fullfile(this_dir, 'demo_outputs');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

cfg.mode = 29;
cfg.interp_steps = 300000;
cfg.evolve_steps = 18000;
cfg.interp_dt = 0.025;
cfg.evolve_dt = 0.05;
cfg.convergence_tol = 1e-8;
cfg.lambda_zero_cleanup_steps = 18000;
cfg.homotopy_perturb_sigma = 1e-3;
cfg.checkerboard_tol = 0.05;
cfg.continuation_full_zero_tol = 1e-6;
cfg.lambda_start = 1.0;
cfg.lambda_end = 0.0;
cfg.ds_init = 0.02;
cfg.ds_min = 1e-4;
cfg.ds_max = 0.08;
cfg.max_steps = 600;
cfg.newton_tol = 1e-5;
cfg.newton_max_iter = 10;
cfg.rcond_thresh = 1e-12;
cfg.seed_count = 3;
cfg.seed_dlambda = 0.02;
cfg.seed_tol = 1e-7;
cfg.seed_max_iter = 20;
cfg.jump_dstep_factor = 10;
cfg.stab_tol = 1e-5;
cfg.pos_tol = 1e-6;
cfg.zero_tol = 1e-4;
cfg.expected_neutral_max = 1;

[A_sparse, neighbor_idx, neighbor_count] = build_adjacency(cfg.mode);
optim_options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 150, 'OptimalityTolerance', 5e-4, 'StepTolerance', 1e-5);

cases(1).name = 'continuation_demo';
cases(1).source = 'Encoding / TJ007 / trial 1752';
cases(1).expected_method = 1;
cases(1).phase8 = [ ...
    0.000000000000 3.278542876244 1.272622048855 2.432067751884 2.301191687584 2.472756743431 1.511756032705 2.466531157494; ...
    0.701341986656 1.828285574913 1.952021956444 3.198242068291 1.347963199019 2.834732413292 2.636923551559 2.244858860970; ...
    5.010563556348 1.609973788261 0.571813464165 2.934186339378 5.337545816098 2.366902709007 3.115253448486 2.368743896484; ...
    0.347340166569 1.513491392136 6.084806744252 0.472124457359 5.384988729154 1.428842663765 2.416456341743 1.966600418091; ...
    1.468787759542 2.425541400909 0.440174877644 5.318435136472 5.815491620694 6.144402925168 5.923288408910 0.402638196945; ...
    6.027084652578 5.568407718335 5.374378625547 5.611715261136 0.408202111721 6.173417035733 5.960379425679 5.616914693509; ...
    4.701527539884 4.668299142514 4.107021689415 5.346117679273 5.475628797208 4.914999190961 4.207947611809 3.711585402489; ...
    4.293663382530 3.652017474174 3.650933146477 4.923508826886 0.817004054785 5.061159078275 3.678972125053 3.421165108681];

cases(2).name = 'homotopy_demo';
cases(2).source = 'Encoding / UP021 / trial 956';
cases(2).expected_method = 2;
cases(2).phase8 = [ ...
    0.000000000000 5.982154581939 5.867390666400 3.916676257049 4.381125185882 4.079091523086 3.042434602976 4.043632004653; ...
    5.572228525077 5.791272495185 4.902795527373 3.970339033996 3.921281073486 3.472968488932 2.731422811747 3.091820627451; ...
    4.750219319259 4.844332788383 4.430177185928 4.076894972716 3.863084767257 2.671032339334 2.066346198320 2.080891042948; ...
    5.009368393813 4.691157553588 4.423179839049 4.242075417434 3.755836460982 2.261895924807 1.828173667192 1.789906054735; ...
    5.638535473739 5.081076953803 4.724459383880 4.451061222945 3.902735922729 2.810275703669 1.982541233301 1.763728290796; ...
    5.544913504516 5.391682718192 5.029732678329 4.627364371215 4.059655163680 3.779387209808 4.711497280990 1.304193526506; ...
    5.728028867637 5.558654282485 5.288744423782 4.988025520240 4.664705727492 4.326100562011 4.930108402167 0.541357472539; ...
    5.721482131873 5.792619262134 5.382454131042 5.315390203391 5.192284677421 4.867781851684 4.663597557937 2.015659838915];

results = repmat(struct(), numel(cases), 1);

for i = 1:numel(cases)
    phase8 = cases(i).phase8;
    target29 = interpolate_29(phase8, cfg.interp_steps, cfg.interp_dt, cfg.convergence_tol, cfg.checkerboard_tol);

    anchor_diff = angle(exp(1i * (target29(1:4:end, 1:4:end) - phase8)));
    if max(abs(anchor_diff(:))) > 1e-10
        error('Anchor drift detected in %s.', cases(i).name);
    end

    [c_val, d_val, k_val, opt_temp] = optimize_coupling(target29, optim_options, neighbor_idx, neighbor_count);
    if opt_temp >= cfg.stab_tol
        error('Unstable target in %s.', cases(i).name);
    end

    target_eigs = jacobian_eigs(target29, c_val, d_val, k_val, neighbor_idx, neighbor_count);
    [is_stable, ~, ~, ~] = check_spectrum_stability(target_eigs, cfg.stab_tol, cfg.pos_tol, cfg.zero_tol, cfg.expected_neutral_max);
    if ~is_stable
        error('Target spectrum failed in %s.', cases(i).name);
    end

    method = 0;
    homotopy_record = [];

    [cont_success, branch_jump, cont_record, ~] = run_continuation( ...
        target29, c_val, d_val, k_val, cfg.mode, A_sparse, ...
        cfg.lambda_start, cfg.lambda_end, cfg.ds_init, cfg.ds_min, cfg.ds_max, cfg.max_steps, ...
        cfg.newton_tol, cfg.newton_max_iter, cfg.rcond_thresh, ...
        cfg.seed_count, cfg.seed_dlambda, cfg.seed_tol, cfg.seed_max_iter, ...
        cfg.jump_dstep_factor, ...
        neighbor_idx, neighbor_count, cfg.stab_tol, cfg.pos_tol, cfg.zero_tol, cfg.expected_neutral_max, ...
        cfg.continuation_full_zero_tol, cfg.checkerboard_tol);

    if cont_success && ~branch_jump
        cont_eigs = jacobian_eigs(cont_record(:, :, end), c_val, d_val, k_val, neighbor_idx, neighbor_count);
        [cont_ok, ~, ~, ~] = check_spectrum_stability(cont_eigs, cfg.stab_tol, cfg.pos_tol, cfg.zero_tol, cfg.expected_neutral_max);
        if cont_ok
            method = 1;
            homotopy_record = cont_record;
        end
    end

    if method == 0
        [homo_success, homo_record, ~] = run_homotopy( ...
            target29, c_val, d_val, k_val, cfg.mode, ...
            cfg.stab_tol, cfg.pos_tol, cfg.zero_tol, cfg.expected_neutral_max, ...
            neighbor_idx, neighbor_count, cfg.evolve_steps, cfg.evolve_dt, cfg.convergence_tol, ...
            cfg.lambda_zero_cleanup_steps, cfg.homotopy_perturb_sigma, cfg.checkerboard_tol);

        if homo_success
            homo_eigs = jacobian_eigs(homo_record(:, :, end), c_val, d_val, k_val, neighbor_idx, neighbor_count);
            [homo_ok, ~, ~, ~] = check_spectrum_stability(homo_eigs, cfg.stab_tol, cfg.pos_tol, cfg.zero_tol, cfg.expected_neutral_max);
            if homo_ok
                method = 2;
                homotopy_record = homo_record;
            end
        end
    end

    if method ~= cases(i).expected_method
        error('Method mismatch in %s: got %d, expected %d.', cases(i).name, method, cases(i).expected_method);
    end

    results(i).name = cases(i).name;
    results(i).source = cases(i).source;
    results(i).method = method;
    results(i).phase8 = phase8;
    results(i).target29 = target29;
    results(i).homotopy_record = homotopy_record;
    results(i).c = c_val;
    results(i).d = d_val;
    results(i).k = k_val;
    results(i).opt_temp = opt_temp;

    f = figure('Visible', 'off', 'Position', [100 100 1080 360]);
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    imagesc(phase8);
    axis square;
    colormap('hsv');
    clim([0 2*pi]);
    colorbar;
    title('Original 8x8');

    nexttile;
    imagesc(target29);
    axis square;
    colormap('hsv');
    clim([0 2*pi]);
    colorbar;
    title('Interpolated 29x29');

    nexttile;
    imagesc(homotopy_record(:, :, end));
    axis square;
    colormap('hsv');
    clim([0 2*pi]);
    colorbar;
    if method == 1
        title('Final continuation');
    else
        title('Final homotopy');
    end

    sgtitle(sprintf('%s | %s', cases(i).name, cases(i).source));
    exportgraphics(f, fullfile(out_dir, [cases(i).name '.png']), 'Resolution', 180);
    close(f);
end

save(fullfile(out_dir, 'demo_two_clean_cases.mat'), 'results', 'cfg');
disp(struct2table(rmfield(results, {'phase8', 'target29', 'homotopy_record'})));
