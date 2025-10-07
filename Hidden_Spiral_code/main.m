%% ----------------------------- Configuration -----------------------------
% Output directory
out_dir = 'C:\Users\Gavin\Desktop\presenttaion\Codes\data\Encoding\H\figure';

% Grid size we simulate on (target 15×15)
grid_size = 15;

% Homotopy / simulation controls
num_homotopy_steps = 81;   % p = 1..81
num_steps_per_homotopy = 10000;
dt = 0.1;

% Coupling architecture weights (keep same behavior as your code)
w_second_axis  = 0;  % 2-step axial neighbors
w_diagonals    = 0;  % diagonal neighbors

% -------- Thresholds (requested to be set at the beginning of the script) -
opt_tol = 1.0e-5;    
eig_tol = 1.0e-5;    

% ------------------------------ Input Phase -------------------------------
% Provide raw phase matrix M before running this script (8×8 or 15×15).
% If M is 8×8, we upsample to 15×15 using your sine_interpolation. If it is
% already 15×15, we use it directly.
assert(exist('M','var')==1, 'Please provide raw phase matrix M in workspace.');
if ~ismatrix(M)
    error('M must be a 2D matrix of phases.');
end

switch num2str(size(M))
    case '8  8'
        phase_init = sine_interpolation(M);  
    case '15  15'
        phase_init = M;
    otherwise
        error('M must be 8x8 or 15x15. Got %dx%d.', size(M,1), size(M,2));
end

if size(phase_init,1)~=grid_size || size(phase_init,2)~=grid_size
    error('Interpolated phase is not %dx%d.', grid_size, grid_size);
end

%% ----------------------- Identify coupling parameters --------------------
[coeff_c, coeff_d, coeff_k, opt_val] = minimize_eigenvalue(phase_init);

if opt_val >= opt_tol
    fprintf('Not able to stabilize at the beginning (opt=%.3e).\n', opt_val);
    return
end

hprime = @(x) sin(x) ...
            + coeff_c*(cos(x)-1) ...
            + coeff_d*(cos(2*x)-1) ...
            + coeff_k*sin(2*x);

%% -------------------------- Base frequency field -------------------------
base_omega = build_base_omega(phase_init, hprime, w_second_axis, w_diagonals);

%% -------------------------- Homotopy simulation --------------------------
homotopy_phase = zeros(grid_size, grid_size, num_homotopy_steps);
homotopy_phase(:,:,1) = phase_init;
max_eig_record = zeros(1, num_homotopy_steps);
order_param    = zeros(1, num_homotopy_steps);
for p = 2:num_homotopy_steps
    freq_scale = 1.0125 - 0.0125*p;
    omega_t    = base_omega * freq_scale;
    phases = homotopy_phase(:,:,p-1);
    if p == 2
        phases = phases + 0.01*rand(grid_size);  % your small kick
    end
    for step = 1:num_steps_per_homotopy
        phases = step_phases(phases, omega_t, hprime, dt, w_second_axis, w_diagonals);
    end
    phases = mod(phases - phases(1,1), 2*pi);
    homotopy_phase(:,:,p) = phases;

    eigs_here = Jacobian(phases, coeff_c, coeff_d, coeff_k);
    max_eig_record(p) = eigs_here(end);
    if max_eig_record(p) > eig_tol
        fprintf('Not able to stabilize at homotopy step p=%d (λ_max=%.3e).\n', p, max_eig_record(p));
        return
    end
    dphi = mod(homotopy_phase(:,:,p) - homotopy_phase(:,:,1), 2*pi);
    order_param(p) = compute_order_parameter(dphi);
end

%% ------------------------------- Plotting --------------------------------
final_img = figure('Visible','off'); 
imagesc(homotopy_phase(:,:,num_homotopy_steps));
axis image off
colormap(hsv);
clim([0, 2*pi]);
cb = colorbar; cb.Label.String = 'phase (rad)';
title(sprintf('Homotopy step p=%d', num_homotopy_steps), 'Interpreter','none');
out_heat = fullfile(out_dir, sprintf('homotopy_heatmap_p%02d.png', num_homotopy_steps));
exportgraphics(gcf, out_heat, 'Resolution', 300);
close(gcf);

op_fig = figure('Visible','off'); 
plot(1:num_homotopy_steps, order_param, 'LineWidth', 1.5, 'Marker', '.');
xlabel('homotopy step p');
ylabel('order parameter');
title('Order parameter over homotopy');
grid on
out_op = fullfile(out_dir, 'homotopy_order_parameter.png');
exportgraphics(gcf, out_op, 'Resolution', 300);
close(gcf);
fprintf('Saved:\n  %s\n  %s\n', out_heat, out_op);

%% ============================== Functions ================================
function omega = build_base_omega(P, hprime, w2, wd)
    n = size(P,1);
    omega = zeros(n,n);

    for i = 1:n
        for j = 1:n
            
            up    = P(max(i-1,1), j);
            down  = P(min(i+1,n), j);
            left  = P(i, max(j-1,1));
            right = P(i, min(j+1,n));

          
            up2    = P(max(i-2,1), j);
            down2  = P(min(i+2,n), j);
            left2  = P(i, max(j-2,1));
            right2 = P(i, min(j+2,n));

            
            ul = P(max(i-1,1), max(j-1,1));
            ur = P(max(i-1,1), min(j+1,n));
            dl = P(min(i+1,n), max(j-1,1));
            dr = P(min(i+1,n), min(j+1,n));

         
            self = P(i,j);
            omega(i,j) = ...
               - hprime(left - self)  - hprime(right - self) ...
               - hprime(up   - self)  - hprime(down - self) ...
               - w2*( hprime(up2   - self) + hprime(down2 - self) ...
                    + hprime(left2 - self) + hprime(right2 - self) ) ...
               - wd*( hprime(ur - self) + hprime(dr - self) ...
                    + hprime(dl - self) + hprime(ul - self) );
        end
    end
end

function phases_next = step_phases(phases, omega_t, hprime, dt, w2, wd)
    n = size(phases,1);
    phases_next = phases; 

    for i = 1:n
        for j = 1:n
            self = phases(i,j);

            
            up    = phases(max(i-1,1), j);
            down  = phases(min(i+1,n), j);
            left  = phases(i, max(j-1,1));
            right = phases(i, min(j+1,n));

           
            up2    = phases(max(i-2,1), j);
            down2  = phases(min(i+2,n), j);
            left2  = phases(i, max(j-2,1));
            right2 = phases(i, min(j+2,n));

           
            ul = phases(max(i-1,1), max(j-1,1));
            ur = phases(max(i-1,1), min(j+1,n));
            dl = phases(min(i+1,n), max(j-1,1));
            dr = phases(min(i+1,n), min(j+1,n));

            total_force = ...
                hprime(up    - self) + hprime(down  - self) + ...
                hprime(left  - self) + hprime(right - self) + ...
                w2*( hprime(up2   - self) + hprime(down2 - self) + ...
                    hprime(left2 - self) + hprime(right2 - self) ) + ...
                wd*( hprime(ur - self) + hprime(dr - self) + ...
                    hprime(dl - self) + hprime(ul - self) );

            phases_next(i,j) = self + dt*(omega_t(i,j) + total_force);
        end
    end
end

function r = compute_order_parameter(dphi)
% |(1/N) sum exp(i * Δφ)| with Δφ relative to initial phase
    N = numel(dphi);
    r = abs( sum(exp(1i*dphi), 'all') ) / N;
end

%To check if the flattening process is on the same branch during process,
%simply reverse the process and check final phase with interpolated data
%phase.You may use the order parameter plot as reference. 