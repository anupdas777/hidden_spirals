function ready = sine_interpolation(Phase)
% SINE_INTERPOLATION  Upsamples an N×N phase grid to (2N–1)×(2N–1) and evolves it
%   ready = sine_interpolation(Phase)
%
%   - Phase: an N×N matrix of initial phases.
%   - The function upsamples by nearest-neighbor to size M = 2N–1,
%     then integrates nearest-neighbor sin-coupling for num_steps iterations.
%   - Returns the final phase relative to the (1,1) node, mod 2π.

%—— PARAMETERS ——
num_steps = 20000;
dt        = 0.1;

%—— 1) Determine input size and output grid size ——
[Nr, Nc] = size(Phase);
assert(Nr == Nc, 'Phase must be square');
Ni   = Nr;
M    = 2*Ni - 1;       % output grid dimension

%—— 2) Upsample by nearest-neighbor ——
[r_raw, c_raw, v_raw] = find(Phase);
r_norm = (r_raw - 1) / (Ni - 1);
c_norm = (c_raw - 1) / (Ni - 1);
r_up   = round(r_norm * (M - 1)) + 1;
c_up   = round(c_norm * (M - 1)) + 1;

Phase_up = zeros(M);
for k = 1:numel(r_up)
    Phase_up(r_up(k), c_up(k)) = v_raw(k);
end

%—— 3) Mark “fixed” nodes to re-impose each step ——
ini_up = Phase_up;
ini_up(1,1) = Phase_up(1,1);  % ensure the (1,1) node stays fixed
[fixed_r, fixed_c] = find(ini_up);

%—— 4) Pre-allocate arrays ——
phases       = Phase_up;
newphase     = zeros(M, M);
phase_record = zeros(M, M, num_steps);

%—— 5) Define coupling function ——
dhp = @(x) sin(x);

%—— 6) Time-stepping loop ——
for step = 1:num_steps
    for i = 1:M
        for j = 1:M
            up    = phases(max(i-1,1), j);
            down  = phases(min(i+1,M), j);
            left  = phases(i, max(j-1,1));
            right = phases(i, min(j+1,M));
            
            totalForce = dhp(up    - phases(i,j)) + ...
                         dhp(down  - phases(i,j)) + ...
                         dhp(left  - phases(i,j)) + ...
                         dhp(right - phases(i,j));
            
            newphase(i,j) = phases(i,j) + dt * totalForce;
        end
    end
    
    % Re-impose all fixed nodes explicitly
    for k = 1:numel(fixed_r)
        newphase(fixed_r(k), fixed_c(k)) = Phase_up(fixed_r(k), fixed_c(k));
    end
    
    phases(:,:,1) = [];  % clear old if needed (not strictly necessary)
    phases        = newphase;
    phase_record(:,:,step) = phases;
end

%—— 7) Compute relative phase at final time ——
final_grid = phase_record(:,:,end);
ready = mod( final_grid - final_grid(1,1), 2*pi );
end
