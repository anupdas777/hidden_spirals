function [is_stable, num_pos, num_neutral, max_other] = check_spectrum_stability(eigs, stab_tol, pos_tol, zero_tol, expected_max)
    neutral = abs(eigs) < zero_tol;
    num_neutral = sum(neutral);
    other = eigs(~neutral);

    if isempty(other)
        max_other = -Inf;
        num_pos = 0;
    else
        max_other = max(other);
        num_pos = sum(other > pos_tol);
    end

    is_stable = num_neutral <= expected_max && num_pos == 0 && max_other <= stab_tol;
end

