function Phase_Trialt = unwrap_phases(trial_number, Phases)
    Phase_Trialt = squeeze(Phases(:,:,trial_number,:) - Phases(1,1,trial_number,:));
    temp_size = size(Phases);
    reshaped = reshape(Phase_Trialt, 64, temp_size(4));
    for i = 1:64
        reshaped(i,:) = unwrap(reshaped(i,:));
    end
    Phase_Trialt = reshape(reshaped, 8, 8, temp_size(4));
end

