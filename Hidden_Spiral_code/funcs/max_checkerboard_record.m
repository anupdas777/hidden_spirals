function max_score = max_checkerboard_record(record)
    max_score = 0;
    [ii, jj] = ndgrid(1:size(record, 1), 1:size(record, 2));
    alt = ones(size(record, 1), size(record, 2));
    alt(mod(ii + jj, 2) == 1) = -1;

    for kk = 1:size(record, 3)
        frame = record(:, :, kk);
        score = abs(mean(exp(1i * frame(:)) .* alt(:)));
        max_score = max(max_score, score);
    end
end

