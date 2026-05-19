function [trial, start_time, end_time] = decode_trial(x, DataCoordinate)
    trial = DataCoordinate(1,x);
    start_time = DataCoordinate(2,x);
    end_time = DataCoordinate(3,x);
end

