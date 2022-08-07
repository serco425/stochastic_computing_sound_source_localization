function [param] = calcMismatch(gain_std, d_std, M, param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rows = length(param.f);

if gain_std ~= 0
    diff = gain_std;
    while diff > gain_std/10
        gain = gain_std*randn(M,1);
        diff = abs(gain_std-std(gain));
    end

    param.gain_dev = NaN(rows,M);

    for channel = 1:M;
        param.gain_dev(:,channel) = gain(channel)*ones(rows,1);
    end
else
    param.gain_dev = zeros(rows,M);
end

if d_std ~= 0
    diff = d_std;
    while diff > d_std/10
        d = d_std*randn(M,1);
        diff = abs(d_std-std(d));
    end

    rows = length(param.f);
    param.place_dev = NaN(rows,M);

    for channel = 1:M;
        param.place_dev(:,channel) = d(channel)*ones(rows,1);
    end
else
    param.place_dev = zeros(rows,M);
end