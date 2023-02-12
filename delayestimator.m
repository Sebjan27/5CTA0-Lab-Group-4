function tau = delayestimator(x,y)
    l = length(x);
    r = xcorr(x,y,l);
    [r_max tau] = max(r);
    tau = tau - l - 1;
end