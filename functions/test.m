for n = 1:size (Y,2);
    y = Y(:,n);
    bmin= prctile (y, 0.05);
    [b1, sn1] = estimate_baseline_noise(y, bmin);
    b (:,n)= b1(:)';
    sn (:,n)= sn1(:)';
end