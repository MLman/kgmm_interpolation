function T = myisreal(X)
    TF_finite = real(isfinite(X));
    TF_finite = prod(TF_finite(:));
    TF_nan = real(~isnan(X));
    TF_nan = prod(TF_nan(:));
    T = isreal(X)&& TF_nan && TF_finite;
end