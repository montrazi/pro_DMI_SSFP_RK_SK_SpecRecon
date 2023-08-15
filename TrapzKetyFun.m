function ct = TrapzKetyFun(Cp,pk,t)
    kT = pk(1);
    kep = pk(2);

    stepSize = t(2)-t(1);

    % Initialize output vector
    ct = zeros(size(Cp));

    % Calculate concentration at each time point
    for k=1:length(t)
        ct(k) = kT*trapz(Cp(1:k).*exp(kep*(t(1:k)-t(k))));
    end

    % Scale by temporal resolution, and transpose (transpose is personal preference)
    ct = stepSize * ct';
end
