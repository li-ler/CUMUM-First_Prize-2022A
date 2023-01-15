function v = V(xf)
    r = 1;
    v = 0;
    if xf >= 0.8
        v = 0;
    elseif xf > 0
        v = pi * (r ^ 2) * ((0.8 - xf) ^ 3) / (0.8 ^ 2) / 3;
    elseif xf > -3
        v = 0.8 * pi / 3 - pi * (r ^ 2) * xf;
    else
        v = 9.8 * pi / 3;
    end
end

