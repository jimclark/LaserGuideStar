function [E] = ecc_from_mean(M,e)
    format long
    E = M;
    k = 1;
    err = 1e-10;
    while (max(k)>err)
        y = E - e*sin(E) - M;
        dy = 1 - e*cos(E);
        k = abs(y./dy);
        E = E-(y./dy);
    end
end