function [v_sh, M] = beamResponse(x, L_sh, L, EI, forcePos, forceVal, momentPos, momentVal, q_sup)
% Constant EI. x is positive leftward from the fixed support.
% Forces/deflection upward positive; applied couples CCW positive.
% Cantilever conditions: v(0) = 0, v'(0) = 0. M is sagging-positive.
    M = q_sup .* (L-x).^2/2;
    for k = 1:numel(forceVal)
        M = M + forceVal(k).*max(forcePos(k)-x, 0);
    end
    for k = 1:numel(momentVal)
        M = M - momentVal(k).*(x <= momentPos(k));
    end
    P0 = particular(0);
    dP0 = -q_sup*L^3/(6*EI);
    for k = 1:numel(forceVal)
        dP0 = dP0 - forceVal(k)*forcePos(k)^2/(2*EI);
    end
    for k = 1:numel(momentVal)
        dP0 = dP0 + momentVal(k)*momentPos(k)/EI;
    end
    v_sh = particular(L_sh)-P0-dP0*L_sh;

    function P = particular(z)
        P = q_sup.*(L-z).^4/(24*EI);
        for j = 1:numel(forceVal)
            P = P + forceVal(j).*max(forcePos(j)-z, 0).^3/(6*EI);
        end
        for j = 1:numel(momentVal)
            P = P - momentVal(j).*max(momentPos(j)-z, 0).^2/(2*EI);
        end
    end
end
