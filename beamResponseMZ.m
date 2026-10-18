function [v_sh, M] = beamResponseMZ(x, L_sh, L, EI, forcePos, forceVal, momentPos, momentVal, q_sup, L_g, delta)
% Constant EI; upward forces/deflection positive, couples CCW positive.
% v(0)=0 and v(L_g)=delta [m]. The caller may pass L_d as L_g.
% Supply the reactions from the corresponding pinned/brace model.
    if nargin < 11, delta = 0; end
    assert(L_g > 0 && L_g <= L, 'Second displacement location must be in (0,L].');
    M = q_sup.*(L-x).^2/2;
    for k = 1:numel(forceVal)
        M = M + forceVal(k).*max(forcePos(k)-x, 0);
    end
    for k = 1:numel(momentVal)
        M = M - momentVal(k).*(x <= momentPos(k));
    end
    P0 = particular(0);
    C1 = (delta+P0-particular(L_g))/L_g;
    v_sh = particular(L_sh)-P0+C1*L_sh;

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
