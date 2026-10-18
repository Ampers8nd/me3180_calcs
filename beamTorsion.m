function [theta, tau_max, T_internal, J] = beamTorsion(x_locs, x_theta, torquePos, torqueVal, width, height, G, th_support)
% Uniform closed rectangular tube; thin-wall Saint-Venant model.
% Fixed at x=0. Positive signed torques use the right-hand rule about +x.
% theta is evaluated once at x_theta; remaining arrays follow x_locs.
    assert(isscalar(x_theta) && x_theta >= 0, 'x_theta must be nonnegative and scalar.');
    assert(numel(torquePos) == numel(torqueVal), 'Torque arrays must have matching lengths.');
    assert(all([width height G th_support] > 0), 'Dimensions and G must be positive.');
    assert(2*th_support < min(width,height), 'Thickness must leave a hollow interior.');
    A_mid = (width-th_support)*(height-th_support);
    perimeter_mid = 2*(width+height-2*th_support);
    J = 4*A_mid^2*th_support/perimeter_mid;
    theta = 0;
    T_internal = zeros(size(x_locs));
    for k = 1:numel(torqueVal)
        theta = theta + torqueVal(k)*min(x_theta,torquePos(k))/(G*J);
        T_internal = T_internal + torqueVal(k).*(x_locs < torquePos(k));
    end
    tau_max = abs(T_internal)/(2*A_mid*th_support);
end
