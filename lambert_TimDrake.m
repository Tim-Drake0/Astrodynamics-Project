function [semimajor, eccen, true_anomaly_1, true_anomaly_2] = lambert_TimDrake(mu, r1, r2, theta, time)
% Tim Drake Astrodynamics Project 3
% Computes semimajor, eccentricity, and 2 true anomaly's of an orbit.
% INPUTS:
%   mu = [km^3/s^2] gravitational constant of the Central body
%   r1 = [km] the radial distance from the central body at the initial time
%   r2 = [km] the radial distance from the central body at the final time
%   theta = [degrees] the angle between positions 1 and 2 in radians; if theta > π, then the desired trajectory goes the "long" way around the central body
%   time = [seconds] the travel time from position 1 to position 2. Can be
%          scalar or vector.
% OUTPUTS:
%   semimajor = [km] the semimajor axis of the elliptical orbit that satisfies the Lambert problem
%   eccen = the orbit eccentricity of the elliptical orbit that satisfies the Lambert problem
%   true_anomaly_1 = [rads] the true anomaly of the satellite in radians at position 1 (where inbound values are negative)
%   true_anomaly_2 = [rads] the true anomaly of the satellite in radians at position 2 (where inbound values are negative)

c = sqrt(r1^2 + r2^2 - 2 * r1 * r2 * cos(theta));
s = .5 * (r1 + r2 + c);
a_m = .5 * s;
beta_m = 2 * asin(sqrt((s - c) / (2 * a_m)));
t_m = sqrt(a_m^3/mu) * (pi - beta_m + sin(beta_m));
% Check if on upper or lower
if time > t_m
    low = a_m;
    high = 200 * a_m;
else
    low = 200 * a_m;
    high = a_m;
end
eVal = 100;
% Bisection method to find semimajor
while (abs(eVal) > 1e-6)
    a_guess = (high + low)/2;
    alpha = 2 * asin(sqrt(s ./ (2 * a_guess)));
    beta = 2 * asin(sqrt((s - c) ./ (2 * a_guess)));
    if time > t_m
        alpha = (2 * pi) - alpha;
    end
    if theta > pi
        beta = -1 * beta;
    end
    time_calculated = sqrt(a_guess.^3/mu) .* (alpha - beta - (sin(alpha) - sin(beta)));
    eVal = time_calculated - time;
    if (eVal > 0)
        % replace the high value
        high = a_guess;
    else
        low = a_guess;
    end
end
semimajor = a_guess;
p = 4 * ((semimajor .* (s - r1) .* (s - r2)) ./ c^2) .* (sin((alpha + beta) ./ 2).^2);
eccen = sqrt(1 - (p ./ semimajor));
true_anomaly_1 = acos((1 ./ eccen) .* (((semimajor .* (1 - eccen.^2)) ./ r1) - 1));
true_anomaly_2 = true_anomaly_1 + theta;
% keep it between 0 and 2pi
if true_anomaly_2 - 2*pi > 0
    true_anomaly_2 = true_anomaly_2 - 2*pi;
end