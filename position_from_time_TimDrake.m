function [radius, true_anomaly, eccentric_anomaly] = position_from_time_TimDrake(mu, semimajor, eccen, time)
% Tim Drake Astrodynamics Project 2
% Computes radius, true anomaly, and eccentric anomaly.
% INPUTS:  
%   mu = [km^3/s^2] gravitational constant of the Central body
%   semimajor = [km] the semimajor axis of the elliptical orbit
%   eccen = the orbit eccentricity of the elliptical orbit
%   time = [seconds] the travel time from periapse to the given radius (where inbound values are negative)
%          (can be a column vector)
% OUTPUTS: 
%   radius = [km] the radial distance from the central body at the input time.
%   true_anomaly = [rad] the true anomaly of the satellite at the given position (where inbound values are negative)
%   eccentric_anomaly = [rad] the true anomaly of the satellite at the given position (where inbound values are negative)

M = sqrt(mu ./ (semimajor .^ 3)) .* time;
eccentric_anomaly = getU(M, eccen, M);
true_anomaly = 2 * atan(sqrt((1 + eccen) / (1 - eccen)) * tan(eccentric_anomaly / 2));
radius = (semimajor .* (1 - eccen.^2)) ./ (1 + (eccen .* cos(true_anomaly)));
end

function u_new = getU(u, eccen, M)
% Computes the eccentric anomaly using recursion.
% INPUTS: u = [rad] current eccentric anomaly
%     eccen = the orbit eccentricity of the elliptical orbit
%         M = [rad] mean anomaly
% OUTPUT: u_new : [rad] new eccentric anomaly
F = u - (eccen .* sin(u)) - M;
F_p = 1 - (eccen .* cos(u));
u_new = u - (F ./ F_p);
if abs(F) <= 10^-5
    u_new = u; return
end
u_new = getU(u_new, eccen, M);
end

