function [time, true_anomaly, eccentric_anomaly] = time_from_position_TimDrake(mu, a, eccen, radius, is_outbound)
% Tim Drake Astrodynamics Project 1
% Computes true and eccentric anomaly  and the time from periapse for the given inputs.
% INPUTS:  
%   mu = [km^3/s^2] gravitational constant of the Central body
%   a = [km] semimajor axis of the elliptical orbit
%   eccen = the orbit eccentricity of the elliptical orbit
%   radius = [km] the radial distance from the central body at a given
%            time (Can be an array to compute multiple at once)
%   is_outbound = [bool] true if the orbit is moving away from periapse, false if it
%                 is moving towards periapse (If computing multple radii, this is an array
%                 of T/F values)
% OUTPUTS:
%   time = [seconds] the travel time from periapse to the given radius (where inbound values are negative)
%   true_anomaly = [rad] the true anomaly of the satellite at the given position (where inbound values are negative)
%   eccentric_anomaly = [rad] the true anomaly of the satellite at the given position (where inbound values are negative)
true_anomaly = acos((1 ./ eccen) * (((a .* (1 - (eccen .^ 2))) ./ radius) - 1));
if ~is_outbound % f is negative for inbound
    true_anomaly = true_anomaly .* -1;
end
eccentric_anomaly = 2 * atan(sqrt((1 - eccen) / (1 + eccen)) * tan(true_anomaly ./ 2));
M_1 = eccentric_anomaly - (eccen .* sin(eccentric_anomaly));
time = sqrt((a .^ 3) ./ mu) .* M_1;
end
