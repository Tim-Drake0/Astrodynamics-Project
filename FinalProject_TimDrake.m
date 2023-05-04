function [times, DVs, Apollo_data, Earth_data, Transfer_data] = FinalProject_TimDrake(t1, t2, t3, t4)
    % Input times are in seconds
    times = 0; % array of each time input
    DVs = [0 0 0 0]; % array of delta V's for N stages 
    Apollo_data = []; % 3x2 matrix of apollo data. for graphing, not really needed
    Earth_data = []; % 3x2 matrix of Earth data. for graphing, not really needed
    Transfer_data = []; % 3x2 matrix of transfer data. for graphing, not really needed
    
    % Earth Constants
    mu_Earth = 3.986e5; % [km^3/s^2]
    r_Earth = 6378; % [km]
    a_Earth = 1.496e+8; % [km]
    eccen_Earth = 0.0167;

    % Apollo Constansts
    mu_Apollo = 4.5e-4; % [km^3/s^2]
    r_Apollo = 10; % [km]
    a_Apollo = 1.5109e+8; % [km]
    eccen_Apollo = 0.0200;
    f_Apollo = 1.5708; % [rad]
    mu_Sun = 1.327e11; % [km^3/s^2]

    % Mission:
    % Time offset of Apollo from Earth at t1
    [Apollo_offset, ~, ~] = time_from_anomaly_TimDrake(mu_Sun, a_Apollo, eccen_Apollo, f_Apollo, true);
    
    % Calculate Apollo and Earth positions at t1
    [r_Apollo_t1, f_Apollo_t1, ~] = position_from_time_TimDrake(mu_Sun, a_Apollo, eccen_Apollo, t1+Apollo_offset);
    [r_Earth_t1, f_Earth_t1, ~] = position_from_time_TimDrake(mu_Sun, a_Earth, eccen_Earth, t1);

    % Apollo and Earth positions at t2
    [r_Apollo_t2, f_Apollo_t2, ~] = position_from_time_TimDrake(mu_Sun, a_Apollo, eccen_Apollo, t2+Apollo_offset);
    [r_Earth_t2, f_Earth_t2, ~] = position_from_time_TimDrake(mu_Sun, a_Earth, eccen_Earth, t2);
    
    % Velocity of Earth at t1
    V_Earth_t1 = sqrt(((2 * mu_Sun) ./ r_Earth_t1) - ((mu_Sun / a_Earth))); % [km/s]
    
    % Transfer orbit from Earth(t1) to Apollo(t2)
    [aT1, eT1, f_Transfer1_t2, f_Transfer2_t2, isConverging] = lambert_TimDrake(mu_Sun, r_Earth_t1, r_Apollo_t2, f_Earth_t1, f_Apollo_t2, t2-t1);
    if isConverging == false %if lambert isnt converging then go to next times
        return
    end
    Vinf1 = calcVinf(mu_Sun, r_Earth_t1, aT1, a_Earth, eT1, eccen_Earth, V_Earth_t1, true, true); % [km/s]
    DV1 = sqrt(((2 * mu_Earth) ./ r_Earth) + Vinf1.^2); % [km/s] DV req to leave Earth SOI

    % Velocity of Apollo at t2
    V_Apollo_t2 = sqrt(((2 * mu_Sun) ./ r_Apollo_t2) - ((mu_Sun / a_Apollo))); % [km/s]

    Vinf2 = calcVinf(mu_Sun, r_Apollo_t2, aT1, a_Apollo, eT1, eccen_Apollo, V_Apollo_t2, false, false);
    DV2 = sqrt(((2 * mu_Apollo) ./ r_Apollo) + Vinf2.^2);  % DV req to land on Apollo

    % Calculate Apollo and Earth positions at t3
    [r_Apollo_t3, f_Apollo_t3, ~] = position_from_time_TimDrake(mu_Sun, a_Apollo, eccen_Apollo, t3+Apollo_offset);
    [r_Earth_t3, f_Earth_t3, ~] = position_from_time_TimDrake(mu_Sun, a_Earth, eccen_Earth, t3);
    
    % Calculate Apollo and Earth positions at t4
    [r_Apollo_t4, f_Apollo_t4, ~] = position_from_time_TimDrake(mu_Sun, a_Apollo, eccen_Apollo, t4+Apollo_offset);
    [r_Earth_t4, f_Earth_t4, ~] = position_from_time_TimDrake(mu_Sun, a_Earth, eccen_Earth, t4);
    
    % Velocity of Apollo at t3
    V_Apollo_t3 = sqrt(((2 * mu_Sun) ./ r_Apollo_t3) - ((mu_Sun / a_Apollo))); % [km/s]

    % Transfer orbit from Apollo(t3) to Earth(t4)
    [aT2, eT2, f_Transfer1_t4, f_Transfer2_t4, isConverging] = lambert_TimDrake(mu_Sun, r_Apollo_t3, r_Earth_t4, f_Apollo_t3, f_Earth_t4, t4-t3);
    if isConverging == false %if lambert isnt converging then go to next times
        return
    end
    Vinf3 = calcVinf(mu_Sun, r_Apollo_t3, aT2, a_Apollo, eT2, eccen_Apollo, V_Apollo_t3, true, true); % [km/s]
    DV3 = sqrt(((2 * mu_Apollo)/ r_Apollo) + Vinf3.^2) - 0; % [km/s] DV req to leave Apollo SOI

    % Velocity of Earth at t4
    V_Earth_t4 = sqrt(((2 * mu_Sun) ./ r_Earth_t4) - ((mu_Sun / a_Earth)));% [km/s]
    
    % Rendezvous with the ISS
    r_ISS = 6778;% [km]
    V_ISS = sqrt(((2 * mu_Earth) ./ r_ISS) - ((mu_Earth / r_ISS)));% [km/s]
    
    Vinf4 = calcVinf(mu_Sun, r_Earth_t4, aT2, a_Earth, eT2, eccen_Earth, V_Earth_t4, true, true);% [km/s]
    DV4 = abs(V_ISS - sqrt(((2 * mu_Earth)/ r_ISS) + Vinf4.^2)); % [km/s] DV req to enter ISS orbit

    times = [t1 t2 t3 t4]./86400; % [days]
    DVs = [DV1 DV2 DV3 DV4]; % [km/s] 

    % For graphing orbits:
    Apollo_data = [r_Apollo_t1 (f_Apollo_t1); 
        r_Apollo_t2 (f_Apollo_t2); 
        r_Apollo_t3  (f_Apollo_t3); 
        r_Apollo_t4  (f_Apollo_t4)];
    Earth_data = [r_Earth_t1 (f_Earth_t1); 
        r_Earth_t2 (f_Earth_t2); 
        r_Earth_t3 (f_Earth_t3);  
        r_Earth_t4 (f_Earth_t4)];
    Transfer_data = [aT1 eT1 (f_Transfer1_t2) (f_Transfer2_t2); 
                    aT2 eT2 (f_Transfer1_t4) (f_Transfer2_t4)];

end

function Vinf = calcVinf(mu, radius_obj, a_transfer, a_obj, eccen_transfer, eccen_obj, V_obj, obj_inbound, transfer_inbound)
    V_transfer = sqrt(((2 * mu) ./ radius_obj) - ((mu ./ a_transfer)));
    h_transfer = sqrt(mu .* a_transfer .* (1 - eccen_transfer.^2));
    phi_transfer = acos(h_transfer ./ (radius_obj .* V_transfer));
    if transfer_inbound
        phi_transfer = -phi_transfer;
    end
    h_obj = sqrt(mu .* a_obj .* (1 - eccen_obj.^2));
    phi_obj = real(acos(h_obj ./ (radius_obj .* V_obj)));
    if obj_inbound
        phi_obj = -phi_obj;
    end
    Vinf = sqrt((V_transfer^2) + (V_obj^2) - (2 * V_transfer * V_obj * cos(phi_transfer-phi_obj)));
end

function [semimajor, eccen, f1, f2, isConverging] = lambert_TimDrake(mu, r1, r2, f1, f2, time)
% Tim Drake Astrodynamics Project 3
% Computes semimajor, eccentricity, and 2 true anomaly's of an orbit.
% INPUTS:
%   mu = [km^3/s^2] gravitational constant of the Central body
%   r1 = [km] the radial distance from the central body at the initial time
%   r2 = [km] the radial distance from the central body at the final time
%   theta = [rad] the angle between positions 1 and 2 in radians; if theta > π, then the desired trajectory goes the "long" way around the central body
%   time = [seconds] the travel time from position 1 to position 2. Can be
%          scalar or vector.
% OUTPUTS:
%   semimajor = [km] the semimajor axis of the elliptical orbit that satisfies the Lambert problem
%   eccen = the orbit eccentricity of the elliptical orbit that satisfies the Lambert problem
%   true_anomaly_1 = [rads] the true anomaly of the satellite in radians at position 1 (where inbound values are negative)
%   true_anomaly_2 = [rads] the true anomaly of the satellite in radians at position 2 (where inbound values are negative)
% Angle between Apollo(t2) and Earth(t1)
theta = f2-f1;
if theta < 0
    theta = 2*pi + theta;
end
count = 0;
isConverging = true;
c = sqrt(r1.^2 + r2.^2 - 2 .* r1 .* r2 .* cos(theta));
s = .5 .* (r1 + r2 + c);
a_m = .5 .* s;
beta_m = 2 * asin(sqrt((s - c) / (2 * a_m)));
if theta > pi
    beta_m = -beta_m;
end
t_m = sqrt(a_m.^3/mu) .* (pi - beta_m + sin(beta_m));
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
        beta = -beta;
    end
    time_calculated = sqrt(a_guess.^3/mu) .* (alpha - beta - (sin(alpha) - sin(beta)));
    eVal = time_calculated - time;
    if (eVal > 0)
        % replace the high value
        high = a_guess;
    else
        low = a_guess;
    end
    count = count + 1;
    
    if count > 1000
    %    %error("Lambert is not converging")
        isConverging = false;
        semimajor = 0;
        eccen = 0;
        f1 = 0;
        f2 = 0;
        return
    end
end
semimajor = a_guess;
p = 4 .* ((semimajor .* (s - r1) .* (s - r2)) ./ c.^2) .* (sin((alpha + beta) ./ 2).^2);
eccen = sqrt(1 - (p ./ semimajor));

f1 = acos((1/eccen)*(semimajor*(1-eccen^2)/r1-1));
timeTest = time_from_anomaly_Check(f1, f1+theta, semimajor, eccen, mu);
if abs(timeTest - time) > 1e-3
    % Flip and hope it works
    f1 = -f1;
    timeTest = time_from_anomaly_Check(f1, f1+theta, semimajor, eccen, mu);
    if abs(timeTest - time) > 1e-3
        fprintf('Error on getting the true anomaly right: f1=%f and time=%f but time should=%f\n', ...
            f1, timeTest, time);
    end
end
f2 = f1 + theta;
end

function elapsedTime = time_from_anomaly_Check(f1, f2, a, e, mu)
      sqe = sqrt((1-e)/(1+e));
      u1 = 2*atan(sqe*tan(f1/2));
      u2 = 2*atan(sqe*tan(f2/2));
      M1 = u1 - e*sin(u1);
      M2 = u2 - e*sin(u2);
      while (M2-M1 < 0)
            M2 = M2 + 2*pi;
      end
      elapsedTime = sqrt(a^3/mu)*(M2-M1);
end

function [m_initial, m_propellant, m_structural, deltaV] = multistage_sameMR_TimDrake(deltaVtotal, m_payload, Isp, lambda, maxDV, maxMR)
% Tim Drake Astrodynamics Project 4b
% Calculates the masses of all the elements of a multistage rocket, assuming that each stage has the same mass ratio.
% INPUTS: deltaVtotal = [m/s] required velocity change for the mission (meters/sec)
        %   m_payload = [kg] payload mass on the top of the rocket
        %         Isp = [seconds] specific impulse of each stage's engine (seconds) [vector]
        %      lambda = structural mass fraction of each stage [vector]
% OUTPUTS:  m_initial = [kg] initial mass of the each stage [vector]
%        m_propellant = [kg] propellant mass of each stage [vector]
%        m_structural = [kg] structural mass of each stage [vector]
%              deltaV = [m/s] velocity change caused by each stage (meters/sec) [vector]
[m_initial, m_propellant, m_structural] = deal(zeros(1, length(Isp)));
MR = exp(deltaVtotal / (9.81 * sum(Isp)));
deltaV = 9.81 .* Isp .* log(MR);
for i = length(Isp):-1:1
    if MR > maxMR(i)
        return
    end
    if deltaV(i) > maxDV(i)
        return
    end
    [m_initial(i), m_propellant(i), m_structural(i)] = oneStageMass_TimDrake(deltaV(i), m_payload, Isp(i), lambda(i));
    m_payload = m_initial(i);
end
end

function [m_initial, m_propellant, m_structural] = oneStageMass_TimDrake(deltaV, m_payload, Isp, lambda)
% Tim Drake Astrodynamics Project 4a
% Calculates the intitial, propellant, and structural masses of a rocket stage of the given values.
% INPUTS: 
%   deltaV = [m/s] velocity change caused by this stage (meters/sec)
%   m_payload = [kg] payload mass carried by this stage
%   Isp = [seconds] specific impulse of this stage's engine (seconds)
%   lambda = structural mass fraction of this stage
% OUTPUTS:
%   m_initial = [kg] initial mass of the entire stage
%   m_propellant = [kg] propellant mass of this stage
%   m_structural = [kg] structural mass of this stage

MR = exp((deltaV) ./ (9.81 .* Isp));
propellant_mass_frac = (MR - 1) ./ MR;
payload_ratio = 1 - (propellant_mass_frac ./ (1 - lambda));
m_initial = m_payload ./ payload_ratio;
m_propellant = propellant_mass_frac .* m_initial;
m_structural = propellant_mass_frac .* m_initial .* (lambda ./ (1 - lambda));
end

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

function [time, true_anomaly, eccentric_anomaly] = time_from_anomaly_TimDrake(mu, a, eccen, true_anomaly, is_outbound)
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
%true_anomaly = real(acos((1 ./ eccen) * (((a .* (1 - (eccen .^ 2))) ./ radius) - 1)));
if ~is_outbound % f is negative for inbound
    true_anomaly = -true_anomaly;
end
eccentric_anomaly = 2 * atan(sqrt((1 - eccen) / (1 + eccen)) * tan(true_anomaly ./ 2));
M_1 = eccentric_anomaly - (eccen .* sin(eccentric_anomaly));
time = sqrt((a .^ 3) ./ mu) .* M_1;
end


















