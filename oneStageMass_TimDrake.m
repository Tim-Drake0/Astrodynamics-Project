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

MR = exp(deltaV ./ (9.81 .* sum(Isp)));
propellant_mass_frac = (MR - 1) ./ MR;
payload_ratio = 1 - (propellant_mass_frac ./ (1 - lambda));
m_initial = m_payload ./ payload_ratio;
m_propellant = propellant_mass_frac .* m_initial;
m_structural = propellant_mass_frac .* m_initial .* (lambda ./ (1 - lambda));
end