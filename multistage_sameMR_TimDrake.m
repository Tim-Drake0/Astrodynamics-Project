function [m_initial, m_propellant, m_structural, deltaV] = multistage_sameMR_TimDrake(deltaVtotal, m_payload, Isp, lambda)
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

MR = exp(deltaVtotal / (9.81 * sum(Isp)));
deltaV = 9.81 .* Isp .* log(MR);
[m_initial, m_propellant, m_structural] = deal([0 0 0]);
for i = length(Isp):-1:1
    [m_initial(i), m_propellant(i), m_structural(i)] = oneStageMass_TimDrake(deltaV(i), m_payload, Isp(i), lambda(i));
    m_payload = m_initial(i);
end
end
