% Mission Profile
%t1 = 20; % [days] Departure time from Earth 
%t2 = 282; % [days] Arrival time at the NEO 
%t3 = 500; % [days] Departure time from the NEO 
%t4 = 915; % [days]  Arrival time at the ISS 
%DV1 = 11.616; % liftoff from Earth 
%DV2 = 3.0735; % landing on the NEO,
%DV3 = 1.9541; % liftoff from the NEO
%DV4 = 3.3284; % rendezvous with the ISS.
%m_liftOff = 3,232,295.25; % liftoff mass of the rocket
%num_stages = 4; % number of stages of the rocket
%cost = $ 11.3396 billion
% DV total = 19.972 km/s

format long g
% Earth Constants
mu_Earth = 3.986e5; % km^3/s^2
r_Earth = 6378; % km
a_Earth = 1.496e+8; % km
eccen_Earth = 0.0167;
SOI_Earth = 0.929e6/a_Earth; %km

% Apollo Constants
mu_Apollo = 4.5e-4; % km^3/s^2
r_Apollo = 10; % km
a_Apollo = 1.5109e+8; % km
eccen_Apollo = 0.0200;

mu_Sun = 1.327e11; %km^3/s^2

minDV = inf;
graphOrbits = true;

for t1 = 20
    for t2 = 282
        for t3 = 500
            for t4 = 915
                if t4 > t3 && t3 > t2 && t2 > t1
                    [arTimes, arDVs, Apollodata, Earthdata, Transferdata] = FinalProject_TimDrake(t1*86400, t2*86400, t3*86400, t4*86400);
                    if sum(arDVs) < minDV && sum(arDVs) ~= 0
                        minDV = sum(arDVs);
                        minTime = arTimes(4)-arTimes(1);
                        aTimes = arTimes;
                        aDVs = arDVs;
                        Apollo_data = Apollodata;
                        Earth_data = Earthdata;
                        Transfer_data = Transferdata;
                    end
                end
            end
        end
    end
end

% Rocket Design
N_stage = 10;
minCost = inf;
m_PL = 5000 + 10 * minTime; % [kg]
Isp = 450 .* ones(1, N_stage); % [seconds]
lam = 0.15 .* ones(1, N_stage-5);
lambda = cat(2, lam, [0.15 0.12 0.1 0.08 0.05]); % structural mass fraction 
for numStages = 2:N_stage    
    [maxDV, maxMR] = getMaxDV(Isp(N_stage+1-numStages:N_stage), lambda(N_stage+1-numStages:N_stage)); % max DV per stage no payload
    [m_i, ~, ~, stage_DV] = multistage_sameMR_TimDrake(sum(aDVs)*1000, m_PL, Isp(N_stage+1-numStages:N_stage), lambda(N_stage+1-numStages:N_stage), maxDV, maxMR);
    if all(m_i > 0)
        cost = 5 * (minTime) + (m_i(1) / 500) + (100 * numStages); % in millions of dollars
        if cost < minCost
            minCost = cost;
            m_initial = m_i;
            mPL = m_PL;
            num_stages = numStages;
            stageDV = stage_DV;
        end
    end
end

unit = ' million\n';
if minCost > 1000
    minCost = minCost/1000;
    unit = ' billion\n';
end

% Print results
fprintf("Cost: $ " + string(minCost) + unit)
fprintf("Time: " + string(minTime) + ' Days \n')
fprintf("Delta V: " + string(minDV) + ' km/s \n')
fprintf("Liftoff mass: " + string(m_initial(1)) + ' kg \n')
stagePrint = "Number of stages: " + string(num_stages) + "\n";
for i = 1:num_stages
        stagePrint = stagePrint + "    S" + string(i) + ": " + string(m_initial(i)) + " kg (" + string(stageDV(i)/1000) + " km/s)\n";
end
fprintf(stagePrint + "   mPL: " + string(mPL) + " kg\n")
fprintf("Input times: \n    t1: " + string(aTimes(1)) + ' days \n' + ... 
                       "    t2: " + string(aTimes(2)) + ' days \n' + ...
                       "    t3: " + string(aTimes(3)) + ' days \n' + ...
                       "    t4: " + string(aTimes(4)) + ' days \n') 
fprintf("Delta V's: \n    DV1: " + string(aDVs(1)) + ' km/s \n' + ...
                     "    DV2: " + string(aDVs(2)) + ' km/s \n' +...
                     "    DV3: " + string(aDVs(3)) + ' km/s \n' +...
                     "    DV4: " + string(aDVs(4)) + ' km/s \n') 
                 
if graphOrbits == true % Graph results
    hold on
    plot(0, 0, 'r*') % focus
    graphOrbit(a_Earth, eccen_Earth, 0, 2*pi, 1, false, 0)
    graphOrbit(a_Apollo, eccen_Apollo, 0, 2*pi, 1, false, 0)
    
    aT1 = Transfer_data(1,1);
    eT1 = Transfer_data(1,2);
    aT2 = Transfer_data(2,1);
    eT2 = Transfer_data(2,2);
    
    E1=Earth_data(1,2);
    E2=Earth_data(2,2);
    E3=Earth_data(3,2);
    E4=Earth_data(4,2);
    ET1=Transfer_data(1,3);
    ET4=Transfer_data(2,3);
    
    A1=Apollo_data(1,2);
    A2=Apollo_data(2,2);
    A3=Apollo_data(3,2);
    A4=Apollo_data(4,2);
    AT2=Transfer_data(1,4);
    AT3=Transfer_data(2,4);
    
    compareE = [E1 E4; ET1 ET4];
    compareA = [A2 A3; AT2 AT3];
    
    graphOrbit(aT1, eT1, 0, AT2-ET1, 3, E1, a_Earth-aT1)
    graphOrbit(aT2, eT2, 0, AT3-ET4, 3, A3, a_Apollo-aT2)
    
    circle(cos(Earth_data(4,2)),sin(Earth_data(4,2)),SOI_Earth)
    
    graphPlanet(a_Apollo, eccen_Apollo, A1, 'A1')
    graphPlanet(a_Earth, eccen_Earth, E1, 'E1')
    
    graphPlanet(a_Apollo, eccen_Apollo, A2, 'A2')
    graphPlanet(a_Earth, eccen_Earth, E2, 'E2')
    
    graphPlanet(a_Apollo, eccen_Apollo, A3, 'A3')
    graphPlanet(a_Earth, eccen_Earth, E3, 'E3')
    
    graphPlanet(a_Apollo, eccen_Apollo, A4, 'A4')
    graphPlanet(a_Earth, eccen_Earth, E4, 'E4')
    
    plot([0 cos(Earth_data(1,2))], [0 sin(Earth_data(1,2))],':')
    plot([0 cos(Apollo_data(2,2))], [0 sin(Apollo_data(2,2))],':')
    plot([0 cos(Earth_data(4,2))], [0 sin(Earth_data(4,2))],':')
    plot([0 cos(Apollo_data(3,2))], [0 sin(Apollo_data(3,2))],':')
    
    axis([-1.5 1.5 -1.5 1.5])
    axis square
    grid on
    legend('Sun', 'Earth Orbit', 'Apollo Orbit', 'Transfer t1-->t2', 'Transfer t3-->t4')
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

function [x, y] = rotate(x, y, theta)
    v = [x;y];
    for i=1:numel(x)
        R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
        vR = R*v(:,i);
        x(:,i) = vR(1);
        y(:,i) = vR(2);
    end
end

function [maxDV, maxMR] = getMaxDV(Isp, lambda)
    maxMR = 1 ./ lambda;
    maxDV = -9.81 .* Isp .* log(lambda);
end

function graphOrbit(semimajor, eccen, f1, f2, width, rot, shift)
   
    shift = shift/1.496e+8;
    f = linspace(f1, f2);
    [c, ~, xp, yp] = getPos(semimajor, eccen, f);
    
    if rot ~= 0
        [xp, yp] = rotate(xp+c+shift, yp, rot);
        plot(xp, yp,'LineWidth',width); % Orbit
        return
    end
    
    
    plot(xp+c+shift, yp,'LineWidth',width); % Orbit
end

function [c, radius, xp, yp] = getPos(semimajor, eccen, f)
    c = (semimajor * eccen)/1.496e+8;
    radius = (semimajor .* (1 - eccen .^2)) ./ (1 + (eccen .* cos(f)));
    xp = (radius .* cos(f))/1.496e+8;
    yp = (radius .* sin(f))/1.496e+8;
end

function graphPlanet(semimajor, eccen, f, name)
    [c, ~, xp, yp] = getPos(semimajor, eccen, f);
    plot(xp+c, yp,'o') % planet
    text(xp+c, yp, string(name))
end

function circle(x,y,r)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit);
end
