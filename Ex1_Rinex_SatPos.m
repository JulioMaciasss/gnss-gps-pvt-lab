%% ETSIT MIT GPS Lab 2025
%
% Script 1: lectura de archivos RINEX y calculo posicion satelites 
%
% @Omar Garcia Crespillo
%-------------------------------------------------------------------------------

%clear all
%close all
clearvars;

addpath(genpath('../'));
addpath('Data');
addpath(genpath('Functions'));
addpath('Rinex');

filename_nav = 'Data/RinexParnav.16n';
filename_obs = 'Data/RinexParobs.16o';

%% CONSTANTS
c = 0.299792458e9;  %% Speed of light
omega_edot = 7.2921151467e-5;  %% Earth rate

%% Read navigation RINEX file(s)
if(~exist('Eph')),
[Eph, iono] = RINEX_get_nav(filename_nav, []);

%% Read observation RINEX file(s)
[obs, ~, tow, week_R, date_R, pos_R, interval, antoff_R, antmod_R] = ...
    load_RINEX_obs(filename_obs, []);
end

%keyboard; %% Analizar la estructura de las variables devueltas con las Ephemerides 'Eph', con las medidas 'obs' y 'tow'.


Nepochs = 1;
Nepochs = size(date_R,1);  %% Uncomment

%% Main Loop over time
for ii=1:Nepochs,
    
    % Find visible satellites
    svIdx = find(~isnan(obs.C1(:,ii)));
    
    % Extract the pseudorange measurements
    pseudoranges = obs.C1(svIdx,ii);
    
    %% Computation of satellite position. Loop over each satellite in view with index jj.
    for jj=1:length(svIdx),
        
    % Compute the time of transmission
    %% COMPLETAR AQUI (Solution)
    
    timeTx(jj) = tow(ii) - pseudoranges(jj)/c;
    
    %% Choose the closes Ephemeris
    % % En la posicion 1 de Eph, tenemos el Numero PRN del satellite, en la
    % posicion 18 tenemos el toe (time of Ephemeris)
    idMatch = find(Eph(1,:)==svIdx(jj) & Eph(18,:)>=timeTx(jj));  
    [~, k]=min(abs(Eph(18,idMatch)-timeTx(jj)));
    ephSelect = idMatch(k);
    
    %% Compute Satellite position
    % Necesita el tiempo en que queremos calcular la posicion del satelite
    % (time of transmision), todos los parametros de las Ephemerides (la
    % Eph elegida) y el Numero PRN del satelite del que queremos la
    % posicion.
    [satpos(jj,:), satvel(jj,:), ~] = satellite_orbits_clock(timeTx(jj), Eph(:,ephSelect), svIdx(jj)); 
 
    
    %% Convertir de ECEF a LLH la posicion de los satellites (Cuidado, satpos esta en: latitude[rad], longitud[rad], altitud[m])
    %% COMPLETAR AQUI
    satposLLH(jj,1:3) = cart2geod(satpos(jj,:)) ; 
    satposLLHdeg(jj,:) = [180/pi * (satposLLH(jj,1:2)),satposLLH(jj,3)] ; % Y aqui en [deg,deg,m]
   
    end
    satpos_all(ii,:,:) = satpos;

    %% Calculo velocidad de satelites
    %% COMPLETAR AQUI
    satVelTotal = sqrt(satvel(:,1).^2+satvel(:,2).^2+satvel(:,3).^2); 

    fprintf('Epoch %d ================= \n',ii);
    disp('Satellite positions -------------');
    disp('Lat long in [deg]');
    disp(satposLLHdeg(:,1:2));
    disp('Altitude [m]');
    disp(satposLLHdeg(:,3));

    disp('Pseudoranges [m] --------------');
    disp(pseudoranges);
    disp('Doppler [Hz] -------------');
    disp(obs.D1(svIdx,ii));

    disp('Velocidad Satelites [m/s]');
    disp(satVelTotal);
end

%% Uncomment when all epochs are considered
plotSatelliteViewer;