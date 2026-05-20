%% ETSIT MIT GPS Lab 2025
%
% Script 2: Script principal de lectura RINEX, calculo position y exportar
% a Google Earth
%
% @Omar Garcia Crespillo
%-------------------------------------------------------------------------------
%clear all
%close all
clearvars;

addpath(genpath('../'));

filename_nav = 'Data/RinexParnav.16n';
filename_obs = 'Data/RinexParobs.16o';

%% CONSTANTS
c = 0.299792458e9;  %% Speed of light
omega_edot = 7.2921151467e-5;  %% Earth rate


%% Read navigation RINEX file(s)
if(~exist('Eph')),
    [Eph, iono] = RINEX_get_nav(filename_nav, []);
    
    %% Read observation RINEX file(s)
    [obs, time_R, tow, week_R, date_R, pos_R, interval, antoff_R, antmod_R] = ...
        load_RINEX_obs(filename_obs, []);
end

%% Total number of Epochs extracted from data
Nepochs = size(obs.C1,2);


%% Reserve some memory for final variables
solution = NaN(4,Nepochs);
solutionLLH = NaN(3,Nepochs);
solutionLLHdeg = NaN(3,Nepochs);

disp('Starting Loop Epochs .......');
%% Main Loop over time
for ii=1:Nepochs,
    
    % Find visible satellites
    svIdx = find(~isnan(obs.C1(:,ii)));
    
    % Extract the pseudorange measurements
    pseudoranges = obs.C1(svIdx,ii);
    
    satpos = NaN(length(svIdx),3);
    satvel = NaN(length(svIdx),3);
    dtS = NaN(length(svIdx),1);
    %% Computation of satellite position
    for jj=1:length(svIdx),
        
        % First guess of Time of transmission
        %% COMPLETAR COMO EN Ex1_Rinex_SatPos
        timeTx(jj) = tow(ii) - pseudoranges(jj)/c ; 
        
        % Choose the closes Ephemeris
        idMatch = find(Eph(1,:)==svIdx(jj) & Eph(18,:)>=timeTx(jj));  %Coincide
        [~, k]=min(abs(Eph(18,idMatch)-timeTx(jj)));
        ephSelect = idMatch(k);
        
        % Compute Satellite position
        if ~isempty(ephSelect), %% to prevent if no ephemeris found
            [satpos(jj,:), satvel(jj,:), dtS(jj)] = satellite_orbits_clock(timeTx(jj), Eph(:,ephSelect), svIdx(jj));
        end
    end
    % If we are missing Ephemeris for a visible satellite, we eliminate the satellite
    nanIdx = find(isnan(satpos(:,1)));
    pseudoranges(nanIdx) = [];
    svIdx(nanIdx) = [];
    timeTx(nanIdx) = [];
    satpos(nanIdx,:) = [];
    satvel(nanIdx,:) = [];
    dtS(nanIdx,:) = [];


    % Approximate initial solution and initial clock guess (e.g. center of
    % Earth)
    apcoords = [0 0 0 0]; 
    
    if length(svIdx) > 3,  % we need 4 sat mininum to compute the solution% Least-Squares Nonlinear solver
        [X0,Y0,Z0,dtu,PDOP(ii),TDOP(ii),GDOP(ii)] = solveLS(svIdx,pseudoranges,satpos,dtS, apcoords,Eph, tow(ii), iono);

        solution(:,ii) = [X0, Y0, Z0, dtu];
        
        %% COMPLETAR COMO EN Ex1_Rinex_SatPos
        solutionLLH(:,ii) = cart2geod(solution(1:3,ii));
        solutionLLHdeg(:,ii) = [180/pi * (solutionLLH(1:2,ii)) ; solutionLLH(3,ii)] ; 
        
        fprintf('Current position in LLH [deg,deg,m]: [%10f,%10f,%6f] \n',solutionLLHdeg(1,ii),solutionLLHdeg(2,ii),solutionLLHdeg(3,ii));
        fprintf('PDOP: %d, TDOP: %d, GDOP: %d \n',PDOP(ii),TDOP(ii),GDOP(ii));
    end
end
disp('Finished all Epochs');

% Remove Epochs without position (for example because <4 satellites in view)
nanIdx = find(isnan(solutionLLHdeg(1,:)));
solutionLLHdeg(:,nanIdx) = [];

%% Export to KML (Google Earth)
kmlFileName = 'output.kml';  % Cambia este nombre cada vez que avances en algo el codigo
generatorKML(solutionLLHdeg',kmlFileName,1);

