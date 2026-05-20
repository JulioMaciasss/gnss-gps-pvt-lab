%% LS solver
% ETSIT MIT GPS Lab 2025
%
% @Omar Garcia Crespillo, 2025
%------------------------------------------------------------------------------------
function [X0,Y0,Z0,dtu,PDOP,TDOP,GDOP]=solveLS(svIdx, pseudoranges, satpos, dtS, apcoords, Eph, tow, iono)


%% Constants
c = 0.299792458e9; % speed of light in m/s
omega_edot = 7.2921151467e-5; % Earth rotation constant


%% For the iterative non linear solution
Nit_max= 15; % max number of iterations
Nit = 1; % iterations counter

% Variable Initialization
new_userPos = apcoords;
userPos(1:3) = satpos(1,:);
Npseudo = length(pseudoranges);
dtS = zeros(length(svIdx),1);

%% Main linearization least-squares LOOP ======================================
while (norm(userPos(1:3)-new_userPos(1:3))>1e-9) && Nit <= Nit_max
    userPos =  new_userPos;
    
    %% 1. Iterar por cada satelite y calcular su posicion --------------        
    for jj=1:length(svIdx)

        %% TIP: Util para seleccionar la Ephemeris correcta
        t_tx(jj) = tow - pseudoranges(jj)/c - dtS(jj); 
        t_rx(jj) = tow - userPos(4)/c;
        travelTime(jj) = t_rx(jj)-t_tx(jj);
       
        % satposLLH(jj,1:3) = cart2geod(satpos(jj,:)) ; 
        % satposLLHdeg(jj,:) = [180/pi * (satposLLH(jj,1:2)),satposLLH(jj,3)] ; % Y aqui en [deg,deg,m]
        
        idMatch = find(Eph(1,:)==svIdx(jj) & Eph(18,:)>=t_tx(jj)); 
        [~, k]=min(abs(Eph(18,idMatch)-t_tx(jj)));
        ephSelect = idMatch(k);

        [satpos(jj,:), satvel(jj,:), dtS(jj)] = satellite_orbits_clock(t_tx(jj), Eph(:,ephSelect), svIdx(jj)); 

        omega_matr = [cos(omega_edot * travelTime(jj)) , sin(omega_edot * travelTime(jj)), 0; 
            -sin(omega_edot * travelTime(jj)), cos(omega_edot * travelTime(jj)),0;
            0 ,0 , 1 ];
        satpos(jj,:) = omega_matr * satpos(jj,:)'; 

    %% PARTE 2 (ignorar para la parte 1) ====================
        
         % A. Sagnac Correction   
         % Escribir aqui para cada satelite 
         % B. Errores atmosfericos (descomentar en parte 2)
        [Az(jj), El(jj), D(jj)] = topocent(userPos(1:3),satpos(jj,:));
        llh = cart2geod(userPos(1:3));
        dIono(jj) = iono_error_correction(llh(1)*180/pi, llh(2)*180/pi, Az(jj), El(jj), t_tx(jj), iono, []);
        dTropo(jj) = tropo_error_correction(El(jj),llh(3));
          
        R_k = norm(satpos(jj,:) - userPos(1:3)) 

        H(jj,:) =  [-[satpos(jj,1)  - userPos(1), satpos(jj,2)  - userPos(2), satpos(jj,3)  - userPos(3)]./R_k, 1];

        z(jj,1) = pseudoranges(jj,:) - R_k +c*dtS(jj);% -dIono(jj) -dTropo(jj);

    end   

    x_ls = (H' * H)^(-1) * H' * z
    new_userPos(1) = userPos(1) + x_ls(1)
    new_userPos(2) = userPos(2) + x_ls(2)
    new_userPos(3) = userPos(3) + x_ls(3)
    new_userPos(4) = userPos(4) + x_ls(4)

       
      
        %% =========================================================
        
    %% ------------------------------------------------------------
    
    %% 2. Solucion al problema de optimizacion LS --------------------
    
    %% --------------------------------------------------------------
    fprintf('Iteracion LS Numero %d \n',Nit);
    Nit = Nit+1;
end

% Asignar variables de salida
X0 = new_userPos(1);
Y0 = new_userPos(2);
Z0 = new_userPos(3);
dtu = new_userPos(4);

%% PARTE 2 (Ignorar para parte 1) ====
PDOP = 0;
GDOP = 0;
TDOP = 0;
%% ===================================
