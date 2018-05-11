% Übung 1, Aufgabe 3b, Numerik 1 
% Theresa Kiszler
% Abgabe am 23.04.2018

function fixpunkt = run_3_1()
%RUN_3_1() berechnet einen Fixpunkt durch Fixpunktiteration
  % Initialisiere Anfangswerte
  x_1old = 0; x_2old = 0;
  
  % Teste unterschiedliche Startwerte 
  for x_1 = 0.1:0.5:3;
  
    fprintf('\n -------------------######----------------------- \n');
    fprintf('Suche Fixpunkt mit Startwert x = (%d, 1) \n', x_1);
    x_2 = 1;
    x_k       = [x_1, x_2];
    x_kminus1 = [x_1old, x_2old];
    % Initialisiere Genauigkeit
    converge = 1;
    
    % While loop laeuft solange bis Genauigkeit von 10^-8 garantiert ist.
    while (converge > 10^(-8))
       
      x_kplus1 = phi(x_k);
      theta_k = (norm(x_kplus1 - x_k) / norm(x_k - x_kminus1));      
      converge = theta_k/(1-theta_k) * norm(x_kplus1 - x_k);
      
      % Unterbreche Iteration falls kein Erfolg moeglich.
      if (theta_k >= 1) 
        fail = true;
        disp('Kein Erfolg moeglich!');
        x_k = [NA, NA];
        break
      endif
      
      % Variablen einen Schritt weiter setzen.
      x_minus1 = x_k;
      x_k = x_kplus1;
    endwhile
    
    fixpunkt = x_k;
    fprintf(' Gefundener Fixpunkt: (%d, %d) \n', x_k(1), x_k(2));
  
  endfor

endfunction

function phi = phi(x_k)
%PHI(x_k) benutzt F um den naechsten Schritt zu ermitteln
  phi = [(1/3 * x_k(2)^2 + 1/8), (1/4*x_k(1)^2 - 1/6)];
endfunction
