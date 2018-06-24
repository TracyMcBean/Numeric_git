% Numerik 1 SS18
% Uebung 10, Aufgabe 1
% David und Tracy

function [] = run_10_1()
% Test multilevel Quadratur Funktion
  a=0; b=1; 
  n = (0:10:500);
  Nmax = 1000;
  tol = 10^(-10);
  % Test 1 y = x
%  disp('--------------------------------------')
%  disp('Test 1: y = x');
%  y1 = multilevel_quad(@testfunktion1, a, b, n, Nmax, tol)
%  
  % Test 2 y = sqrt(x)
  disp('--------------------------------------')
  disp('Test 2: sqrt(x)');
  for i = 1:length(n)
    y_sqrt(i) = multilevel_quad(@sqrt, a, b, n(i), Nmax, tol);
    error_sqrt(i) = norm(y_sqrt(i) - 2/3);
  end
  
  
  % Test 3 y = sin(x)
  a=0; b=pi; 
  disp('--------------------------------------')
  disp('Test 3: sin(x)');
  for i = 1:length(n)
    y_sin(i) = multilevel_quad(@sin, a, b, n(i), Nmax, tol);
    error_sin(i) = norm(y_sin(i) - 2);
  end 
  
  % Test 4 arctan
  a=-1; b=1;
  gamma=[1 10 100 500]; 
  disp('--------------------------------------')
  disp('Test 4');
  for i = 1:length(n)
    for j = 1:4
      g = @(x) testfunktion4(x, gamma(j));
      y_arctan(i,j) = multilevel_quad(g, a, b, n(i), Nmax, tol);
      error_arctan(i,j) = norm(y_arctan(i,j) - 2/gamma(j)*atan(gamma(j)) ); 
    end
  end
  
  % Plots
  figure
  semilogy(n, error_sqrt, 'bo-');
  
  figure 
  semilogy(n, error_sin, 'bo-');
  
  figure 
  semilogy(n, error_arctan(:,1), 'bo-');
  hold on
  semilogy(n, error_arctan(:,2), 'ro-');
  semilogy(n, error_arctan(:,3), 'go-');
  semilogy(n, error_arctan(:,4), 'ko-');
  hold off
end

function y = multilevel_quad(fun, a, b, n, Nmax, tol)
% Funktion zur Integration  anhand von Multilevel Quadratur
%  fun : Zu integrierende Funktion
%  a   : obere Intervallgrenze
%  b   : untere Intervallgrenze
%  n   : Anzahl der aequidistanten Intervalle zu Beginn
%  Nmax: Anzahl der maximal zulässigen Intervalle
%  tol : Toleranzgrenze

  % Array fuer lokale Fehler
  e_loc = zeros(1);
  % Initialisiere globalen Fehler
  e_glob = 1;
  % Stuetzstellen
  z = linspace(a,b,n);
  % Initialisiere Trapez und Simpson Summe
  Trapez = 0;
  Simpson = 0;
  
  k = 2;  
  % Fehler berechnung fuer jedes Teilintervall
    while (k <= length(z)) && (n <= Nmax) 
      
      % Maximumsstrategie fuer Schwellwert der lokalen Fehler
      nu = max(abs(e_loc(:)));
      
      h_k = z(k) - z(k-1);
      % I_k = [z(k),z(k-1)];
      % Trapez-Regel
      Trapez(k) = 1/2*( fun(z(k-1)) + fun(z(k)) )*h_k;
    
      % Simpson-Regel
      z_halb = 1/2*( z(k-1) + z(k));
      Simpson(k) = h_k/6*( fun(z(k-1)) + 4*fun(z_halb) + fun(z(k)));
    
      % lokale a posteriori Fehlerschaetzung
      e_loc(k) = Simpson(k) - Trapez(k);
   
      % Einhaltung des Schwellwerts fuer globalen Fehler
      if (abs(e_loc(k)) ) > nu && (abs(e_glob) > 1/2*tol)
        z(k+1:end+1) = z(k:end);
        z(k) = z_halb;
        n = n+1;
      else 
        k = k+1;
      end
    % globaler Fehler
    e_glob = sum(e_loc(:));
    end % while loop    
  
  %disp('Anzahl von Stuetzstellen am Ende:'); n
  y = sum(Simpson);
  
end

function y = testfunktion1(x)
  y = x;
end

function y = testfunktion4(x, gamma)
  y = 1 / (1 + (gamma*x)^2);
end

% Tatsächliche Stammfunktionen zum Vergleich der Fehler
function y = stamm_sqrt(x)
  y = 2/3*x^(3/2)
end

function y = stamm_gamma(x, gamma)
  y = 2/gamma * arctan(gamma)
end