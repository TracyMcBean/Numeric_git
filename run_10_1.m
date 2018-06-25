% Numerik 1 SS18
% Uebung 10, Aufgabe 1
% David und Tracy

function [] = run_10_1()
% Test multilevel Quadratur Funktion
  a=0; b=1; 
  n = 10;
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
 
  [y_sqrt,simpsqrt] = multilevel_quad(@sqrt, a, b, n, Nmax, tol);
  error_sqrt = abs(simpsqrt - 2/3);
  nsq = n:(n+length(simpsqrt)-1);
  disp(y_sqrt)
  
  
  % Test 3 y = sin(x)
  a=0; b=pi; 
  disp('--------------------------------------')
  disp('Test 3: sin(x)');
  
    [y_sin,simpsin] = multilevel_quad(@sin, a, b, n, Nmax, tol);
    error_sin = abs(simpsin - 2);
    nsi = n:(n+length(simpsin)-1);
    disp(y_sin)
  
  % Test 4 arctan
  a=-1; b=1;
  gamma=[1 10 100 500]; 
  disp('--------------------------------------')
  disp('Test 4');
  
  
      g1 = @(x) testfunktion4(x, gamma(1));
      [y_arctan1,simptan1] = multilevel_quad(g1, a, b, n, Nmax, tol);
      error_arctan1 = abs(simptan1 - 2/gamma(1)*atan(gamma(1)) );
      nta1=n:(n+length(simptan1)-1);
      disp(y_arctan1)
      
      g2 = @(x) testfunktion4(x, gamma(2));
      [y_arctan2,simptan2] = multilevel_quad(g2, a, b, n, Nmax, tol);
      error_arctan2 = abs(simptan2 - 2/gamma(2)*atan(gamma(2)) );
      nta2=n:(n+length(simptan2)-1);
      disp(y_arctan2)
      
      g3 = @(x) testfunktion4(x, gamma(3));
      [y_arctan3,simptan3] = multilevel_quad(g3, a, b, n, Nmax, tol);
      error_arctan3 = abs(simptan3 - 2/gamma(3)*atan(gamma(3)) );
      nta3=n:(n+length(simptan3)-1);
      disp(y_arctan3)
      
      g4 = @(x) testfunktion4(x, gamma(4));
      [y_arctan4,simptan4] = multilevel_quad(g4, a, b, n, Nmax, tol);
      error_arctan4 = abs(simptan4 - 2/gamma(4)*atan(gamma(4)) );
      nta4=n:(n+length(simptan4)-1);
      disp(y_arctan4)
      
      
  
  
  
  % Plots
  figure
  semilogy(nsq, error_sqrt, 'bo-');
  xlabel('n');
  ylabel('Fehler');
  title('f(x)=sqrt(x), in [0,1]')
  figure 
  semilogy(nsi, error_sin, 'bo-');
  xlabel('n');
  ylabel('Fehler'); 
  title('f(x)=sin(x) in [0,pi]')  
  figure 
  semilogy(nta1, error_arctan1, 'bo-');
  hold on
  semilogy(nta2, error_arctan2, 'ro-');
  semilogy(nta3, error_arctan3, 'go-');
  semilogy(nta4, error_arctan4, 'ko-');
  hold off
  xlabel('n');
  ylabel('Fehler');
  title('f(x)=1/(1+(gamma*x)^2) in [-1,1]')
  legend('Gamma=1','Gamma=10','Gamma=100','Gamma=500');
end

function [y,Simp] = multilevel_quad(fun, a, b, n, Nmax, tol)
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
  Simp = n:Nmax;
  
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
      Simpson(k) = h_k*( fun(z(k-1)) + 4*fun(z_halb) + fun(z(k)))/6;
    
      % lokale a posteriori Fehlerschaetzung
      e_loc(k) = Simpson(k) - Trapez(k);
   
      % Einhaltung des Schwellwerts fuer globalen Fehler
      if (abs(e_loc(k)) ) > nu && (abs(e_glob) > 1/2*tol)
        z(k+1:end+1) = z(k:end);
        z(k) = z_halb;
        n = n+1;
      else 
        Simp(k)=sum(Simpson);  
        k = k+1;
      end
      z = linspace(a,b,n);
    % globaler Fehler
    e_glob = sum(e_loc(:));
    Simp = Simp(1:k-1);
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