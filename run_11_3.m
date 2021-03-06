% Numerik 1, SS18
% Uebung 11 Aufgabe 3
% Tracy, David

function [Y_Euler]=run_11_3()
  % Teste Euler:
  tspan = [0., 2*pi]; % wegen exakter LSG gekuerzt
  x0 = [1,0];
  
  i = 1;
  for N = 5:10:1000
    dt(i) = (tspan(2)-tspan(1)) / N; 
    
    Y_Euler = EulerEx(@f, tspan, x0, N);
  
    % Teste Taylor Entwicklung 2. Ordnung
    Df = [ -1, 1];
    Y_Taylor = Taylor(@f, Df, tspan, x0, N);
  
    % Exakte Loesung:
    Y_Exakt = Exakt(@F_exact, tspan, x0, N);
  
    error_Euler(i,:) = mean( abs(Y_Exakt - Y_Euler) );
    error_Taylor(i,:) = mean( abs(Y_Exakt - Y_Taylor) );
  
    i = i+1;
  end
  figure
  %plot (Y_Euler(:,1), Y_Euler(:,2));
  hold on
  %plot (Y_Euler(:,1), Y_Euler(:,3));
  plot (dt(:), error_Taylor(:,2), 'g--');
  plot (dt(:), error_Euler(:,2), 'r--');
  xlabel('tau')
  ylabel('Fehler')
  legend('Fehler nach Taylor','Fehler nach Euler')
  hold off
end

% (a) Implementiere Euler -----------------------------------------------------%

function [Y]= EulerEx(f, tspan, x0, N)
  % Funktion zur Berechnung von der zeitlichen Entwicklung der Funktion f mit 
  % der Zeit.
  % f    : Funktion
  % tspan: Zeitspanne
  % x0   : Anfangswert fuer x
  % N    : Anzahl der Schritte
  
  % Ueberpruefe ob N ein integer ist
  if (floor(N)==N)
    dt = (tspan(2)-tspan(1)) / N; % Tau
  else
    disp('Schritte muss ein int sein.');
    return;
  end
  
  % Initialisiere f(x0)
  x(1,:) = x0;
  t(1) = tspan(1);
  k = 2;
  
  % Berechne X-Werte für Zeitschritte
  while (t(k-1) <= tspan(2))
    t(k) = t(1) + k*dt;
    x(k,:) =  x(k-1,:).' + dt * f(x(k-1,:), t(k-1));
    k = k + 1;
  end
  
  t = t.';
  % Merge arrays
  Y = [t(:), x(:,1), x(:,2)];
  
end

% (b) Implementiere Taylor ----------------------------------------------------%

function [Y] = Taylor(f, Df, tspan, x0, N)
  % Berechnung der Funktionsentwicklung mit der Zeit anhand der Taylor Verfahren 
  % 2. Ordnung.
  % f    : Funktion
  % Df   : Ableitung von f
  % tspan: Zeitspanne
  % x0   : Anfangswert fuer x
  % N    : Anzahl der Schritte
  
  % Ueberpruefe ob N ein integer ist
  if (floor(N)==N)
    dt = (tspan(2)-tspan(1)) / N; % Tau
  else
    disp('Schritte muss ein int sein.');
    return;
  end
  
  % Initialisiere f(x0)
  x(1,:) = x0;
  t(1) = tspan(1);
  k = 2;
  
  % Berechne X-Werte für Zeitschritte
  while (t(k-1) <= tspan(2))
    t(k) = t(1) + k*dt;
    x(k,:) = x(k-1,:).' + dt * f(x(k-1,:), t(k-1)) + dt^2/2 * Df(:).'*f(x(k-1,:), t(k-1));
    k = k + 1;
  end
  
  t = t.';
  Y = [t(:), x(:,1), x(:,2)];
end

% Exakte Loesung --------------------------------------------------------------%
function Y = Exakt(f, tspan, x0, N)
  
  % Ueberpruefe ob N ein integer ist
  if (floor(N)==N)
    dt = (tspan(2)-tspan(1)) / N; % Tau
  else
    disp('Schritte muss ein int sein.');
    return;
  end
  
  x(1,:) = x0;
  t(1) = tspan(1);
  k = 2;
  
  % Berechne X-Werte für Zeitschritte
  while (t(k-1) <= tspan(2))
    t(k) = t(1) + k*dt;
    x(k,:) = f(t(k));
    k = k + 1;
  end
  
  t = t.';
  Y = [t(:), x(:,1), x(:,2)];
end

% Test funktion ---------------------------------------------------------------%

function dx = f(x, t)
  dx = zeros(2,1);
  dx(1) = -x(2);
  dx(2) = x(1);
end

function x = F_exact(t)
  x = zeros(2,1);
   x(1) = cos(t);
   x(2) = sin(t);
end