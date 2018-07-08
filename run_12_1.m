% Numerik 1, SS18
% Uebung 12 Aufgabe 1b
% Tracy, David

function []=run_12_1()
  
  tspan = [0,1000]; 
  x0 = [0,1];
  dt = 0.01;      %tau
    
  % Initialisiere f(x0)
  EulerExp(1,:) = x0; EulerImp(1,:) = x0;
  t(1) = tspan(1);
  k = 2;
  
  Df = [0,-1; 1,0];
  
% Berechne X-Werte für Zeitschritte
  while (t(k-1) <= tspan(2))
    t(k) = t(1) + k*dt;
% Expliziter Euler
    EulerExp(k,1) = EulerExp(k-1,1) + dt * f(EulerExp(k-1,:), t(k-1))(1);
    EulerExp(k,2) = EulerExp(k-1,2) + dt * f(EulerExp(k-1,:), t(k-1))(2);
    
% Impliziter Euler
% Aufgrund des Rechenaufwands und der gleichzeitig schnellen Konvergenz 
% gegen 0 werden hier nur die ersten 20 Schritte berechnet.
    if (t(k-1) < 21)
      x = EulerImp(k-1,:); 
      [EulerImp(k,:), fval, info] = fsolve(@(xNext) x + f(xNext, t(k)) - xNext, x); 
    end
    
    k = k + 1;
  end
  
  % Fuer Beide in einem Plot
%  bothEulerx(:,1) = EulerExp(:,1);
%  bothEulerx(:,2) = EulerImp(:,1);
%  bothEulery(:,1) = EulerExp(:,2);
%  bothEulery(:,2) = EulerImp(:,2);
%  tnew = [t(:),t(:)];
  
  figure
  % plot3 (tnew, bothEulerx, bothEulery);
  plot3 (t(:), EulerExp(:,1), EulerExp(:,2));
  xlabel ("Time t");
  ylabel ("x");
  zlabel ("y");
  title ("ODE solved with explicit Euler (timestep: 0.01)")

  figure 
  plot3 (t(1:20), EulerImp(1:20,1), EulerImp(1:20,2));
  xlabel ("Time t");
  ylabel ("x");
  zlabel ("y");
  title ("ODE solved with implicit Euler (timestep: 0.01)")


end

function dx = f(x, t)
  dx = zeros(2,1);
  dx(1) = -x(2);
  dx(2) = x(1);
end
