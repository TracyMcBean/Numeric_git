% Uebung 2 Aufgabe 1, Numerik 1 SS18
% Theresa Kiszler, David
% Abgabe: 30.04

% Function to run the iteration with x_0 and y_0
function x = run_1_2()
  %RUN_1_2() use the newton method for rootfinding 
  % define roots
  root1 = 0;
  root2 = pi;
  
  x_0 = -1;
  fprintf('first run with x_0: %d \n', x_0);
  newton_method(x_0, root1);
  disp(' --------------------------------');
  newton_method(x_0, root2);
  
  y_0 = 4;
  disp(' --------------------------------');
  fprintf('second run with y_0: %d \n', y_0);
  newton_method(y_0, root1);
  disp(' --------------------------------');
  newton_method(y_0, root2);
   
endfunction 

% Newton method
function x_n = newton_method(var_in, x_star)
  %NEWTON_METHOD(var_in) define what var_in is
  
  % Initialize
  converge = 0;
  x_n = var_in;
  n = 1;
  error_vec = [0];
  
  
  % check convergence criteria
  while (converge == 0)
    % Error calculation
    error = abs(x_n - x_star);
    error_vec(n) = error;
    
    n = n + 1;
    
    f = x_n*(x_n - pi)^2; 
    Df = (x_n - pi)^2 - 2*x_n*(x_n - pi);
    dx_n = -f/Df;
    x_nplus1 = x_n + dx_n;
   
    % check break criteria
    dx_n_quer = - (x_nplus1*(x_nplus1 - pi)^2)/ Df; 
    % fprintf(' Values: \n f = %d, Df = %d, x_n = %d, x_nplus1 = %d, dx_n = %d , dx_n_quer: %d \n', f, Df, x_n, x_nplus1, dx_n, dx_n_quer);
    if (norm(dx_n_quer) > norm(dx_n) | n > 20)
      fprintf('Success not possible. \n');
      x_n = NA;
      break
    endif
    
    % set new convergence
    if (norm(x_n-x_star) <= norm(dx_n))
      converge = 1;
    endif
    
    % go one step up
    x_n = x_nplus1;

  endwhile
  fprintf(' Solution is: %d \n', x_n); 

  % Plotting the error
  figure
  semilogy([1:length(error_vec)], error_vec(:))
  
endfunction


