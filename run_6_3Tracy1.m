function err = run_6_3()
  x = linspace(-1,1,5);
  err = plotspline(x);
  
end

function err = plotspline(x)
  svec = zeros(20, length(x));
  fvec = zeros(20, length(x));
  y_all = zeros(20, length(x));
  
  i = 1;
  while (i < length(x))
    y = linspace(x(i), x(i+1), 20);
    y_all(:,i) = y;
    % Teilintervalle
    for j = 1:length(y)
      svec(j,i) = func(x(i)) * (y(j)-x(i+1))/(x(i) - x(i+1)) + func(x(i+1))* (y(j) - x(i))/(x(i+1)-x(i));
      fvec(j,i) = func(y(j));  
    end
    err(i) = max(abs(svec(:,i)-fvec(:,i)));
    i = i + 1;
  end %while
   err = max(err);
  
  figure
  plot( y_all(:,1), svec(:,1), y_all(:,1) , fvec(:,1));
  hold on;
  for i = 2:length(x)
    plot( y_all(:,i), svec(:,i), y_all(:,i) , fvec(:,i));
  end
  hold off;
  
end

function y = func(x)
   y = 1/(1+25*x**2);
end