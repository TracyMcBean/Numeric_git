% Uebung 4, Numerik 1
% David, Tracy
% Aufgabe 1

function result = run_1_4()
    %RUN_1_4() loest lineares Ausgleichsproblem fuer verschieden Gitter
    
    disp('Test 1, Gitter mit m = 2 und n = 1:')
    result = givens_rotation(2,1);
    fprintf('Test result: %s \n', result)
    
    disp('----------------------------------------------')
    disp('Test 2, Gitter mit m = 4 und n = 2:')
    result =  givens_rotation(4,2);
    fprintf('Test result:\n')
    disp(result)
    
    disp('----------------------------------------------')
    disp('Test 3, Gitter mit m = 10 und n = 3:')
    result =  givens_rotation(10,3);
    fprintf('Test result:\n')
    disp(result)
    
end %function

function result = givens_rotation(m, n)
    %GIVENS_ROTATION() fuehrt Givens Rotation aus um lin. Ausgleichsproblem
    %zu loesen.
    % n: Grad vom Polynom
    % m: Gitteraufloesung
    
    % Initialisiere
    x = zeros(m, 1);
    b = zeros(m, 1);
    xcooef = zeros(m,1);  % fuer matrix A gebraucht
    A = zeros(m,n);
    x(m)= 2*pi;
    xcooef(:) = 1;
    
    % Falls mehr als 2 Gitterpunkte setze die Abstaende entsprechend.
    if ( m > 2)
        step = 2*pi / m;
        for i = 2:(m-1)
            x(i) = x(i-1) + step;
        end %for
    end % if
 
    for i = 1:m
        b(i) = sin(x(i));
    end 
    
    % Matrix A auffuellen mit t
    for i = 1:n
        A(:, i) = xcooef;
        xcooef = xcooef.* x;
    end 
    
    % Erstelle 3D-Matrix mit Einheitsmatrix
    Q_single = eye(n)
    for k = 1:n
        Q(:,:,k) = Q_single
    end 
    
    % Berechne die
    
    result = x;
end %function