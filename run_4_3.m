cd /home/davidkamm/Numerik/A4

U = qr(rand(30));
V = qr(rand(30));
S = diag(2.^(-1:-1:-30));
A = U*S*V;

%Anwendung der Funktionen:
[Qhouse,Rhouse] = householder(A);
[Qgram,Rgram] = gs(A);

%Bestimmung der Fehler ||Qt*Q-I|| und ||QR-A||
FehlerQhouse = norm(transpose(Qhouse)*Qhouse-eye(length(A)))
FehlerQRhouse = norm(Qhouse*Rhouse-A)
FehlerQgram = norm(transpose(Qgram)*Qgram-eye(length(A)))
FehlerQRgram = norm(Qgram*Rgram-A)


function [Q,R] = householder(A)
%Die Funktion Zerlegt eine gegebene Matrix A in Qt und R, sodass gilt:
%Qt*A=R. Diese Funktion wirkt auf nxn Matrizen, fuer mxn muss sie
%geringfuegig modifiziert werden.

%Initialisiren
n = length(A);
B = A;

%Reservierung von Speicher fuer die n-1 Faktoren aus Qt und vorfestlegung
%dieser als Einheitsmatrizen
Qfaktoren = zeros(n,n,n-1);
for j =1:n-1
    Qfaktoren(:,:,j)= eye(n);
end


%Iterative bestimmung der n-1 Faktoren. Hier findet die eigentliche QR
%Zerlegung durch Householder-Reflexion statt. 
for i = 1:n-1
    I = eye(length(B));
    v = B(:,1)-norm(B(:,1))*I(:,1);
    Q = I - 2*(v*transpose(v))/(transpose(v)*v);
    Qfaktoren(i:n,i:n,i) = Q;
    C = Q*B;
    B = C(2:end,2:end);
end

%Multiplikation der Faktoren zur Bestimmung der Transponierten von Q
for k = 2:n-1
    Qfaktoren(:,:,k) = Qfaktoren(:,:,k)*Qfaktoren(:,:,k-1);
end

%Setzung der Ausgabewerte Qt und R 
Qt = Qfaktoren(:,:,k);
R = Qt*A;
Q = transpose(Qt);


%Der naechste Schritt setzt Werte die kleiner als eine vom Benutzer
%gewuenschte genauigkeit sind = 0. So werden Fehler durch die
%Maschienengenauigkeit als 0 dargestellt.

% genauigkeit = 10^(-8);
% 
% for g = 1:n
%     for  h = 1:n
%         if R(g,h) < genauigkeit
%             R(g,h)= 0;
%         end
%     end
% end

end

function [Q,R] = gs(A)
% Berechnet die QR-Zerlegung von A mittels des klassischen Gram-Schmidt-Verfahrens.

Q = zeros(size(A));
R = zeros(min(size(A)));

for j = 1:size(A,2)
    v = A(:,j);
    
    %for i = 1:j-1
    %    R(i,j) = Q(:,i)'*A(:,j);
    %    v = v - R(i,j)*Q(:,i);
    %end
    R(1:j-1,j) = Q(:,1:j-1)'*A(:,j);
    v = v - Q(:,1:j-1)*R(1:j-1,j);

    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);    
end
end
