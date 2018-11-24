function [a, s2, R2, dv] = ma_reg(X,y)
%function [a, s2, R2, dv] = ma_reg(X,y)
% permet d’intégrer les calculs faits précedemment
%a = vecteur des coefficients du modèle
%s2 = variance des erreurs
%R2 = coefficient détermination
%dv = tableau contenant résidu, levier et contribution

[n,p] = size(X);
[m,q] = size(y);
if n ~= m
    error('X et y doivent avoir le même nombre de lignes')
end
a = (X'*X)\(X'*y);

e = y-X*a; 
s2 = e'*e/(n-p); 
R2 = 1 - e'*e/sum((y-mean(y)).^2);
h =diag(X*((X'*X)\(X'))); 
c = h./(1-h).^2/p.*e.^2/(e'*e);

dv = [e h c]; 
end