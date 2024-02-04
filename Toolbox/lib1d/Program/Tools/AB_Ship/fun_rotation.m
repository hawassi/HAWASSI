function [Xnew, Znew]=fun_rotation(X,Z,Sthetac,Sxc)
X=X-Sxc;
Xnew=Sxc+X.*cos(Sthetac)-Z.*sin(Sthetac);
Znew=X.*sin(Sthetac)+Z.*cos(Sthetac);