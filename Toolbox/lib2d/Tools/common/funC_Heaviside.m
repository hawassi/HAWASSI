function Heav=funC_Heaviside(input)
Heav = zeros(size(input));
Heav(input >= 0) = 1;
