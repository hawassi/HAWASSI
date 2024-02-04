% HybridSpatialSpectral operator; 
% from fourier function z to fourier function Az; a is space interp function
% iftz = ifft(z,'symmetric'); requires 3 FT's
function Az = HSS(A,a,z,iftz)
    Az = ( fft(a.*Ifft(A.*z)) + A.*fft(a.*iftz) )/2;
end