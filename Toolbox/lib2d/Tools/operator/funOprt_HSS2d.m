% HybridSpatialSpectral operator; 
% from fourier function z to fourier function Az; a is space interp function
% iftz = ifft(z,'symmetric'); requires 3 FT's
function Az = funOprt_HSS2d(A,a,z,iftz)
    Az = ( fft2(a.*funC_ifft2(A.*z)) + A.*fft2(a.*iftz) )/2;
end