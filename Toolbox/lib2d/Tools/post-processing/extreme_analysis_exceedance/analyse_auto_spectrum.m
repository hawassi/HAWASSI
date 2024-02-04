function [omega,S]=analyse_auto_spectrum(time,x,domega)

% syntax:       [omega,S]=analyse_auto_spectrum(time,x,domega)
% example:      [omega,S]=analyse_auto_spectrum(time,x)
% description:  Function analyse_auto_spectrum estimates the
%               spectal density of the signal x(time), using
%               Welch method (half overlapping hanning windows). The actual
%               work is outsourced to analyse_cross_spectrum
%
% input:        - the time axis of signal x(time) [s]
%               - the signal x
%               - optional, the desired domega stepp. If not specified
%               0.02rad/s will be used
%
% output        - omega [rad/s], the frequency axis
%               - S, the cross spectral density
%
% revisions:    Nr.     Date        By              Version
%                0    May 28 2010   JvdB               1.0 


if nargin < 3, domega= 0.02; end

[omega,S]=analyse_cross_spectrum(time,x,x,domega);
S=real(S); %imaginary part should already be zero, but there may be some small residual