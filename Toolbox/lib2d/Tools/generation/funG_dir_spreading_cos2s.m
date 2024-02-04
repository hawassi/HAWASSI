function [Rtheta,StdTheta,cos2pdf,theta]=funG_dir_spreading_cos2s(Nomg,s,mDir_rad)
% ========================================================================%
% function for directional-spreading 'cos^2s'
% Homogeneous sea-state is aimed : 
%  - one wave component will only has one wave direction.
%  - the direction are chosen randomly based on 'cos^2s' distribution
% ========================================================================%

theta  = linspace(-pi/2,pi/2, Nomg ); 
Ntheta = length(theta);

% probability density function (PDF)
% the integral of the function should be 1
cos2pdf  = (gamma(s+1)./( sqrt(pi)*gamma(s+0.5))).* cos(theta).^(2*s) ; %% eq 6.3.24 Holthuijsen

% int_cos2pdf = trapz(theta,cos2pdf);    %the integration of the directional spreading function should be = 1
% disp(['the integral of the directional spreading function is ', num2str(int_cos2pdf)]);

% cummulative density function (CDF)
cos2cdf = zeros(size(cos2pdf));
total = 0;
for i=1:Ntheta
    total = total + cos2pdf(i);
    cos2cdf(i) = total;
end
% normalization
cos2cdf = cos2cdf./max(cos2cdf);

%% added function : with a very small sloping curve, to avoid problem where one variable gives more than one value
yadd = (exp(-12)/pi).*theta + exp(-12)/2;
cos2cdf = cos2cdf + yadd;

% inverse of CDF
% now we generate a uniform random number [0,1] then generate
% random number from cos2s-distribution :
NR = Nomg;
rng('shuffle');pause(0.1); %%shuffle random generator; 
RN = rand([NR,1]);

% Rcos2s  = interp1(cos2cdf,theta- mDir_rad,RN,'spline');  %random phase which has distribution cos^2
Rcos2s  = interp1(cos2cdf,theta,RN,'spline');  %random phase which has distribution cos^2

%% calculating standard deviation
std_cos2pdf  = std(Rcos2s);
mean_cos2pdf = mean(Rcos2s)+mDir_rad;

% saving output
Rtheta   = Rcos2s'  ;           %random theta \in [-pi/2, pi/2];
StdTheta = std_cos2pdf;         % standard deviation of the distribution in [Rad/s]
DTheta   = cos2pdf';            % PDF 

% Ntheta   = length(theta);       % should be the same with Nt
% dtheta   = theta(2) - theta(1); % delta theta