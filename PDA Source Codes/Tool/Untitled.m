% 1、直线拟合的matlab代码

% Fitting a best-fit line to data, both noisy and non-noisy
x = rand(1,10);
n = rand(size(x)); % Noise
y = 2*x + 3; % x and y satisfy y = 2*x + 3
yn = y + n; % x and yn roughly satisfy yn = 2*x + 3 due to the noise
% Determine coefficients for non-noisy line y=m1*x+b1
Xcolv = x(:); % Make X a column vector
Ycolv = y(:); % Make Y a column vector
Const = ones(size(Xcolv)); % Vector of ones for constant term
Coeffs = [Xcolv Const]\Ycolv; % Find the coefficients
m1 = Coeffs(1);
b1 = Coeffs(2);
% To fit another function to this data, simply change the first 
% matrix on the line defining Coeffs
% For example, this code would fit a quadratic 
% y = Coeffs(1)*x^2+Coeffs(2)*x+Coeffs(3)
% Coeffs = [Xcolv.^2 Xcolv Const]\Ycolv; 
% Note the .^ before the exponent of the first term
% Plot the original points and the fitted curve
figure
plot(x,y,'ro')
hold on
x2 = 0:0.01:1;
y2 = m1*x2+b1; % Evaluate fitted curve at many points
plot(x2, y2, 'g-')
title(sprintf('Non-noisy data: y=%f*x+%f',m1,b1))
% Determine coefficients for noisy line yn=m2*x+b2
Xcolv = x(:); % Make X a column vector
Yncolv = yn(:); % Make Yn a column vector
Const = ones(size(Xcolv)); % Vector of ones for constant term
NoisyCoeffs = [Xcolv Const]\Yncolv; % Find the coefficients
m2 = NoisyCoeffs(1);
b2 = NoisyCoeffs(2);
% Plot the original points and the fitted curve
figure
plot(x,yn,'ro')
hold on
x2 = 0:0.01:1;
yn2 = m2*x2+b2;
plot(x2, yn2, 'g-')
title(sprintf('Noisy data: y=%f*x+%f',m2,b2))