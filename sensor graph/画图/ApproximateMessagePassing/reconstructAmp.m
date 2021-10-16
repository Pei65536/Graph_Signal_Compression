function xhat = reconstructAmp(A, y, T, tol)
% RECONSTRUCTAMP recovers a sparse vector x from few linear measurements y.
%
% xhat = reconstructAmp(A, y, T, tol, x, verbose)
%
%   Arguments:
%       A - measurement matrix
%       y - measurements
%       T - max number of iterations (optional)
%       tol - stopping criteria (optional)
%       x - original vector used to print progress of MSE (optional)
%       verbose - print progress optional
%
% Ulugbek Kamilov, LTHC, EPFL, 2010.

% Set some parameters
if(nargin < 3)
    T = 500;
end
if(nargin < 4)
    tol = 0.0001;
end

% Length of the original signal
N = size(A, 2);

% Length of the measurement vector
n = size(A, 1);

% Initial estimate
xhat = zeros(N, 1);
z = y;

% Start estimation
for t = 1:T
    % Pre-threshold value
    gamma = xhat + A'*z;

    % Find n-th largest coefficient of gamma
    threshold = largestElement(abs(gamma), n);

    % Estimate the signal (by soft thresholding)
    xhat = eta(gamma, threshold);

    % Update the residual
    z = y - A*xhat + (z/n)*sum(etaprime(gamma, threshold));

    % Stopping criteria
    if(norm(y - A*xhat)/norm(y) < tol)
        break;
    end
end