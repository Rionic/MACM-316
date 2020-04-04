% GERANDOM:  Gaussian Elimination for a random matrix.
% Compute the mean solution error over 'Mtry' trials for
% the system A*x = b, where A is a random NxN matrix and 
% x is an N-vector (containing only 1's).

% YOU NEED MODIFY THIS CODE FOR COMPUTING ASSIGNMENT 4

N    = 40000;             % Matrix size
Mtry = 100;            % Number of trials -- CHANGE THIS!
errs = zeros(Mtry, 1); % Vector of errors
x = ones(N,1);         % Exact solution
for i = 1 : Mtry
  A = 2*rand(N,N)-1;   % Random NxN matrix with entries in [-1, 1] 
  b = A*x;             % Right-hand side vector
  y = A \ b;           % Approximate solution from GE
  errs(i) = max(abs(y-x)); % Max-norm error in y
end

% Compute mean and standard deviation of the error
mean_err = mean(errs)
sdev_err = sqrt(var(errs));
median_err = median(errs)
% Plot the errors on a semi-log scale
semilogy(1:Mtry, errs)
title(['Error from ' num2str(Mtry) ' solves with a ', ...
       num2str(N) 'x' num2str(N) ' random matrix'])
xlabel('M (trial number)')
ylabel('Error')
grid on, shg