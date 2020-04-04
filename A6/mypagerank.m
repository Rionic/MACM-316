% Estimate the PageRank for a small web connectivity
% matrix using the power method.

p = 0.98;  % "teleport" probability

% Set up the connectivity matrix using index
% vectors that define the from-to connections
% between pages numbered 1 to n.
n = 18;      % number of pages

ii = [5 11 2 3 1 7 9 6 10 8 4 5 14 1 8 9 8 10 12 13 13 9 6 12 8 7 7 13 14 18 17 11 12 15 12 13 16 13 12 13 17]; % which page jj is connecting to
jj = [13 1 1 2 3 3 3 4 4 4 5 6 6 7 7 8 9 9 9 9 10 10 10 11 11 11 12 12 13 14 14 15 15 16 16 16 17 17 17 18 18]; % #outlinks
G = sparse(ii, jj, 1, n, n);

c = full(sum(G));  % column sums
e = ones(n,1);
k = find(c~=0);    % indices of zero columns 
D = sparse(k, k, 1./c(k), n, n);
z = ((1-p)*(c~=0) + (c==0)) / n;

A = p*G*D + e*z;   % PageRank matrix

x = e/n;           % initial guess
xold = zeros(n,1);
niter = 0;
while norm(x-xold) > 0.0001,
  niter = niter + 1;
  xold = x;
  x = A*x;
end
x, niter
bar(x), shg    % plot ranks as a bar plot
title("p = 0.98 with extra link")