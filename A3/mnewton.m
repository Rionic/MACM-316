function [root, iter, xlist] = mnewton( func, pfunc, ppfunc, xguess, tol )
if nargin < 4
  fprintf(1, 'NEWTON: must be called with at least four arguments' );
  error( 'Usage:  [root, niter, xlist] = newton( func, pfunc, xguess, [tol], [mult] )' );
end
if nargin < 5, tol  = 1e-6; end
func = fcnchk( func );
pfunc= fcnchk( pfunc );
x    = xguess;
fx   = feval( func,  x );
fpx  = feval( pfunc, x );
fppx = feval( ppfunc, x);
if( fx == 0 | fpx == 0 ) 
  error( 'NEWTON: both f and f'' must be non-zero at the initial guess' );
end
xlist= [ x ];
done = 0;
iter = 0;
while( ~done )
  x0  = x;
  x   = x0 - (fx * fpx) / (fpx^2 - (fx*fppx));
  fx  = feval( func,  x );
  fpx = feval( pfunc, x );
  fppx = feval(ppfunc, x);
  if( abs(x-x0) < tol )     % absolute tolerance on x
    done = 1;
  else
    xlist = [ xlist; x ];   % add to the list of x-values 
    iter  = iter + 1;
  end
end

root = x;
%END newton.
