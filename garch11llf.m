% log likelihood function
% GARCH( 1 , 1 )
function llf = garch11llf( para, r, T, eps, m, sigma )
% c = para( 1 );
c = 100 * para( 1 );
rho = [0 0];%para( 2 : 3 );
% w = para( 4 ); 
w = 100^2 * para( 4 );
g = para( 5 ) ;
beta1 = para( 6 ) ;
% eps = zeros( T + 3 , 1 ) ;
eps( 3 , : ) = para( 7 ) ;
% sigma = zeros( T + 2 , 1 ) ;
sigma( 1 : 2 ) = para( 8 : 9 ) ; 
if length( r ) == T
    r = [ 0 ; para( 10 ) ; para( 11 ) ; r ] ;
else
    r = [ 0 ; para( 10 ) ; para( 11 ) ; r( 4 : end ) ] ;
end
% pv = para(12);
pv = 0.01 * para(12);
% m = zeros( T + 2 , 1 ) ;
for t = 1 : 1 : ( T + 2 )
    if t < 3 
        m( t ) = c + pv * sigma( t );
    else
        sigma( t ) =min(1e5, w + g * eps( t )^2 + beta1 * sigma( t - 1 ));
        m( t ) = c + pv * sigma( t ) ;
        eps( t + 1 ) = r( t + 1 ) - m( t ) - rho( 1 ) * ( r( t ) - m( t - 1 ) ) ...
            - rho( 2 ) * ( r( t - 1 ) - m( t - 2 ) );
    end
end
% strip off the initial values
% eps = eps( 3 : end , : ) ;
% sigma = sigma( 3 : end , : ) ;
Loglikelihoods = log( exp( -0.5 * ( eps( 4 : end , : ).^2 ) ./ sigma( 3 : end , : ) ) ./ sqrt( 2 * pi * sigma( 3 : end , : ) )) ;
% Loglikelihoods = -0.5 * ( eps.^2 ) ./ sigma - 0.5 * log( ( 2 * pi * sigma ) ) ;
llf = -sum(max(-1e10, Loglikelihoods));
end
