%Performs one step of Runge-Kutta.  4th order is the default.

function U = rk( t, U, k, odeFun, stages )

if stages == 1
    U = U + k*odeFun(t,U);                       %1st order Runge-Kutta (Forward Euler)
elseif stages == 2
    s1 = odeFun( t, U );                         %2nd order Runge-Kutta, step1
    s1 = odeFun( t+k/2, U+k/2*s1 );              %step 2
    U  = U + k*s1;                               %Combine for new value
elseif stages == 3
    s1 = odeFun( t, U );                         %3rd order Runge-Kutta, step 1
    s2 = odeFun( t+k/3, U+k/3*s1 );              %step 2
    s2 = odeFun( t+2*k/3, U+2*k/3*s2 );          %step 3
    U  = U + k/4 * ( s1 + 3*s2 );                %Combine for new value
else
    s1 = odeFun( t, U );                         %4th order Runge-Kutta, step 1
    s2 = odeFun( t+k/2, U+k/2*s1 );              %step 2
    s3 = odeFun( t+k/2, U+k/2*s2 );              %step 3
    s4 = odeFun( t+k, U+k*s3 );                  %step 4
    U  = U + k/6 * ( s1 + 2*s2 + 2*s3 + s4 );    %Combine for new value
end