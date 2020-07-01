function [A,rhs] = fd1d_bvp_test02 ( x )

  [A,rhs] = prob_assemble ( @a2, @a2prime, @b2, @c2, @f2, x );

  return
  
end

function value = a2 ( x )

  value = 1.0 + x .* x;

  return
end


function value = a2prime ( x )

  value = 2.0 * x;

  return
end

function value = b2( x )

value = 0.0;

return 
end

function value = c2 ( x )

  value = 2.0;

  return
end

function value = f2 ( x )

  value = x .* ( 5.0 - x ) .* exp ( x );

  return
end

function value = exact2 ( x )

  value = ...
      ( x .* ( 1.0 - x ) .* exp ( x ) )          .* ( x <= 2.0 / 3.0 ) ...
    + ( x .* ( 1.0 - x )  * exp ( 2.0 / 3.0 ) )  .* (      2.0 / 3.0 < x );

  return
end