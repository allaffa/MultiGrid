function [A,rhs] = fd1d_bvp_test03 ( x )
 
  [A,rhs] = prob_assemble ( @a3, @a3prime, @ b3, @c3, @f3, x );

  return
  
end

function value = a3 ( x )

  value = ...
      100 * ( 1.0 + x .* x )   .*  ( x <= 1.0 / 3.0 ) ...
    + 0.01 * ( x + 7.0 / 9.0 )  .*  (      1.0 / 3.0 < x );

  return
end

function value = a3prime ( x )

  value = ...
      100 * ( 2.0 * x )  .*  ( x <= 1.0 / 3.0 ) ...
    + 0.01 * ( 1.0 )      .*  (      1.0 / 3.0 < x );

  return
end

function value = b3 ( x )

  value = 0.0;

  return
end

function value = c3 ( x )

  value = 2.0 * x;

  return
end

function value = f3 ( x )

  value = - x .* ( 2.0 * x .* x - 3.0 * x - 3.0 ) .* exp ( x );

  return
end

function value = exact3 ( x )

  value = ...
      ( x .* ( 1.0 - x ) .* exp ( x ) )  .* ( x <= 2.0 / 3.0 ) ...
    + ( x .* ( 1.0 - x ) )               .*      ( 2.0 / 3.0 < x );

  return
end

