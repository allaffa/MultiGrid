function [A,rhs] = fd1d_bvp_test01 ( x )

  [A,rhs] = prob_assemble ( @a1, @a1prime, @b1, @c1, @f1, x );

  return
  
end



function value = a1 ( x )

  value = 1.0;

  return
end


function value = a1prime ( x )


  value = 0.0;

  return
end

function value = b1 ( x )

  value = 0.0;

  return
end


function value = c1 ( x )

  value = 0.0;

  return
end


function value = f1 ( x )

  value = x .* ( x + 3.0 ) .* exp ( x );

  return
end


function value = exact1 ( x )

  value = x .* ( 1 - x ) .* exp ( x );

  return
end







