function [A,rhs] = fd1d_bvp_test04 ( x )
 
  [A,rhs] = prob_assemble ( @a4, @a4prime, @b4, @c4, @f4, x );

  return
  
end



function value = a4 ( x )

  value = 1.0;

  return
end


function value = a4prime ( x )


  value = 0.0;

  return
end

function value = b4 ( x )

  value = 10^(4);

  return
end


function value = c4 ( x )

  value = 0.0;

  return
end


function value = f4 ( x )

  value = 1.0;

  return
end









