function [A,rhs] = fd2d_heat_steady_test03 ( x, y )

%*****************************************************************************80
%
%% FD2D_HEAT_STEADY_TEST01 demonstrates the use of FD2D_HEAT_STEADY.
%
%  Specify the spatial grid.
%
%  Solve the finite difference approximation to the steady 2D heat equation.
%
  [A,rhs] = prob_assemble ( @d, @b1, @b2, @c, @f, x, y );
%

  return
end


function value = d ( x, y )

%*****************************************************************************80
%
%% D evaluates the heat conductivity coefficient.

  value = 1.0 + 0.5*x + 0.5*y;

  return
end

function value = b1 ( x, y )

%*****************************************************************************80
%
%% B1 evaluates the convective field in the x-direction.

  value = 0.0 * x + 0.0 * y;

  return
end

function value = b2 ( x, y )

%*****************************************************************************80
%
%% B2 evaluates the convective field in the x-direction.

  value = 0.0 * x + 0.0 * y;
  
  return
end

function value = c ( x, y )

%*****************************************************************************80
%
%% C evaluates the reactive coefficient.

  value = 0.0 + 0.0*x + 0.0*y;

  return
end



function value = f ( x, y )

%*****************************************************************************80
%
%% F evaluates the heat source term.

  value = 1.0 + 0.0*x + 0.0*y;

  return
end
