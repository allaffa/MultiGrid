%% Define the elliptic problem

function [A,rhs] = exact_solution_test ( x, y )

  nx = size(x,1);
  ny = size(y,1);
  
  x1 = x(1);
  x2 = x(end);
  y1 = y(1);
  y2 = y(end);

  %Define the coordinates for the centers of the gaussians

  %top-left
  x_c1 = 1/4*( x2-x1 );
  y_c1 = 3/4*( y2-y1 );

  %bottom-right
  x_c2 = 3/4*( x2-x1 );
  y_c2 = 1/4*( y2-y1 );

  %standard deviation
  dev1 = 0.02;
  dev2 = 0.02;

  g1 = @(x,y) exp( -( (x-x_c1)^2 + (y-y_c1)^2)/(dev1^2) );
  g2 = @(x,y) exp( -( (x-x_c2)^2 + (y-y_c2)^2)/(dev2^2) );
  
  %Define the exact solution as a linear combination of the gaussians

  %coefficients for linear combination of gaussians
  c1 = 2;
  c2 = 3;

  uex = @(x,y) c1*g1(x,y) + c2*g2(x,y);
  
  %Definition of the lambda function for the right hand side

  function value = f ( x, y )
           value = - 2 * c1 * (-2 * (x-x_c1)^2 + dev1)/( dev1^2 ) * g1(x,y) + ...
           - 2 * (-2 * c1 * (y-y_c1)^2 + dev1)/( dev1^2 ) * g1(x,y) + ...
           - 2 * (-2 * c2 * (x-x_c2)^2 + dev2)/( dev2^2 ) * g2(x,y) + ...
           - 2 * (-2 * c2 * (y-y_c2)^2 + dev2)/( dev2^2 ) * g2(x,y);
          
           value = - value;
       return;
  end
  
  %
  %  Solve the finite difference approximation to the steady 2D heat equation.
  %
  [A,rhs] = prob_assemble ( @d, @b1, @b2, @c, @f, x, y );
  %

  return
end

%% D evaluates the heat conductivity coefficient.

function value = d ( x, y )

  value = 1.0 + 0.0*x + 0.0*y;

  return
end

%% B1 evaluates the convective field in the x-direction

function value = b1 ( x, y )

  value = 0.0 * x + 0.0 * y;

  return
end

%% B2 evaluates the convective field in the x-direction

function value = b2 ( x, y )

  value = 0.0 * x + 0.0 * y;
  
  return
end

%% C evaluates the reactive coefficient.

function value = c ( x, y )

  value = 0.0 + 0.0*x + 0.0*y;

  return
end
