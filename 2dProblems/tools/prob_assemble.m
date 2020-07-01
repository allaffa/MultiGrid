function [A,rhs] = prob_assemble ( d, b1, b2, c, f, x, y )

%*****************************************************************************80
%
%% FD2D_HEAT_STEADY solves the steady 2D heat equation.
%
%  Discussion:
%
%    Nodes are assigned a singled index K, which increases as:
%
%    (NY-1)*NX+1  (NY-1)*NX+2  ...  NY * NX
%           ....         ....  ...    .....
%           NX+1         NX+2  ...   2 * NX
%              1            2  ...       NX
%
%    Therefore, the neighbors of an interior node numbered C are
%
%             C+NY
%              |
%      C-1 --- C --- C+1
%              |
%             C-NY
%
%    Nodes on the lower boundary satisfy:
%      1 <= K <= NX
%    Nodes on the upper boundary satisfy:
%      (NY-1)*NX+1 <= K <= NY * NX
%    Nodes on the left boundary satisfy:
%      mod ( K, NX ) = 1
%    Nodes on the right boundary satisfy:
%      mod ( K, NX ) = 0
%
%    If we number rows from bottom I = 1 to top I = NY
%    and columns from left J = 1 to right J = NX, we have
%      K = ( I - 1 ) * NX + J
%    and
%      J = 1 + mod ( K - 1, NX )
%      I = 1 + ( K - J ) / NX
%      
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    23 July 2013
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer NX, NY, the number of grid points in X and Y.
%
%    Input, real X(NX), Y(NY), the coordinates of grid lines.
%
%    Input, function pointer @D(X,Y), evaluates the thermal conductivity.
%
%    Input, function pointer @F(X,Y), evaluates the heat source term.
%
%    Output, real U(NX,NY), the approximation to the solution at the grid points.
%

%
%  Set the total number of unknowns.
%
  nx = size(x,1);
  ny = size(y,1);
  n = nx * ny;
%
%  Set up the matrix as a sparse quantity.
%
%  A = sparse ( [], [], [], n, n, 5 * n );
  A = sparse ( [], [], [], n, n, 9 * n );
  rhs = zeros ( n, 1 );
%
%  Define the matrix at interior points.

  [ A, rhs ] = interior9_adaptive ( nx, ny, x, y, d, b1, b2, c, f, A, rhs );
%
%  Handle boundary conditions.
%
  [ A, rhs ] = boundary ( nx, ny, x, y, A, rhs );


  return
end
