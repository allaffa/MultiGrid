function [ A, rhs ] = interior5 ( nx, ny, x, y, d, b1, b2, c, f, A, rhs )

%*****************************************************************************80
%
%% INTERIOR sets up the matrix and right hand side at interior nodes.
%
%  Discussion:
%
%    Nodes are assigned a single index K, which increases as:
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
%    If we number rows from bottom I = 1 to top I = NY
%    and columns from left J = 1 to right J = NX, then the relationship
%    between the single index K and the row and column indices I and J is:
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
%    27 August 2013
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
%    Input, function pointer @B1(X,Y), evaluates the convective field in x-direction.
%
%    Input, function pointer @B2(X,Y), evaluates the convective field in y-direction.
%
%    Input, function pointer @C(X,Y), evaluates the reactive coefficient.
%
%    Input, function pointer @F(X,Y), evaluates the heat source term.
%
%    Input, real sparse A(N,N), the system matrix, without any entries set.
%
%    Input, real RHS(N), the system right hand side, without any entries set.
%
%    Output, real sparse A(N,N), the system matrix, with the entries for the
%    interior nodes filled in.
%
%    Output, real RHS(N), the system right hand side, with the entries for the
%    interior nodes filled in.
%

%
%  For now, assume X and Y are equally spaced.
%
  dx = x(2) - x(1);
  dy = y(2) - y(1);

  for ic = 2 : ny - 1
    for jc = 2 : nx - 1

      in = ic + 1;
      is = ic - 1;
      je = jc + 1;
      jw = jc - 1;

      kc = ( ic - 1 ) * nx + jc;
      ke = kc + 1;
      kw = kc - 1;
      kn = kc + nx;
      ks = kc - nx;

      dce = d ( 0.5 * ( x(jc) + x(je) ),         y(ic) );
      dcw = d ( 0.5 * ( x(jc) + x(jw) ),         y(ic) );
      dcn = d (         x(jc),           0.5 * ( y(ic) + y(in) ) );
      dcs = d (         x(jc),           0.5 * ( y(ic) + y(is) ) );

      A(kc,kc) = ( dce + dcw ) / dx / dx + ( dcn + dcs ) / dy / dy;
      A(kc,ke) = - dce         / dx / dx;
      A(kc,kw) =       - dcw   / dx / dx;
      A(kc,kn) =                           - dcn         / dy / dy;
      A(kc,ks) =                                 - dcs   / dy / dy;
      
      
      b1c = b1( x(jc),                            y(ic) ) ;
      b2c = b2( x(jc),                            y(ic) ) ;
      
      if( b1c>0 )         
              
          A(kc, kc) = A(kc, kc) + b1c / dx;
          A(kc, kw) = A(kc, kw) - b1c / dx;
          
      else
          
          A(kc, kc) = A(kc, kc) - b1c / dx;
          A(kc, ke) = A(kc, ke) + b1c / dx;
   
      end
      
      if( b2c>0 )         
              
          A(kc, kc) = A(kc, kc) + b1c / dx;
          A(kc, ks) = A(kc, ks) - b1c / dx;
          
      else
          
          A(kc, kc) = A(kc, kc) - b1c / dx;
          A(kc, kn) = A(kc, kn) + b1c / dx;
   
      end    
      
      
      A(kc,kc) = A(kc,kc) + c( x(jc), y(ic) );

      rhs(kc,1) = f ( x(jc), y(ic) );

    end
  end

  return
end