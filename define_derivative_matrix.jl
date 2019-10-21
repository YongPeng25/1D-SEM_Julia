# Gauss Lobatto Legendre points and integration weights,
# which can be found in Bernhard Schuberth's (2003) thesis:
# The Spectral Element Method for Seismic Wave Propagation
# Theory, Implementation and Comparison to Finite Difference Methods
#
# The functions for the derivative of the Lagrange polynomials at the GLL points
# are modified from the Fortran code in Specfem1D pakage (Dimitri Komatitsch and Jeroen Tromp)
#
function define_derivative_matrix(N)
    # GLL_Points in [-1, 1]
    if N==3
        xigll = [-1, 0, 1]
    elseif N==4
        xigll = [-1, -0.4472135954999579, 0.4472135954999579, 1]
    elseif N==5
        xigll = [-1, -0.6546536707, 0, 0.6546536707, 1]
    elseif N==6
        xigll = [-1, -0.7650553239, -0.2852315164, 0.2852315164, 0.7650553239,
                     1]
    elseif N==7
        xigll = [-1, -0.8302238962, -0.4688487934, 0, 0.4688487934,
                     0.8302238962, 1]
    elseif N==8
        xigll = [-1, -0.8717401485, -0.5917001814, -0.2092992179, 0.2092992179,
                     0.5917001814, 0.8717401485, 1]
    else
        println("N is supposed to meet the condition: 3<= N <=8. ")
        return
    end
    # Weights on GLL_Points
    if N==3
        wgll = [0.3333333333333333, 1.3333333333333333, 0.3333333333333333]
    elseif N==4
        wgll = [0.1666666666666667, 0.8333333333333334, 0.8333333333333334, 0.1666666666666667]
    elseif N==5
        wgll = [0.1, 0.5444444444, 0.7111111111, 0.5444444444, 0.1]
    elseif N==6
        wgll = [0.0666666666, 0.3784749562, 0.5548583770, 0.5548583770,
                     0.3784749562, 0.0666666666]
    elseif N==7
        wgll = [0.0476190476, 0.2768260473, 0.4317453812, 0.4876190476,
                     0.4317453812, 0.2768260473, 0.0476190476]
    elseif N==8
        wgll = [0.0357142857, 0.2107042271, 0.3411226924, 0.4124587946,
                     0.4124587946, 0.3411226924, 0.2107042271, 0.0357142857]
    else
        println("N is supposed to meet the condition: 3<= N <=8. ")
        return
    end
# calculate derivatives of the Lagrange polynomials
# and precalculate some products in double precision
# hprime(i,j) = h'_i(xigll_j) by definition of the derivative matrix
  hprime = zeros(N,N);
  for i=1:N
    for j=1:N
      hprime[j, i] = lagrange_derivative_gll(i,j,xigll)
    end
  end
#  println(hprime)
  xigll, wgll, hprime
end

"""
Calculates the values of the derivative of the Lagrange polynomials
at the GLL points
"""
function lagrange_derivative_gll(i,j,xigll)
#------------------------------------------------------------------------
#
#     Compute the value of the derivative of the I-th
#     Lagrange interpolant through the
#     N-1 Gauss-Lobatto Legendre points xiGLL at point xiGLL(j)
#
#------------------------------------------------------------------------
  N = length(xigll)
  degpoly = N - 1
  if (i == 1 && j == 1)
    deriv_GLL = - (degpoly)*((degpoly)+1.0)/ 4.0
  elseif (i == N && j == N)
    deriv_GLL = (degpoly)*((degpoly)+1.0)/ 4.0
  elseif (i == j)
    deriv_GLL = 0.0
  else
    deriv_GLL = pnleg(xigll[j],degpoly)/(pnleg(xigll[i],degpoly)*(xigll[j]-xigll[i]))
              + (1.0-xigll[j]*xigll[j])*pndleg(xigll[j],degpoly)/((degpoly)*((degpoly)+1.0)*
              pnleg(xigll[i],degpoly)*(xigll[j]-xigll[i])*(xigll[j]-xigll[i]))
  end
  deriv_GLL
end

function pnleg(Z, N)
#------------------------------------------------------------------------
#
#       Compute the value of the Nth order Legendre polynomial at Z.
#       Based on the recursion formula for the Legendre polynomials.
#
#------------------------------------------------------------------------
  P1  = 1.0
  P2  = Z
  P3  = P2

  for K = 1:N-1
    FK  = K
    P3  = ((2.0*FK+1.0)*Z*P2 - FK*P1)/(FK+1.0)
    P1  = P2
    P2  = P3
  end

  P3
end

function pndleg(Z,N)
#------------------------------------------------------------------------
#
#     Compute the derivative of the Nth order Legendre polynomial at Z.
#     Based on the recursion formula for the Legendre polynomials.
#
#------------------------------------------------------------------------
  P1   = 1.0
  P2   = Z
  P1D  = 0.0
  P2D  = 1.0
  P3D  = 1.0

  for K = 1:N-1
    FK  = K
    P3  = ((2.0*FK+1.0)*Z*P2 - FK*P1)/(FK+1.0)
    P3D = ((2.0*FK+1.0)*P2 + (2.0*FK+1.0)*Z*P2D - FK*P1D) / (FK+1.0)
    P1  = P2
    P2  = P3
    P1D = P2D
    P2D = P3D
  end

  P3D
end
