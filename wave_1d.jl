# 1D wave propagation modeling with Spectral Element Method (SEM)
# This code can be used for studying the basic conception of SEM
# @author: Peng Yong, China  (yongpeng2013@gmail.com)
# Date: 20/08/2019     
include("define_derivative_matrix.jl")
include("source.jl")
using GR

##############################################
# number of spectral elements
  NSPEC = 500
# number of GLL points (polynomial degree plus one)
  NGLL = 4
# Courant-Friedrichs-Lewy (CFL) stability value
  CFL = 0.4
# number of timesteps
  nt = 3000 #1641
# fixed boundary conditions
  fix_bc = false
# plane wave source
  PlaneWave = false
# model parameters
  Length   = 50.0e+03  # m
  Density  = 1.5e+03   # kg/m^3
  Rigidity = 3.0e+10   # Pa
# source and receivers
  ispec_source = 1 #Int64(NSPEC/2)
  i_source = 1
  source_amp = 1.0e+7
##############################################
# number of global points
  NGLOB = (NGLL-1) * NSPEC + 1
# parameters for grids
  v  = sqrt(Rigidity/Density)
  dh =  Length/(NGLOB-1)
  dt = CFL*dh/v
  hdur = 82.035395797065604*dt
# define GLL point (xigll), weights (wgll),
# and polynomial derivatives (hprime)
  (xigll, wgll, hprime) = define_derivative_matrix(NGLL)
# evenly spaced anchors between 0 and 1
  x1 = zeros(NSPEC)
  x2 = zeros(NSPEC)
  for ispec = 1:NSPEC
    x1[ispec] = Length*(ispec-1)/NSPEC
    x2[ispec] = Length*(ispec  )/NSPEC
  end
# Set up the mesh properties
  rho      = zeros(NGLL, NSPEC)
  mu       = zeros(NGLL, NSPEC)
  dxi_dx   = zeros(NGLL, NSPEC)
  jacobian = zeros(NGLL, NSPEC)
  for ispec = 1:NSPEC
    for i = 1:NGLL
      rho[i,ispec] = Density
      mu[i,ispec] = Rigidity
      dxi_dx[i,ispec] = 2.0/(x2[ispec]-x1[ispec]) # This is d(xi) / dx
      jacobian[i,ispec] = (x2[ispec]-x1[ispec])/2.0
    end
  end

# Set up local to global numbering
  ibool = zeros(Int, NGLL, NSPEC)
  iglob = 1
  for ispec = 1:NSPEC
    for i = 1:NGLL
      if i > 1
        global iglob = iglob + 1
      end
      ibool[i, ispec] = iglob
    end
  end
  # println(iglob)
# Compute the position of the global grid points
  x = zeros(NGLOB)
  for ispec = 1:NSPEC
    for i = 1:NGLL
        iglob = ibool[i,ispec]
        x[iglob] = 0.5*(1.0-xigll[i])*x1[ispec]+0.5*(1.0+xigll[i])*x2[ispec]
    end
  end

# calculate the global mass matrix 'mass_global'
  mass_global = zeros(NGLOB)
  for ispec = 1:NSPEC
    for i = 1:NGLL
      mass_local = wgll[i]*(rho[i,ispec])*jacobian[i,ispec]
      iglob = ibool[i,ispec]
      mass_global[iglob] += mass_local
    end
  end


  deltat = dt
  deltatover2 = dt/2.0
  deltatsqover2 = dt*dt/2.0
  iglob_source = ibool[i_source,ispec_source]
# Initialize the fields to zero
  displ = zeros(NGLOB)
  veloc = zeros(NGLOB)
  accel = zeros(NGLOB)
  if PlaneWave
      for iglob = 1:NGLOB
            displ[iglob] = sin(pi*x[iglob]/Length)
      end
  end
  temp = zeros(NGLL)
# Main time loop
  for it = 1:nt
      # println(sum(abs.(displ)))
      for iglob = 1:NGLOB
          displ[iglob] += deltat*veloc[iglob] + deltatsqover2*accel[iglob]
          veloc[iglob] += deltatover2*accel[iglob]
          accel[iglob]  = 0.0
      end
      # println(sum(abs.(displ)))
      for ispec = 1:NSPEC
          for k = 1:NGLL
            #  Compute d(u) / d(xi)
            du_dxi = 0.0
            for i = 1:NGLL
                iglob = ibool[i,ispec]
                du_dxi = du_dxi + displ[iglob]*hprime[k,i]
            end
            epsilon = du_dxi*dxi_dx[k,ispec]
            sigma = mu[k,ispec]*epsilon
            temp[k] = jacobian[k,ispec]*sigma*dxi_dx[k,ispec]
          end # first loop over the GLL points
          for j = 1:NGLL
              templ = 0.0
              for k = 1:NGLL
                templ += temp[k]*hprime[k,j]*wgll[k]
              end
              iglob = ibool[j,ispec]
              accel[iglob] -= templ
          end # Second loop over the GLL points
      end # End loop over all spectral elements
    # println(accel[500])

    # add source
    if (!PlaneWave)
        accel[iglob_source] += source_amp*source_time((it-1)*dt-hdur,hdur)
    end

    # boundary conditions
    if fix_bc
        accel[1] = 0.0
        accel[NGLOB] = 0.0
    end

    for iglob = 1:NGLOB
        # updates acceleration
        accel[iglob] /=  mass_global[iglob]
        # Corrector update velocity
        veloc[iglob] +=  deltatover2*accel[iglob]
    end

  end
plot(x,displ)
