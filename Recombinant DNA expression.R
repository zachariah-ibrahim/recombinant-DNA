# RECOMBINANT DNA EXPRESSION #

library(deSolve)

# INITIAL CONCENTRATIONS #
state.variables <- c(g0 = 5E-10, p0 = 5E-10, xs = 0, x0 = 0, m0 = 0,
                     r0 = 5E-10, ys = 0, y0 = 0, z0 = 0, zs = 0)

# EFFECTIVE PARAMETERS #

# typical cell volume = 3E-12L #
# reaction chamber volume = 1E-6L #
parm.values <- c(k1 = 6E9, k2 = 600, k3 = 60, k4 = 0.9, km = 18,
                 k5 = 6E9, k6 = 135, k7 = 30, k8 = 0.9,
                 kt = 60)

# CONCENTRATION FUNCTION #
con.profile <- function(t, state.variables, parm.values){
  with(as.list(c(state.variables, parm.values)),{
    
    dg0 <- k2*xs + k4*x0 - k1*g0*p0
    dp0 <- k2*xs + k4*x0 - k1*g0*p0
    dxs <- k1*g0*p0 - k2*xs - k3*xs
    dx0 <- k3*xs - k4*x0
    dm0 <- k4*x0 + k6*ys + k8*y0 - km*m0 - k5*m0*r0
    dr0 <- k6*ys + k8*y0 - k5*m0*r0
    dys <- k5*m0*r0 - k6*ys - k7*ys
    dy0 <- k7*ys - k8*y0
    dz0 <- k8*y0 - kt*z0
    dzs <- kt*z0
    
    return(list(c(dg0, dp0, dxs, dx0, dm0,
                  dr0, dys, dy0, dz0, dzs)))
  })
}

# OUTPUT #
out <- ode(y = state.variables,
           parms = parm.values,
           times = seq(0, 60, 0.001),
           func = con.profile)

#plot(out)
