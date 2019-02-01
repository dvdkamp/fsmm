require(deSolve)

run.fsmm = function(fsmm.model.input  = read.csv("fsmm.example.input.csv"),
                    fsmm.model.pars,
                    ## vector of model Pars, Order: A,B,K_s,M_MAX,DR_O
                    RHO_S = 400,
                    # kg m-3 : Stick density  (taken from Nelson, 2000))
                    r = 0.0065,
                    # radius of stick (m)
                    L = 0.41,
                    # length of stick (m)
                    m_i = 0.1,
                    CALC_SHAPE_FACTOR = T) {
  # initial moisture level
  
  ######################################################################
  ## Load in fortran code
  
  system('R CMD SHLIB fsmm.f')
  dyn.load('fsmm.so')
  
  
  ######################################################################
  ## calculate area of shadow cast by stick
  ## from Monteith and Unsworth (2008) pg 105
  
  if (CALC_SHAPE_FACTOR) {
    fsmm.model.input$shape.factor <-
      ifelse(
        fsmm.model.input$Alt_s > 0,
        2 * r * L * (1 / sin(fsmm.model.input$Alt_s)) * sqrt(
          1 - cos(fsmm.model.input$Alt_s) ^ 2 * cos(fsmm.model.input$Azi_s) ^ 2
        ) + pi * r ^ 2 * (1 / tan(fsmm.model.input$Alt_s)) * cos(fsmm.model.input$Azi_s),
        0
      )
    
  } else {
    fsmm.model.input$shape.factor <- 1
    
  }
  
  ######################################################################
  ## Create Time variable (in seconds)
  
  time.out <- seq(1, nrow(fsmm.model.input) * 3600, 3600)
  
  ######################################################################
  ## Create Forcing Variables matrix
  
  forcing.variables <- list(
    matrix(
      ncol = 2,
      data = c(time.out, fsmm.model.input$k.down.dir)
    ),
    matrix(
      ncol = 2,
      data = c(time.out, fsmm.model.input$k.down.diff)
    ),
    matrix(
      ncol = 2,
      data = c(time.out, fsmm.model.input$shape.factor)
    ),
    matrix(
      ncol = 2,
      data = c(time.out, fsmm.model.input$temp)
    ),
    matrix(
      ncol = 2,
      data = c(time.out, fsmm.model.input$L_dn)
    ),
    matrix(
      ncol = 2,
      data = c(time.out, fsmm.model.input$rh)
    ),
    matrix(
      ncol = 2,
      data = c(time.out, fsmm.model.input$wind.speed)
    ),
    matrix(
      ncol = 2,
      data = c(time.out, fsmm.model.input$precip)
    )
  )
  
  ######################################################################
  ## Create state vector with initial values
  
  state <-
    c(
      fuel.temp.mod.core = head(fsmm.model.input$temp, 1),
      fuel.mois.mod.outer = m_i * (pi * r ^ 2 * L * RHO_S),
      fuel.mois.mod.core = m_i * (pi * r ^ 2 * L * RHO_S),
      fuel.temp.mod.core = head(fsmm.model.input$temp, 1)
    )
  
  ######################################################################
  ## Create parameter vector (this includes the radius and length of the stick)
  
  fsmm.model.pars.full <- c(fsmm.model.pars, r, L)
  
  ######################################################################
  ## Run Model. the output 
  
  model.out <- ode(
    y = state,
    times = time.out,
    func = 'fmccore',
    parms = fsmm.model.pars.full,
    verbose = FALSE,
    ynames = TRUE,
    dllname = 'fmcdv',
    initfunc = 'fmcpar',
    nout = 2,
    outnames = c("fuel.mois.mod", "fuel.temp.mod"),
    forcings = forcing.variables,
    initforc = 'fmcforc'
  )
  
  
}
