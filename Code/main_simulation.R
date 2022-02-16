Add_blinking = function(data, parameters, framerate = 25, stds = rgamma(1000,4,4/20), stop_if_unusual = TRUE) {
  # Data is a data-frame with columns (X, Y)
  # Parameters are the blinking parameters (r_F, r_D, r_R, r_B)
  # Framerate is the number of frames per second recorded by the microscope
  # stds is a vector of localisation uncertainties. The uncertainties of simulated blinks are drawn with replacement from stds.
  
  orate = parameters[1]
  bpar = parameters[c(4,2,3)]
  
  if(stop_if_unusual) {
    nb = round(moment.approx(bpar, 1/framerate)[1])
    mean_good = between(nb, 2, 25)
    
    sl = meanapprox(bpar, 1/framerate)*orate
    start_good = (sl <= 0.1)
    
    mstd = mean(stds)
    std_good = mstd < 50
    
    if(!mean_good){
      warn_notice = paste('Parameters would result in an unusual (~', nb, ') number of reappearances. Set stop_if_unusual = F to ignore.', sep = "")
      stop(warn_notice)
    }
    if(!start_good) {
      warn_notice = paste('Parameters would result in very long lived blinking compared to the activation time (ratio ~', round(sl,3), '). Set stop_if_unusual = F to ignore.', sep = "")
      stop(warn_notice)
    }
    if(!std_good) {
      warn_notice = paste('Mean localisation uncertainty is unusually large (', round(mstd,1), '). Set stop_if_unusual = F to ignore.', sep = "")
      stop(warn_notice)
    }
  }
  
  fr = framerate
  gen = pars.to.simulator(bpar, 0, fr, T)
  pargen = function() rexp(1, orate)
  
  W = owin(range(data[,1]), range(data[,2]))
  dpp = ppp(x = data[,1], y = data[,2], window = W)
  
  dtpar = add.cluster.layer(dpp, gen, 0, pargen, fr); dtpar = noise_add(dtpar, stds)$D; dtpar$marks = as.integer(dtpar$marks)
  df = dtpar %$% data.frame(x = x, y = y, time = marks, sds = attributes(.)$s)
  
  return(df)
}
