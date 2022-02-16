MBC = function(data, framerate, alpha = 1, nopt = 15, t.thres = seq(0.03, 0.3, length.out = 50)/100, plot = T) {
  # Data is a data-frame with columns (X, Y, time, localisation uncertainty (standard deviation))
  
  if(plot) {old.pars = par(no.readonly=TRUE); on.exit(par(old.pars))}
  if(plot) {par(mfrow = c(4,4), mai = c(0.6,0.6,0.4,0.2)); layout(matrix(c(1,1,3,3,1,1,3,3,1,1,3,3,2,2,2,2), 4, 4, T))}
  
  X = ppp(x = data[,1], y = data[,2], marks = as.integer(data[,3]), window = owin(range(data[,1]), range(data[,2])))
  unc = data[,4]
  if(plot) {plot(data[,1:2], pch = 19, asp = 1, cex = 0.5, main = "Raw data")}
  
  breakl(); catprint("Estimating rate parameters.")
  rate.estimates = estimate.model.auto(X, framerate, unc, alpha = alpha, cut = NA, nopt = nopt)
  
  ord.params = as.numeric(rate.estimates$par[c(1,3,4,2)])
  names(ord.params) = c("r_F", "r_D", "r_R", "r_B")
  print(signif(ord.params, 2))
  
  breakl(); catprint("Finding most likely clustering.")
  
  eint = (1-alpha)*intensity(X)
  declus = STBHClus_optimize(data, framerate, ord.params, eint, t.thres, plt = plot)
  od = c(declus$loglik, declus$S); names(od) = c("loglik", "Time-dilation")
  print(signif(od, 4))
  
  MBC.data = declus$merged %>% select(-time, -n)
  MBC.data = rename(MBC.data, uncertainty = sd)
  
  if(plot) {plot(MBC.data[,1:2], pch = 19, asp = 1, cex = 0.5, main = "MBC-corrected")}
  
  catprint("Done!")
  catprint("")
  return(MBC.data)
}
