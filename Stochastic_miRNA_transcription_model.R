library(GillespieSSA2)
library(ggplot2)
library(gridExtra)
#wrap this into a function on the weekend so you can do lots of different plots to compare ------
Gillespie_model <- function(parms){
  # Initial values
  r_0 <- parms[["k_R"]] / parms[["gamma_R"]] #r = mRNA free of miRNA - eYFP signal
  r_star_0 <- 0
  #miRNA_0 <- 0
  x0  <- c(r = r_0, r_star = r_star_0,miRNA = parms[["miRNA_0"]])
  #reactions
  reactions <- list(
    #        propensity function     effects                                          name for reaction
    reaction(~k_R,                  c(r = +1),                                      "transcription_of_mRNA"),
    reaction(~k_on * miRNA * r, c(r = -1, miRNA = -1, r_star = +1),                     "formation_of_complex"),
    reaction(~k_off * r_star,        c(r = +1, r_star = -1, miRNA = +1),            "release_of_complex"),
    reaction(~gamma_R * r,          c(r = -1),                                     "natural_breakdown"),
    reaction(~gamma_R_star * r_star, c(r_star = -1, miRNA = +1),                    "breakdown_of_complex")
  )
  tmax <- 500
  
  out_temp <- GillespieSSA2::ssa(initial_state = x0,reactions = reactions,params = parms,final_time = tmax,census_interval = 1)
  data.df <- data.frame("time" = round(out_temp$time,0), out_temp$state)
  data.df <- data.df[200:400,]
  return(data.df)
}
#sensitivity function ------
Sensitivity_func <- function(parms,parameter,values){
  i <- 1
  kr.noise.df <- data.frame()
  for(i in 1:length(values)){
    #loop through all the values for said parameter
    parms[[parameter]] <- values[i]
    kr.noise.df[i,1] <- parms[[parameter]]
    #Run G model using that value and calculate noise
    out <- Gillespie_model(parms = parms)
    SD_r <- sd(x = out$r)
    Mean_r <- mean(x = out$r)
    if(is.na(SD_r)){
      SD_r <- 0
    }
    Noise_r <- SD_r / Mean_r
    if(is.na(Noise_r) | is.nan(Noise_r)){
      Noise_r <- 0
    }
    #Calculating Threshold value for k_r
    lambda <- (parms[["gamma_R_star"]] + parms[["k_off"]])/parms[["k_on"]]
    theta <- (parms[["gamma_R_star"]]/parms[["gamma_R"]]) * parms[["miRNA_0"]]
    k_R.threshold <- theta * parms[["gamma_R"]]
    #only the noise for mRNA
    kr.noise.df[i,2] <- Noise_r
    kr.noise.df[i,3] <- SD_r
    kr.noise.df[i,4] <- Mean_r
    print(i)
  }
  names(kr.noise.df) <- c(paste0(parameter), "Noise","SD","Mean")
  output.list <- list(kr.noise.df,k_R.threshold)
  return(output.list)
}
#k_r plotting function ----
kr_plot <- function(dataset, miRNA_level){
  kr.noise.plot <- ggplot() +
    geom_line(data = dataset[[1]], aes(x = k_R, y = Noise),colour = "blue", linewidth = lnwd) +
    geom_vline(xintercept = dataset[[2]], linetype = "dashed", colour = "blue") +
    xlab("Transcription rate (k_R)") +
    ylab("mRNA noise") +
    ggtitle(paste("mRNA noise when\nmiRNA =",miRNA_level))
  
  kr.mean.plot <- ggplot() +
    geom_line(data = dataset[[1]], aes(x = k_R, y = Mean), colour = "red", linewidth = lnwd) +
    geom_vline(xintercept = dataset[[2]], linetype = "dashed", colour = "red") +
    xlab("Transcription rate (k_R)") +
    ylab("mRNA mean expression") +
    ggtitle(paste("mRNA mean expression when\nmiRNA =",miRNA_level))
  
  grid.arrange(grobs = list(kr.noise.plot,kr.mean.plot), ncol = 2, nrow = 1)
  
}
lnwd <- 1 #width of line

#running variable numbers of miRNA -----
variable_miRNA <- function(parms,levels, parameter, values){
  z <- 1
  return_list <- list()
  for(z in 1:length(levels)){
    parms[["miRNA_0"]] <- levels[z]
    return_list[[z]] <- Sensitivity_func(parms = parms, parameter = parameter, values = values)
    print(paste0("mRNA series:", z))
  }
  return(return_list)
}
z <- 1
#Testing increasing miRNA -----
#i set k_R to 3 which means that it is able to suceed a theta threshold of 300
#equivalent to 10 miRNA with current parameter values
parms10 <- c(k_R = 1, #transcription rate
             k_on = 2.5, #rate at which free miRNAs are removed from the system 
             k_off = 0.05, #unbinding of the miRNA from its target
             gamma_R = 0.5, #rate of intrinsic decay 
             gamma_R_star = 0.3, #destruction of miRNAs target
             miRNA_0 = 10
)

miRNA_seq <- c(0,25,50,100,150,250)
values.seq <- seq(0,100, by = 1)
kr.list <- variable_miRNA(parms = parms10, levels = miRNA_seq, parameter = "k_R", values = values.seq)

kr.noise.plot <- ggplot() +
  geom_line(data = kr.list[[1]][[1]], aes(x = k_R, y = Noise,colour = "0"), linewidth = lnwd) +
  geom_vline(aes(xintercept = kr.list[[1]][[2]],colour = "0"), linetype = "dashed") +
  geom_line(data = kr.list[[2]][[1]], aes(x = k_R, y = Noise,colour = "25"), linewidth = lnwd) +
  geom_vline(aes(xintercept = kr.list[[2]][[2]], colour = "25"), linetype = "dashed") +
  geom_line(data = kr.list[[3]][[1]], aes(x = k_R, y = Noise,colour = "50"), linewidth = lnwd) +
  geom_vline(aes(xintercept = kr.list[[3]][[2]], colour = "50"), linetype = "dashed") +
  geom_line(data = kr.list[[4]][[1]], aes(x = k_R, y = Noise,colour = "100"), linewidth = lnwd) +
  geom_vline(aes(xintercept = kr.list[[4]][[2]], colour = "100"), linetype = "dashed") +
  geom_line(data = kr.list[[5]][[1]], aes(x = k_R, y = Noise,colour = "150"), linewidth = lnwd) +
  geom_vline(aes(xintercept = kr.list[[5]][[2]], colour = "150"), linetype = "dashed") +
  geom_line(data = kr.list[[6]][[1]], aes(x = k_R, y = Noise,colour = "250"), linewidth = lnwd) +
  geom_vline(aes(xintercept = kr.list[[6]][[2]], colour = "250"), linetype = "dashed") +
  xlab("Transcription rate (k_R)") +
  xlim(0,100) +
  ylab("Gene noise") +
  scale_colour_manual(
    values = c("0" = "blue", "25" = "red", "50" = "darkgreen", 
               "100" = "purple", "150" = "orange", "250" = "brown"),
    breaks = miRNA_seq  # Order in the legend
  ) +
  labs(colour = "Total miRNA") +
  ggtitle(paste("The effect of increasing miRNA\non mRNA mean expression"))

kr.mean.plot <- ggplot() +
  geom_line(data = kr.list[[1]][[1]], aes(x = k_R, y = Mean,colour = "0"), linewidth = lnwd) +
  geom_vline(aes(xintercept = kr.list[[1]][[2]],colour = "0"), linetype = "dashed") +
  geom_line(data = kr.list[[2]][[1]], aes(x = k_R, y = Mean,colour = "25"), linewidth = lnwd) +
  geom_vline(aes(xintercept = kr.list[[2]][[2]], colour = "25"), linetype = "dashed") +
  geom_line(data = kr.list[[3]][[1]], aes(x = k_R, y = Mean,colour = "50"), linewidth = lnwd) +
  geom_vline(aes(xintercept = kr.list[[3]][[2]], colour = "50"), linetype = "dashed") +
  geom_line(data = kr.list[[4]][[1]], aes(x = k_R, y = Mean,colour = "100"), linewidth = lnwd) +
  geom_vline(aes(xintercept = kr.list[[4]][[2]], colour = "100"), linetype = "dashed") +
  geom_line(data = kr.list[[5]][[1]], aes(x = k_R, y = Mean,colour = "150"), linewidth = lnwd) +
  geom_vline(aes(xintercept = kr.list[[5]][[2]], colour = "150"), linetype = "dashed") +
  geom_line(data = kr.list[[6]][[1]], aes(x = k_R, y = Mean,colour = "250"), linewidth = lnwd) +
  geom_vline(aes(xintercept = kr.list[[6]][[2]], colour = "250"), linetype = "dashed") +
  xlab("Transcription rate (k_R)") +
  xlim(0,100) +
  ylab("mRNA mean expression") +
  scale_colour_manual(
    values = c("0" = "blue", "25" = "red", "50" = "darkgreen", 
               "100" = "purple", "150" = "orange", "250" = "brown"),
    breaks = miRNA_seq  # Order in the legend
  ) +
  labs(colour = "Total miRNA") +
  ggtitle(paste("The effect of increasing miRNA\non mRNA mean expression"))

grid.arrange(grobs = list(kr.noise.plot,kr.mean.plot), ncol = 2, nrow = 1)
#Effect of natural degradation on expression/ noise of a lowly expressed gene-----
parms.low <- c(k_R = 1, #transcription rate
             k_on = 2.5, #rate at which free miRNAs are removed from the system 
             k_off = 0.05, #unbinding of the miRNA from its target
             gamma_R = 0.01, #rate of intrinsic decay 
             gamma_R_star = 0.3, #destruction of miRNAs target
             miRNA_0 = 10
)

low.gamma_R.list <- variable_miRNA(parms = parms.low, levels = miRNA_seq, parameter = "k_R", values = gamma_R.seq)

low.gamma_R.Mean.plot <- ggplot() +
  geom_line(data = low.gamma_R.list[[1]][[1]], aes(x = gamma_R, y = Mean,colour = "0"), linewidth = lnwd) +
  geom_vline(aes(xintercept = low.gamma_R.list[[1]][[2]],colour = "0"), linetype = "dashed") +
  geom_line(data = low.gamma_R.list[[2]][[1]], aes(x = gamma_R, y = Mean,colour = "10"), linewidth = lnwd) +
  geom_vline(aes(xintercept = low.gamma_R.list[[2]][[2]], colour = "10"), linetype = "dashed") +
  geom_line(data = low.gamma_R.list[[3]][[1]], aes(x = gamma_R, y = Mean,colour = "20"), linewidth = lnwd) +
  geom_vline(aes(xintercept = low.gamma_R.list[[3]][[2]], colour = "20"), linetype = "dashed") +
  geom_line(data = low.gamma_R.list[[4]][[1]], aes(x = gamma_R, y = Mean,colour = "30"), linewidth = lnwd) +
  geom_vline(aes(xintercept = low.gamma_R.list[[4]][[2]], colour = "30"), linetype = "dashed") +
  geom_line(data = low.gamma_R.list[[5]][[1]], aes(x = gamma_R, y = Mean,colour = "40"), linewidth = lnwd) +
  geom_vline(aes(xintercept = low.gamma_R.list[[5]][[2]], colour = "40"), linetype = "dashed") +
  geom_line(data = low.gamma_R.list[[6]][[1]], aes(x = gamma_R, y = Mean,colour = "50"), linewidth = lnwd) +
  geom_vline(aes(xintercept = low.gamma_R.list[[6]][[2]], colour = "50"), linetype = "dashed") +
  xlab("Transcription Rate(gamma_R)") +
  ylab("mRNA mean expression") +
  labs(colour = "miRNA amount") +
  ggtitle(paste("mRNA expression with variable miRNA"))

low.gamma_R.Noise.plot <- ggplot() +
  geom_line(data = low.gamma_R.list[[1]][[1]], aes(x = gamma_R, y = Noise,colour = "0"), linewidth = lnwd) +
  geom_vline(aes(xintercept = low.gamma_R.list[[1]][[2]],colour = "0"), linetype = "dashed") +
  geom_line(data = low.gamma_R.list[[2]][[1]], aes(x = gamma_R, y = Noise,colour = "10"), linewidth = lnwd) +
  geom_vline(aes(xintercept = low.gamma_R.list[[2]][[2]], colour = "10"), linetype = "dashed") +
  geom_line(data = low.gamma_R.list[[3]][[1]], aes(x = gamma_R, y = Noise,colour = "20"), linewidth = lnwd) +
  geom_vline(aes(xintercept = low.gamma_R.list[[3]][[2]], colour = "20"), linetype = "dashed") +
  geom_line(data = low.gamma_R.list[[4]][[1]], aes(x = gamma_R, y = Noise,colour = "30"), linewidth = lnwd) +
  geom_vline(aes(xintercept = low.gamma_R.list[[4]][[2]], colour = "30"), linetype = "dashed") +
  geom_line(data = low.gamma_R.list[[5]][[1]], aes(x = gamma_R, y = Noise,colour = "40"), linewidth = lnwd) +
  geom_vline(aes(xintercept = low.gamma_R.list[[5]][[2]], colour = "40"), linetype = "dashed") +
  geom_line(data = low.gamma_R.list[[6]][[1]], aes(x = gamma_R, y = Noise,colour = "50"), linewidth = lnwd) +
  geom_vline(aes(xintercept = low.gamma_R.list[[6]][[2]], colour = "50"), linetype = "dashed") +
  xlab("Natural degradation (gamma_R)") +
  ylab("mRNA noise") +
  labs(colour = "miRNA amount") +
  ggtitle(paste("mRNA noise with variable miRNA"))

grid.arrange(grobs = list(low.gamma_R.Noise.plot,low.gamma_R.Mean.plot), ncol = 2, nrow = 1)
#Effect of natural degradation on expression / noise on a highly expressed gene -----
parms.high <- c(gamma_R = 15, #transcription rate
               k_on = 2.5, #rate at which free miRNAs are removed from the system 
               k_off = 0.05, #unbinding of the miRNA from its target
               gamma_R = 0.01, #rate of intrinsic decay 
               gamma_R_star = 0.3, #destruction of miRNAs target
               miRNA_0 = 10
)
gamma_R.seq <- seq(0,0.3, by = 0.1)
high.gamma_R.list <- variable_miRNA(parms = parms.high, levels = miRNA_seq, parameter = "gamma_R", values = gamma_R.seq)

high.gamma_R.Mean.plot <- ggplot() +
  geom_line(data = high.gamma_R.list[[1]][[1]], aes(x = gamma_R, y = Mean,colour = "0"), linewidth = lnwd) +
  geom_vline(aes(xintercept = high.gamma_R.list[[1]][[2]],colour = "0"), linetype = "dashed") +
  geom_line(data = high.gamma_R.list[[2]][[1]], aes(x = gamma_R, y = Mean,colour = "10"), linewidth = lnwd) +
  geom_vline(aes(xintercept = high.gamma_R.list[[2]][[2]], colour = "10"), linetype = "dashed") +
  geom_line(data = high.gamma_R.list[[3]][[1]], aes(x = gamma_R, y = Mean,colour = "20"), linewidth = lnwd) +
  geom_vline(aes(xintercept = high.gamma_R.list[[3]][[2]], colour = "20"), linetype = "dashed") +
  geom_line(data = high.gamma_R.list[[4]][[1]], aes(x = gamma_R, y = Mean,colour = "30"), linewidth = lnwd) +
  geom_vline(aes(xintercept = high.gamma_R.list[[4]][[2]], colour = "30"), linetype = "dashed") +
  geom_line(data = high.gamma_R.list[[5]][[1]], aes(x = gamma_R, y = Mean,colour = "40"), linewidth = lnwd) +
  geom_vline(aes(xintercept = high.gamma_R.list[[5]][[2]], colour = "40"), linetype = "dashed") +
  geom_line(data = high.gamma_R.list[[6]][[1]], aes(x = gamma_R, y = Mean,colour = "50"), linewidth = lnwd) +
  geom_vline(aes(xintercept = high.gamma_R.list[[6]][[2]], colour = "50"), linetype = "dashed") +
  xlab("Natural degradation (gamma_R)") +
  ylab("mRNA mean expression") +
  labs(colour = "miRNA amount") +
  ggtitle(paste("mRNA expression with variable miRNA"))

high.gamma_R.Noise.plot <- ggplot() +
  geom_line(data = high.gamma_R.list[[1]][[1]], aes(x = gamma_R, y = Noise,colour = "0"), linewidth = lnwd) +
  geom_vline(aes(xintercept = high.gamma_R.list[[1]][[2]],colour = "0"), linetype = "dashed") +
  geom_line(data = high.gamma_R.list[[2]][[1]], aes(x = gamma_R, y = Noise,colour = "10"), linewidth = lnwd) +
  geom_vline(aes(xintercept = high.gamma_R.list[[2]][[2]], colour = "10"), linetype = "dashed") +
  geom_line(data = high.gamma_R.list[[3]][[1]], aes(x = gamma_R, y = Noise,colour = "20"), linewidth = lnwd) +
  geom_vline(aes(xintercept = high.gamma_R.list[[3]][[2]], colour = "20"), linetype = "dashed") +
  geom_line(data = high.gamma_R.list[[4]][[1]], aes(x = gamma_R, y = Noise,colour = "30"), linewidth = lnwd) +
  geom_vline(aes(xintercept = high.gamma_R.list[[4]][[2]], colour = "30"), linetype = "dashed") +
  geom_line(data = high.gamma_R.list[[5]][[1]], aes(x = gamma_R, y = Noise,colour = "40"), linewidth = lnwd) +
  geom_vline(aes(xintercept = high.gamma_R.list[[5]][[2]], colour = "40"), linetype = "dashed") +
  geom_line(data = high.gamma_R.list[[6]][[1]], aes(x = gamma_R, y = Noise,colour = "50"), linewidth = lnwd) +
  geom_vline(aes(xintercept = high.gamma_R.list[[6]][[2]], colour = "50"), linetype = "dashed") +
  xlab("Natural degradation (gamma_R)") +
  ylab("mRNA noise") +
  labs(colour = "miRNA amount") +
  ggtitle(paste("mRNA noise with variable miRNA"))


grid.arrange(grobs = list(high.gamma_R.noise.plot,high.gamma_R.mean.plot), ncol = 2, nrow = 1)

#Effect of Binding rate with Fixed miRNA ------
parms.k_on <- c(k_R = 1, #transcription rate
             k_on = 2.5, #rate at which free miRNAs are removed from the system 
             k_off = 0.05, #unbinding of the miRNA from its target
             gamma_R = 0.01, #rate of intrinsic decay 
             gamma_R_star = 0.3, #destruction of miRNAs target
             miRNA_0 = 250
)
k_on.seq <- c(0,1,2,5,10)
values.seq <- seq(1,100, by = 1)

x <- 1
k_on.list <- list()
for(x in 1:length(k_on.seq)){
  parms.k_on[["k_on"]] <- k_on.seq[x]
  k_on.list[[x]] <- Sensitivity_func(parms = parms.k_on, parameter = "k_R", values = values.seq)
  
}


k_on.noise.plot <- ggplot() +
  geom_line(data = k_on.list[[1]][[1]], aes(x = k_R, y = Noise,colour = "0"), linewidth = lnwd) +
  geom_vline(aes(xintercept = k_on.list[[1]][[2]],colour = "0"), linetype = "dashed") +
  geom_line(data = k_on.list[[2]][[1]], aes(x = k_R, y = Noise,colour = "1"), linewidth = lnwd) +
  geom_vline(aes(xintercept = k_on.list[[2]][[2]], colour = "1"), linetype = "dashed") +
  geom_line(data = k_on.list[[3]][[1]], aes(x = k_R, y = Noise,colour = "2"), linewidth = lnwd) +
  geom_vline(aes(xintercept = k_on.list[[3]][[2]], colour = "2"), linetype = "dashed") +
  geom_line(data = k_on.list[[4]][[1]], aes(x = k_R, y = Noise,colour = "5"), linewidth = lnwd) +
  geom_vline(aes(xintercept = k_on.list[[4]][[2]], colour = "5"), linetype = "dashed") +
  geom_line(data = k_on.list[[5]][[1]], aes(x = k_R, y = Noise,colour = "10"), linewidth = lnwd) +
  geom_vline(aes(xintercept = k_on.list[[5]][[2]], colour = "10"), linetype = "dashed") +
  xlab("Transcription rate (k_R)") +
  ylab("Gene noise") +
  scale_colour_manual(
    values = c("0" = "blue", "1" = "red", "2" = "darkgreen", 
               "5" = "purple", "10" = "orange"),
    breaks = k_on.seq  # Order in the legend
  ) +
  labs(colour = "rate of complex formation (k_on)") +
  ggtitle(paste("miRNA = 250"))
print(k_on.noise.plot)

k_on.mean.plot <- ggplot() +
  geom_line(data = k_on.list[[1]][[1]], aes(x = k_R, y = Mean,colour = "0"), linewidth = lnwd) +
  geom_vline(aes(xintercept = k_on.list[[1]][[2]],colour = "0"), linetype = "dashed") +
  geom_line(data = k_on.list[[2]][[1]], aes(x = k_R, y = Mean,colour = "1"), linewidth = lnwd) +
  geom_vline(aes(xintercept = k_on.list[[2]][[2]], colour = "1"), linetype = "dashed") +
  geom_line(data = k_on.list[[3]][[1]], aes(x = k_R, y = Mean,colour = "2"), linewidth = lnwd) +
  geom_vline(aes(xintercept = k_on.list[[3]][[2]], colour = "2"), linetype = "dashed") +
  geom_line(data = k_on.list[[4]][[1]], aes(x = k_R, y = Mean,colour = "5"), linewidth = lnwd) +
  geom_vline(aes(xintercept = k_on.list[[4]][[2]], colour = "5"), linetype = "dashed") +
  geom_line(data = k_on.list[[5]][[1]], aes(x = k_R, y = Mean,colour = "10"), linewidth = lnwd) +
  geom_vline(aes(xintercept = k_on.list[[5]][[2]], colour = "10"), linetype = "dashed") +
  xlab("Transcription rate (k_R)") +
  ylab("mRNA mean expression") +
  scale_colour_manual(
    values = c("0" = "blue", "1" = "red", "2" = "darkgreen", 
               "5" = "purple", "10" = "orange"),
    breaks = k_on.seq  # Order in the legend
  ) +
  labs(colour = "rate of complex formation (k_on)") +
  ggtitle(paste("miRNA = 250"))

grid.arrange(grobs = list(k_on.noise.plot,k_on.mean.plot), ncol = 2, nrow = 1)

#Effect of complex degradation on noise----
parms.gamma_R_star <- c(k_R = 1, #transcription rate
                k_on = 2.5, #rate at which free miRNAs are removed from the system 
                k_off = 0.05, #unbinding of the miRNA from its target
                gamma_R = 0.01, #rate of intrinsic decay 
                gamma_R_star = 0.3, #destruction of miRNAs target
                miRNA_0 = 250
)
gamma_R_star.seq <- seq(0,1, by = 0.2)
values.seq <- seq(40,270, by = 5)


x <- 1
gamma_R_star.list <- list()
for(x in 1:length(gamma_R_star.seq)){
  parms.gamma_R_star[["gamma_R_star"]] <- gamma_R_star.seq[x]
  print(parms.gamma_R_star[["gamma_R_star"]])
  gamma_R_star.list[[x]] <- Sensitivity_func(parms = parms.gamma_R_star, parameter = "k_R", values = values.seq)
}


gamma_R_star.noise.plot <- ggplot() +
  geom_line(data = gamma_R_star.list[[1]][[1]], aes(x = k_R, y = Noise,colour = "0"), linewidth = lnwd) +
  geom_vline(aes(xintercept = gamma_R_star.list[[1]][[2]],colour = "0"), linetype = "dashed") +
  geom_line(data = gamma_R_star.list[[2]][[1]], aes(x = k_R, y = Noise,colour = "0.2"), linewidth = lnwd) +
  geom_vline(aes(xintercept = gamma_R_star.list[[2]][[2]], colour = "0.2"), linetype = "dashed") +
  geom_line(data = gamma_R_star.list[[3]][[1]], aes(x = k_R, y = Noise,colour = "0.4"), linewidth = lnwd) +
  geom_vline(aes(xintercept = gamma_R_star.list[[3]][[2]], colour = "0.4"), linetype = "dashed") +
  geom_line(data = gamma_R_star.list[[4]][[1]], aes(x = k_R, y = Noise,colour = "0.6"), linewidth = lnwd) +
  geom_vline(aes(xintercept = gamma_R_star.list[[4]][[2]], colour = "0.6"), linetype = "dashed") +
  geom_line(data = gamma_R_star.list[[5]][[1]], aes(x = k_R, y = Noise,colour = "0.8"), linewidth = lnwd) +
  geom_vline(aes(xintercept = gamma_R_star.list[[5]][[2]], colour = "0.8"), linetype = "dashed") +
  geom_line(data = gamma_R_star.list[[6]][[1]], aes(x = k_R, y = Noise,colour = "1"), linewidth = lnwd) +
  geom_vline(aes(xintercept = gamma_R_star.list[[6]][[2]], colour = "1"), linetype = "dashed") +
  xlab("Transcription rate (k_R)") +
  ylab("Gene noise") +
  xlim(40,275) +
  labs(colour = "Rate of complex degradation\n(gamma_R_star)") +
  ggtitle(paste("miRNA = 250"))

gamma_R_star.mean.plot <- ggplot() +
  geom_line(data = gamma_R_star.list[[1]][[1]], aes(x = k_R, y = Mean,colour = "0"), linewidth = lnwd) +
  geom_vline(aes(xintercept = gamma_R_star.list[[1]][[2]],colour = "0"), linetype = "dashed") +
  geom_line(data = gamma_R_star.list[[2]][[1]], aes(x = k_R, y = Mean,colour = "0.2"), linewidth = lnwd) +
  geom_vline(aes(xintercept = gamma_R_star.list[[2]][[2]], colour = "0.2"), linetype = "dashed") +
  geom_line(data = gamma_R_star.list[[3]][[1]], aes(x = k_R, y = Mean,colour = "0.4"), linewidth = lnwd) +
  geom_vline(aes(xintercept = gamma_R_star.list[[3]][[2]], colour = "0.4"), linetype = "dashed") +
  geom_line(data = gamma_R_star.list[[4]][[1]], aes(x = k_R, y = Mean,colour = "0.6"), linewidth = lnwd) +
  geom_vline(aes(xintercept = gamma_R_star.list[[4]][[2]], colour = "0.6"), linetype = "dashed") +
  geom_line(data = gamma_R_star.list[[5]][[1]], aes(x = k_R, y = Mean,colour = "0.8"), linewidth = lnwd) +
  geom_vline(aes(xintercept = gamma_R_star.list[[5]][[2]], colour = "0.8"), linetype = "dashed") +
  geom_line(data = gamma_R_star.list[[6]][[1]], aes(x = k_R, y = Mean,colour = "1"), linewidth = lnwd) +
  geom_vline(aes(xintercept = gamma_R_star.list[[6]][[2]], colour = "1"), linetype = "dashed") +
  xlab("Transcription rate (k_R)") +
  ylab("mRNA mean expression") +
  xlim(40,275) +
  labs(colour = "Rate of complex formation") +
  ggtitle(paste("miRNA = 250"))

grid.arrange(grobs = list(gamma_R_star.noise.plot,gamma_R_star.mean.plot), ncol = 2, nrow = 1)


#noisy gene-----
parms.noisy <- c(k_R = 1, #transcription rate
             k_on = 2.5, #rate at which free miRNAs are removed from the system 
             k_off = 0.05, #unbinding of the miRNA from its target
             gamma_R = 0.25, #rate of intrinsic decay 
             gamma_R_star = 0.5, #destruction of miRNAs target
             miRNA_0 = 10
)

miRNA_seq <- c(0,25,50,100,150,250)
values.seq <- seq(0,100, by = 2.5)
noisy.list <- variable_miRNA(parms = parms.noisy, levels = miRNA_seq, parameter = "k_R", values = values.seq)


noisy.noise.plot <- ggplot() +
  geom_line(data = noisy.list[[1]][[1]], aes(x = k_R, y = Noise,colour = "0"), linewidth = lnwd) +
  geom_vline(aes(xintercept = noisy.list[[1]][[2]],colour = "0"), linetype = "dashed") +
  geom_line(data = noisy.list[[2]][[1]], aes(x = k_R, y = Noise,colour = "25"), linewidth = lnwd) +
  geom_vline(aes(xintercept = noisy.list[[2]][[2]], colour = "25"), linetype = "dashed") +
  geom_line(data = noisy.list[[3]][[1]], aes(x = k_R, y = Noise,colour = "50"), linewidth = lnwd) +
  geom_vline(aes(xintercept = noisy.list[[3]][[2]], colour = "50"), linetype = "dashed") +
  geom_line(data = noisy.list[[4]][[1]], aes(x = k_R, y = Noise,colour = "100"), linewidth = lnwd) +
  geom_vline(aes(xintercept = noisy.list[[4]][[2]], colour = "100"), linetype = "dashed") +
  geom_line(data = noisy.list[[5]][[1]], aes(x = k_R, y = Noise,colour = "150"), linewidth = lnwd) +
  geom_vline(aes(xintercept = noisy.list[[5]][[2]], colour = "150"), linetype = "dashed") +
  geom_line(data = noisy.list[[6]][[1]], aes(x = k_R, y = Noise,colour = "250"), linewidth = lnwd) +
  geom_vline(aes(xintercept = noisy.list[[6]][[2]], colour = "250"), linetype = "dashed") +
  xlab("Transcription rate (k_R)") +
  xlim(0,100) +
  ylab("Gene noise") +
  scale_colour_manual(
    values = c("0" = "blue", "25" = "red", "50" = "darkgreen", 
               "100" = "purple", "150" = "orange", "250" = "brown"),
    breaks = miRNA_seq  # Order in the legend
  ) +
  labs(colour = "Total miRNA") +
  ggtitle(paste("Noisy gene"))

noisy.mean.plot <- ggplot() +
  geom_line(data = noisy.list[[1]][[1]], aes(x = k_R, y = Mean,colour = "0"), linewidth = lnwd) +
  geom_vline(aes(xintercept = noisy.list[[1]][[2]],colour = "0"), linetype = "dashed") +
  geom_line(data = noisy.list[[2]][[1]], aes(x = k_R, y = Mean,colour = "25"), linewidth = lnwd) +
  geom_vline(aes(xintercept = noisy.list[[2]][[2]], colour = "25"), linetype = "dashed") +
  geom_line(data = noisy.list[[3]][[1]], aes(x = k_R, y = Mean,colour = "50"), linewidth = lnwd) +
  geom_vline(aes(xintercept = noisy.list[[3]][[2]], colour = "50"), linetype = "dashed") +
  geom_line(data = noisy.list[[4]][[1]], aes(x = k_R, y = Mean,colour = "100"), linewidth = lnwd) +
  geom_vline(aes(xintercept = noisy.list[[4]][[2]], colour = "100"), linetype = "dashed") +
  geom_line(data = noisy.list[[5]][[1]], aes(x = k_R, y = Mean,colour = "150"), linewidth = lnwd) +
  geom_vline(aes(xintercept = noisy.list[[5]][[2]], colour = "150"), linetype = "dashed") +
  geom_line(data = noisy.list[[6]][[1]], aes(x = k_R, y = Mean,colour = "250"), linewidth = lnwd) +
  geom_vline(aes(xintercept = noisy.list[[6]][[2]], colour = "250"), linetype = "dashed") +
  xlab("Transcription rate (k_R)") +
  xlim(0,100) +
  ylab("mRNA mean expression") +
  scale_colour_manual(
    values = c("0" = "blue", "25" = "red", "50" = "darkgreen", 
               "100" = "purple", "150" = "orange", "250" = "brown"),
    breaks = miRNA_seq  # Order in the legend
  ) +
  labs(colour = "Total miRNA") +
  ggtitle(paste("Noisy gene"))

grid.arrange(grobs = list(noisy.noise.plot,noisy.mean.plot), ncol = 2, nrow = 1)