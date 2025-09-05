####################################
###    Supplementary Figures     ###
####################################

datpath <- "YourDir"

### Figure S1-S3 ###

#Same script as Figure 2 from the main manuscript 

### Figure S4-S7 ###

# city of interest
cities<-c("abbotsford", "victoria", "lytton")

for (j in 1:length(cities)){
  #j<-1
  city <- cities[[j]]
  
  # Difference in age?
  var<-c("all_mortality_0_64", "all_mortality_65_74", "all_mortality_75")
  label_var<-c("<65 years", "65-74 years", ">74 years")
  
  res_list<-list()
  #par(mfrow = c(2, 2))
  for(i in 1:3){
    #i<-2
    load(paste0(datpath, "results/DLNM/DLNM_Res_summer3dayslag_MMT50to100_", 
                var[i], "_",city, ".RData"))
    res_list[[i]]<-data.frame(temp=pred$predvar, allRRfit=pred$allRRfit, allRRlow=pred$allRRlow, 
                              allRRhigh=pred$allRRhigh, age=rep(label_var[i],length(pred$predvar)), 
                              cen=rep(pred$cen,length(pred$cen)))
    print(pred$cen)
  }
  
  dta_plot<-do.call("rbind", res_list)
  dta_plot$age<-factor(dta_plot$age, levels = label_var)
  
  age_plot<-ggplot(dta_plot, aes(temp, allRRfit, color = age)) + 
    geom_line(linewidth=0.8) +
    geom_ribbon(aes(temp, ymin = allRRlow, ymax = allRRhigh, fill = age), alpha = .2, colour=NA) +
    xlab("Temperature") +
    ylab("RR") +
    ylim(0.5,4) +
    labs(title="Age") +
    theme(plot.title = element_text(size=10), 
          axis.title=element_text(size=10),
          legend.text=element_text(size=10),
          legend.title=element_blank()) +
    scale_fill_colorblind()
  
  #age_plot <- age_plot + guides(fill=guide_legend(title="age (years)"))
  
  # Difference in gender?
  var<-c("all_mortality_male", "all_mortality_female")
  label_var<-c("male", "female")
  
  res_list<-list()
  #par(mfrow = c(2, 2))
  for(i in 1:2){
    #i<-2
    load(paste0(datpath, "results/DLNM/DLNM_Res_summer3dayslag_MMT50to100_", 
                var[i], "_",city, ".RData"))
    res_list[[i]]<-data.frame(temp=pred$predvar, allRRfit=pred$allRRfit, allRRlow=pred$allRRlow, 
                              allRRhigh=pred$allRRhigh, sex=rep(label_var[i],length(pred$predvar)))
  }
  
  dta_plot<-do.call("rbind", res_list)
  
  sex_plot<-ggplot(dta_plot, aes(temp, allRRfit, color = sex)) + 
    geom_line(linewidth=0.8) +
    geom_ribbon(aes(temp, ymin = allRRlow, ymax = allRRhigh, fill = sex), alpha = .2, colour=NA) +
    xlab("Temperature") +
    ylab("RR") +
    ylim(0.5,4) +
    labs(title="Sex") +
    theme(plot.title = element_text(size=10), 
          axis.title=element_text(size=10),
          legend.text=element_text(size=10),
          legend.title=element_blank()) +
    scale_fill_colorblind()
  
  ggsave(
    paste0(datpath, "results/plots/Canada/DLNM/Age_sex_difference_", city, ".png"),
    ggarrange(age_plot, sex_plot+ rremove("ylab"), 
              labels = c("A", "B"),
              ncol = 2, nrow = 1, legend ="bottom",
              font.label = list(size = 8, color = "black", face = "bold", family = NULL))
  )
  #grDevices::dev.set(1)
}

### Figure S8-S9 ###

# same script as Figure 4 and 5 from the main manuscript
