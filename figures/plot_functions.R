#calculates parameter CVs and estimates across runs and scenarios 

library('matrixStats')
library('tidyverse')
library('reshape2')
library('readr')
library(scales)
library(patchwork)
library(tidyverse)
library(stats)
library(readr)

extract_columns <- function(runname, column_names) {

  
  runname = runname
  dir.run=paste0(dir.out,"/",runname)
  run_list = read_csv(paste0(dir.run,'/',runname,'-run_list.csv'))
  
  required_cols <- c("marks", "trial", "run.num", "full.dir")
  if (!all(required_cols %in% names(run_list))) {
    stop("run_list must contain: marks, trial, run.num, full.dir")
  }
  
  run_list %>%
    select(marks, trial, run.num, full.dir) %>%
    rowwise() %>%
    mutate(
      extracted = list({
        file_path <- file.path(full.dir, "tseries.csv")
        
        if (!file.exists(file_path)) {
          warning(paste("File not found:", file_path))
          return(setNames(vector("list", length(column_names)), column_names))
        }
        
        tseries <- read_csv(file_path, show_col_types = FALSE)
        
        # Extract requested columns, returning NA for any missing
        extracted_data <- map(column_names, function(col) {
          if (col %in% names(tseries)) {
            tseries[[col]]
          } else {
            warning(paste("Column", col, "missing in", file_path))
            NA
          }
        })
        
        setNames(extracted_data, column_names)
      })
    ) %>%
    unnest_wider(extracted) %>%
    ungroup()
}


plot_time_series <- function(results, y_columns, color_by_marks = FALSE) {
  required_cols <- c("marks", "trial", "run.num", "tstep")
  if (!all(required_cols %in% names(results))) {
    stop("The results tibble must contain: marks, trial, run.num, and tstep.")
  }
  
  missing_cols <- setdiff(y_columns, names(results))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in results:", paste(missing_cols, collapse = ", ")))
  }
  
  # Unnest everything
  long_data <- results %>%
    select(marks, run.num, tstep, all_of(y_columns)) %>%
    unnest(cols = c(tstep, all_of(y_columns))) %>%
    pivot_longer(cols = all_of(y_columns), names_to = "variable", values_to = "value")
  
  # Summarize: by tstep and variable (+ marks if grouping)
  if (color_by_marks) {
    summary_data <- long_data %>%
      group_by(marks, variable, tstep) %>%
      summarise(
        mean = mean(value, na.rm = TRUE),
        lower = quantile(value, 0.025, na.rm = TRUE),
        upper = quantile(value, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
    
    ggplot(summary_data, aes(x = tstep, y = mean, color = factor(marks), fill = factor(marks))) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
      geom_line(size = 1) +
      facet_wrap(~ variable, scales = "free_y") +
      labs(
        x = "Time Step",
        y = "Mean ± 95% CI",
        color = "Marks",
        fill = "Marks",
        title = "Time Series Mean and 95% CI by Marks"
      ) +
      theme_minimal()
    
  } else {
    # No grouping by marks: pooled across all runs
    summary_data <- long_data %>%
      group_by(variable, tstep) %>%
      summarise(
        mean = mean(value, na.rm = TRUE),
        lower = quantile(value, 0.025, na.rm = TRUE),
        upper = quantile(value, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
    
    ggplot(summary_data, aes(x = tstep, y = mean)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "darkgrey", alpha = 0.3) +
      geom_line(color = "black", size = 1) +
      facet_wrap(~ variable, scales = "free_y") +
      labs(
        x = "Time Step",
        y = "Mean ± 95% CI",
        title = "Time Series Mean and 95% CI Across Runs"
      ) +
      theme_minimal()
  }
}

cv_calc = function(runname = 'update-srfs-tagging-95rep', Nmarks, p.trips, tag_type){
  dir.run=paste0(dir.out,"/",runname)
  
  #import run.list.csv
  run.list = read_csv(paste0(dir.run,'/',runname,'-run_list.csv'))
  
  #imports est results into list
  results = lapply(paste0(run.list$full.dir,'/est_results.csv'),function(i){
    read.csv(i, header = T)
  })
  names(results) = run.list$fold
  
  #extract age estimates 
  res.df = matrix(NA,length(results[[1]]$estimate),length(run.list$marks))
  for(i in 1:length(run.list$marks)){
    res.df[,i] <- as.numeric(results[[i]]$estimate)
  }
  res.df = data.frame(res.df)
  colnames(res.df) = run.list$fold
  
  #extract SE for esimates to check for convergence 
  converge.df = matrix(NA,length(results[[1]]$se),length(run.list$marks))
  for(i in 1:length(run.list$marks)){
    converge.df[,i] <- as.numeric(results[[i]]$se)
  }
  converge.df = data.frame(converge.df)
  #adds parameter names 
  converge.df = cbind(as.character(results[[1]]$X), converge.df)
  colnames(converge.df) = c("formula", run.list$fold)
  
  #can select specific runs that didn't converge or produce unrealistic results 
  res.df = res.df %>% select(-grep("run_52|run_47", names(res.df), value = TRUE))
  
  
  #index estimates based on ntrials and nmarks 
  # Can use the same index to extract the SEs to check for convergence (for percent converge mtext on plot)
  if(tag_type == 'conv'){
    for(f in 1:length(Nmarks)){
      sel_cols = grep(paste0("Nmarks_",Nmarks[f]),colnames(res.df), value = TRUE)
      tmp = res.df[,sel_cols]
      tmp.df = res.df[,sel_cols]
      tmp$mean = rowMeans(tmp.df)# will throw error if less than 2 runs per tag converged
      #tmp$mean = exp(rowMeans(log(tmp.df)))
      tmp$median = apply(tmp.df, 1, median, na.rm=T)
      tmp$sd <- apply(tmp.df, 1, sd, na.rm=TRUE)
      tmp$formula <- as.character(results[[1]]$X)
      tmp$cv <- tmp$sd/tmp$mean
      tmp$ages <- substr(tmp$formula,unlist(gregexpr(pattern ='a',tmp$formula))+1,unlist(gregexpr(pattern ='a',tmp$formula))+2)
      tmp$seasons <- substr(tmp$formula,unlist(gregexpr(pattern ='t',tmp$formula))+1,unlist(gregexpr(pattern ='t',tmp$formula))+2)
      tmp$par <- substr(tmp$formula,1,1)
      #need to adjust given the selected age bins
      #tmp$age_bin <- c(0:6,1:6,0:6,0,0)
      #tmp$age_bin <- ifelse(tmp$ages%in%c(12:24),'[1,2]',ifelse(tmp$ages%in%c(25:60), '(2,6]','(6,25]')) #age 5+
      tmp$age_bin <- ifelse(tmp$ages%in%c(12:24),'[1,2]',ifelse(tmp$ages%in%c(25:72), '(2,6]','(6,25]')) #age 6+
      #need to adjust given the selected time bins
      tmp$time_bin <- ifelse(tmp$seasons%in%c(1:4),'preseason',ifelse(tmp$seasons%in%c(5:8), 'fishing','postseason')) 
      #tmp$time_bin <- ifelse(tmp$seasons%in%c(1:6),'preseason',ifelse(tmp$seasons%in%c(7:9), 'fishing','postseason'))
      #tmp$time_bin = ifelse(tmp$seasons%in%c(1:4,10:12),0,ifelse(tmp$seasons%in%c(7),2,1)) #0 = low pressure, 1 = high pressure, 2= july
      #tmp$time_bin = ifelse(tmp$seasons%in%c(1:4,10:12),0,1) #0 = low pressure, 1 = high pressure
      #tmp$time_bin = ifelse(tmp$seasons%in%c(1:4,9:12, 13:16,21:24),"low",ifelse(tmp$seasons%in%c(7, 19),"july","high")) #2 yr, 0 = low pressure, 1 = high pressure, 2= july
      #tmp$time_bin = ifelse(tmp$seasons%in%c(1:16),"tagging",ifelse(tmp$seasons%in%c(17:20),"fishing","post")) #2 yr, 0 = low pressure, 1 = high pressure, 2= july
      
      assign(paste("nmarks",Nmarks[f],'_CV',sep=""), tmp, envir = .GlobalEnv)
    }
    #fixes name for plotting 
    nmarks0500_CV <<- nmarks500_CV
    rm(nmarks500_CV, envir = .GlobalEnv)
  }
  if(tag_type == "gene"){
    for(f in 1:length(p.trips)){
      sel_cols = grep(paste0("p.trips_",p.trips[f]),colnames(res.df), value = TRUE)
      tmp = res.df[,sel_cols]
      tmp$mean = rowMeans(tmp)# will throw error if less than 2 runs per tag converged 
      tmp$sd <- apply(tmp, 1, sd, na.rm=TRUE)
      tmp$formula <- as.character(results[[1]]$X)
      tmp$cv <- tmp$sd/tmp$mean
      tmp$ages <- substr(tmp$formula,unlist(gregexpr(pattern ='a',tmp$formula))+1,unlist(gregexpr(pattern ='a',tmp$formula))+2)
      tmp$seasons <- substr(tmp$formula,unlist(gregexpr(pattern ='t',tmp$formula))+1,unlist(gregexpr(pattern ='t',tmp$formula))+2)
      tmp$par <- substr(tmp$formula,1,1)
      #need to adjust given the selected age bins
      #tmp$age_bin <- c(0:6,1:6,0:6,0,0)
      tmp$age_bin <- ifelse(tmp$ages%in%c(12:24),'[1,2]',ifelse(tmp$ages%in%c(25:60), '(2,6]','(6,25]')) #age 5+
      #tmp$age_bin <- ifelse(tmp$ages%in%c(12:24),'[1,2]',ifelse(tmp$ages%in%c(25:72), '(3,6]','(6,25]')) #age 6+
      #need to adjust given the selected time bins
      tmp$time_bin <- ifelse(tmp$seasons%in%c(1:4),'preseason',ifelse(tmp$seasons%in%c(5:8), 'fishing','postseason'))
      #tmp$time_bin <- ifelse(tmp$seasons%in%c(1:6),'preseason',ifelse(tmp$seasons%in%c(7:9), 'fishing','postseason'))
      #tmp$time_bin = ifelse(tmp$seasons%in%c(1:4,9:12, 13:16,21:24),0,ifelse(tmp$seasons%in%c(7),2,1)) #0 = low pressure, 1 = high pressure, 2= july
      #tmp$time_bin = ifelse(tmp$seasons%in%c(1:4,10:12),0,1) #0 = low pressure, 1 = high pressure
      assign(paste("p.trips",p.trips[f],'_CV',sep=""), tmp)
    }
  }
  
  #produces data frame of CVs by bin and parameter
  if(tag_type=="conv"){
    #dat_list <- mget(x = ls(pattern = "nmarks"))
    #dat_list = mget(Nmarks, envir=parent.frame(), inherits=TRUE)
    dat_list <<- mget(ls(pattern = "nmarks", envir = .GlobalEnv), envir = .GlobalEnv)
    
    cv1 = matrix(NA,length(results[[1]]$se),length(Nmarks))
    for(r in 1:length(Nmarks)) {
      cv1[,r] = dat_list[[r]][["cv"]]
    }
    cv1 = data.frame(cv1)
    cv1 = cbind(dat_list[[2]][["par"]], dat_list[[2]][["age_bin"]],dat_list[[2]][["time_bin"]],cv1)
    colnames(cv1) = c('par','age','time',Nmarks) #"p.trips or nmarks
  }
  if(tag_type=="gene"){
    dat_list <- mget(x = ls(pattern = "p.trips"))
    cv1 = matrix(NA,length(results[[1]]$se),length(p.trips))
    for(r in 1:length(p.trips)) {
      cv1[,r] = dat_list[[r+1]][["cv"]]
    }
    cv1 = data.frame(cv1)
    cv1 = cbind(dat_list[[2]][["par"]], dat_list[[2]][["age_bin"]],dat_list[[2]][["time_bin"]],cv1)
    colnames(cv1) = c('par','age','time',p.trips) #"p.trips or nmarks
    
  }
  return(cv1)
}

cv_plot = function(cv_dat, para = 'R', tag_type){
  s.cv1 = cv_dat[cv_dat$par %in% para,]
  
  if(para != 'r'){
  if(tag_type=="conv"){
    #CV plots for R and S
    cv_long <- s.cv1 %>%
      pivot_longer(
        cols = c(`500`, `1000`, `1500`, `2000`, `2500`, `3000`),
        names_to = "Nmark",
        values_to = "CV"
      ) %>%
      mutate(
        Nmark = as.numeric(Nmark),
        time = factor(time, levels = c("preseason", "fishing", "postseason")),
        #age = gsub("\\[|\\]", "", age),  # remove brackets
        age = factor(age, levels = c("[1,2]", "(2,6]", "(6,25]"), labels = c("Age Bin: [1,2]", "Age Bin: (2,6]", "Age Bin: (6,25]"))  # relevel & relabel
      )
    
    p.cv = ggplot(cv_long, aes(x = Nmark, y = CV, color = time, group = time)) +
      geom_line(size = 1.2) +
      facet_wrap(~ age, ncol = 1, strip.position = "top") +
      geom_hline(yintercept = c(0.2, 0.3), linetype = "dashed", color = "grey50") +
      scale_color_manual(
        values = c("preseason" = "red", "fishing" = "blue", "postseason" = "black"),
        labels = c("Preseason", "Fishing", "Post Season")
      ) +
      scale_y_continuous(limits = c(0, .75), expand = c(0, 0.01)) +
      labs(
        x = "Number of Marks",
        y = "CV"
      ) +
      theme_classic(base_size = 13) + scale_x_continuous() +
      theme(
        # Remove default grid lines
        panel.grid = element_blank(),
        # Y label same size as ticks
        axis.title.y = element_text(size = 13),
        # Strip (facet title)
        strip.text = element_text(size = 14, face = "bold"),  strip.background = element_blank(),
        # Legend
        legend.position = "bottom",
        legend.title = element_blank(),
        # Plot border
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        # No grid on x-axis
        axis.line = element_blank()
      ) + 
      scale_x_continuous(breaks = unique(cv_long$Nmark))
    #ggsave("figures/pub_figs_v2/update_Scv_plot_rep.png", width = 6, height = 7, dpi = 600)
    
  }
  
  
  if(tag_type=="gene"){
    
    #s.cv1 = cbind(age.bin,s.cv1)
    
    #CV plots for R and S
    cv_long <- s.cv1 %>%
      pivot_longer(
        #cols = c('0.001', '0.005', '0.01'),
        cols = c('0.02', '0.025', "0.03"),
        names_to = "p.trips",
        values_to = "CV"
      ) %>%
      mutate(
        p.trips = as.numeric(p.trips),
        time = factor(time, levels = c("preseason", "fishing", "postseason")),
        #age = gsub("\\[|\\]", "", age),  # remove brackets
        age = factor(age, levels = c("[1,2]", "(2,6]", "(6,25]"), labels = c("Age Bin: [1,2]", "Age Bin: (2,6]", "Age Bin: (6,25]"))  # relevel & relabel
      )
    
    p.cv = ggplot(cv_long, aes(x = p.trips, y = CV, color = time, group = time)) +
      geom_line(size = 1.2) +
      facet_wrap(~ age, ncol = 1, strip.position = "top") +
      geom_hline(yintercept = c(0.2, 0.3), linetype = "dashed", color = "grey50") +
      scale_color_manual(
        values = c("preseason" = "red", "fishing" = "blue", "postseason" = "black"),
        labels = c("Preseason", "Fishing", "Post Season")
      ) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0.01)) +
      labs(
        x = "Participating Trips (%)",
        y = "CV"
      ) +
      theme_classic(base_size = 13) + scale_x_continuous() +
      theme(
        # Remove default grid lines
        panel.grid = element_blank(),
        # Y label same size as ticks
        axis.title.y = element_text(size = 13),
        # Strip (facet title)
        strip.text = element_text(size = 14, face = "bold"),  strip.background = element_blank(),
        # Legend
        legend.position = "bottom",
        legend.title = element_blank(),
        # Plot border
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        # No grid on x-axis
        axis.line = element_blank()
      ) + 
      scale_x_continuous(breaks = unique(cv_long$p.trips), labels = c("2", "2.5", "3"))
    #ggsave("figures/pub_figs_v2/updateRcv_plot_gene.png", width = 6, height = 7, dpi = 600)
    
  }
  }
  if(para == 'r'){
    cv_long <- s.cv1 %>%
      pivot_longer(
        cols = c(`500`, `1000`, `1500`, `2000`, `2500`, `3000`),
        names_to = "Nmark",
        values_to = "CV"
      ) %>%
      mutate(
        Nmark = as.numeric(Nmark),
        time = factor(time, levels = c("fishing")),
        #age = gsub("\\[|\\]", "", age),  # remove brackets
        age = factor(age, levels = c( "(2,6]", "(6,25]"), labels = c( "Age Bin: (2,6]", "Age Bin: (6,25]")) 
      )%>%
      filter(!is.na(age), !is.na(time), CV != 0)
    
    p.cv = ggplot(cv_long, aes(x = Nmark, y = CV, color = time, group = time)) +
      geom_line(size = 1.2) +
      facet_wrap(~ age, ncol = 1, strip.position = "top") +
      geom_hline(yintercept = c(0.2, 0.3), linetype = "dashed", color = "grey50") +
      scale_color_manual(
        values = c("fishing" = "blue"),
        labels = c("Fishing")
      ) +
      scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0.01)) +  # expanded limit to fit your higher CVs
      labs(
        x = "Number of Marks",
        y = "CV"
      ) +
      theme_classic(base_size = 13) +
      scale_x_continuous(breaks = unique(cv_long$Nmark)) +
      theme(
        panel.grid = element_blank(),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank()
      )
    #ggsave("figures/pub_figs_v2/update_lilrcv_plot_95rep.png", width = 6, height = 7, dpi = 600)
    
  }
  
  p.cv
  
}

# Estimate data and plots 

extract_runs_filtered <- function(df, nmarks_label, par_value = "S") {
  df <- df %>% filter(par == par_value)
  
  run_cols <- grep("^Nmarks_", names(df), value = TRUE)
  nmarks_val <- as.numeric(stringr::str_extract(nmarks_label, "\\d+"))
  
  df_long <- tidyr::pivot_longer(
    df,
    cols = all_of(run_cols),
    names_to = "run",
    values_to = "est"
  ) %>%
    dplyr::mutate(
      Nmarks = nmarks_val,
      par = as.character(par),
      age_bin = as.character(age_bin),
      time_bin = as.character(time_bin)
    )
  
  return(df_long)
}

est_plots = function(dat_list, para = 'R', Nmarks,  tag_type = 'conv'){

# Combine all list elements
cv_box_data <- purrr::map_dfr(
  names(dat_list),
  ~ extract_runs_filtered(dat_list[[.x]], .x, par_value = para)
)

if(para != 'r'){
  
cv_box_data <- cv_box_data %>%
  mutate(
    age_bin = recode(age_bin,
                     "[1,2]" = "[1,2]",
                     "(2,6]" = "(2,6]",
                     "(6,25]" = "(6,25]"),
    age_bin = factor(age_bin, levels = c("[1,2]", "(2,6]", "(6,25]"), labels = c("Age Bin: [1,2]", "Age Bin: (2,6]", "Age Bin: (6,25]")),
    time_bin = factor(time_bin, levels = c("preseason", "fishing", "postseason"))
  )

# Convert Nmarks to factor if not already (for discrete x-axis)
cv_box_data$Nmarks <- factor(cv_box_data$Nmarks)
# Calculate y-position for labels — slightly above the max estimate per tag number
label_positions <- cv_box_data %>%
  group_by(Nmarks) %>%
  summarize(y_pos = max(est, na.rm = TRUE) + 0.02 * diff(range(est)))

p.est = ggplot(cv_box_data, aes(x = factor(Nmarks), y = est, fill = time_bin)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~age_bin, ncol = 1)  +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),  
    axis.ticks.y = element_line()       
  ) +  #geom_hline(
  #aes(yintercept = geo_mean_val, color = time_bin),
  #linetype = "dashed", size = 1, inherit.aes = F) + 
  #scale_y_continuous(limits = c(0, .15), expand = c(0, 0.01)) + 
  scale_fill_manual(values = c("preseason" = "red", "fishing" = "blue", "postseason" = "darkgrey"), 
                    labels = c("Preseason", "Fishing", "Post Season")) +
  scale_color_manual(values = c("preseason" = "red", "fishing" = "blue", "postseason" = "darkgrey"), 
                     labels = c("Preseason", "Fishing", "Post Season")) +
  labs(x = "Number of Marks", y = "Parameter Estimate", fill = "Time Bin ") 
#ggsave("figures/pub_figs_v2/update_Sest_plot_95rep_bias.png", width = 6, height = 8, dpi = 600)
}

if(para == 'r'){
# Filter to only fishing and two older age bins
cv_box_data <- cv_box_data %>%
  filter(time_bin == "fishing", age_bin %in% c("(2,6]", "(6,25]")) %>%
  mutate(
    age_bin = factor(age_bin, levels = c("(2,6]", "(6,25]"), labels = c("Age Bin: (2,6]", "Age Bin: (6,25]")),
    Nmarks = factor(Nmarks)
  )

p.est = ggplot(cv_box_data, aes(x = Nmarks, y = est, fill = time_bin)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~ age_bin, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("fishing" = "blue")) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_blank(),           # No background in strip
    strip.text = element_text(face = "bold", size = 13),  # Match line plot facet titles
    axis.title = element_text(size = 13),         # Match axis title size
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Number of Marks",
    y = "Parameter Estimate",
    fill = NULL
  )
#ggsave("figures/pub_figs_v2/update_lilrest_plot_95rep.png", width = 6, height = 8, dpi = 600)
}

p.est

}


#plot functions 
#https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}