#Function to plot alpha diversity

plot_alpha <- function(df, response, x_axis, annotation_df = NULL, est_mean = NULL, 
                       point = "no", violin = "no", box = "no", 
                       fill_var = NA, color_var = NA, textsize = alpha_textsize, annot = NULL, annot_size = 4){
  
  #Convert axes strings to symbols for evaluation
  response_sym <- sym(response)
  x_axis_sym <- sym(x_axis)
  
  #Convert color and fill variables to strings if they are provided
  if(!is.na(fill_var)){
    fill_var_sym <- sym(fill_var)
    
  }else{
    fill_var_sym <- NA
  }
  
  if(!is.na(color_var)){
    color_var_sym <- sym(color_var)
  }else{
    color_var_sym <- NA
  }
  
  #Set up base plot
  plot <- df %>%
    ggplot(aes(x = !!x_axis_sym, y = !!response_sym)) +
    theme_classic() +
    textsize +
    labs(x = "Treatment")+
    scale_fill_manual(values = cbpalette) +
    scale_color_manual(values = cbpalette)
  
  #Add raw data points if indicated
  if (point == "yes") {
    plot <- plot + geom_point(aes(color = !!color_var_sym), size = 3)
  }
  
  #Add violin plot if indicated
  if (violin == "yes") {
    plot <- plot + geom_violin(color = NA, aes(fill = !!fill_var_sym))
  }
  
  #Add estimated marginal means if indicated
  if (!is.null(est_mean)) {
    em_df <- est_mean %>%
      filter(metric == response) %>%
      rename(!!response_sym := emmean)
    
    # #Add horizontal bars representing mean
    # plot <- plot + 
    #   geom_errorbarh(data = em_df, aes(y = !!response_sym, xmin = !!x_axis_sym - 0.3, xmax = !!x_axis_sym + 0.3), color = "black", height = 0)
    
    #Add points representing means
    plot <- plot +
      geom_point(data = em_df, aes(y = !!response_sym, color = !!color_var_sym), size = 3) +
      geom_linerange(data = em_df, aes(ymax = upper.CL, ymin = lower.CL, color = !!color_var_sym))
  }
  
  #Add boxplots if indicated
  
  if(box == "yes") {
    plot <- plot +
      geom_boxplot(aes(color = !!color_var_sym))
  }
  
  #Add annotations if indicated
  if(!is.null(annot)){
    plot <- plot +
      annotate(geom = "text", y = max(df[,response]) + 0.05*max(df[,response]), x = 3.5, 
               label = filter(annot, metric == response)$label, size = annot_size)
  }
  return(plot)
}