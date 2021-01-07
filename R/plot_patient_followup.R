#' Plot patient followup
#' 
#' Create an example plot of patient followup for the manuscript
#' @return A [ggplot2][ggplot2::ggplot2] object.
#' @export
plot_patient_followup <- function(){
  dat <- data.frame(ID = c(1:7), event = c(0, 1, 1, 1, 0, 0, 0),  # c(7, 1, 3, 4, 2, 6, 5)
                    t1 = c(1991, 2008, 1993, 2007, 2012, 2014, 1989),
                    t2 = c(1998, 2009, 2013, 2017, 2019, 2019, 2019), 
                    #t1 = c(1989, 1991, 1993, 2007, 2008, 2012, 2014), 
                   # t2 = c(2019, 1998, 2013, 2017, 2009, 2019, 2019), 
                   censored = c(0, 1, 0, 0, 1, 1, 1))
                   # censored = c(1, 0, 0, 0, 1, 1, 1)) #, .Names = c("ID", "Death", "t1", "t2", "Censored")) #, class = "data.frame", row.names = c(NA, -6L))
  
  # Create event variable
  #dat$Death <- with(dat, ifelse(Death, "Dead", "Censored"))
  
  # Create id.ordered, which is a factor that is ordered by t2
  # This will allow the plot to be ordered by increasing t2, if desired
  dat$id.ordered <- factor(x = dat$ID, levels = c(7,6,5,4,3,2,1)) #order(dat$t2, decreasing = T))
  
  # Use ggplot to plot data from dat object
  ggplot(dat %>% 
           dplyr::mutate(censored=ifelse(censored==1,"Censored","Dead")) %>%
           dplyr::mutate(left_truncated=ifelse(t2<2011,"Left truncated","Included in study")),
         aes(x = id.ordered),group=id.ordered) + 
    # Plot solid line representing non-interval censored time from 0 to t1
    geom_linerange(aes(ymin = t1, ymax = t2, col=left_truncated)) + 
    geom_point(aes(id.ordered, t2,shape=factor(censored),col=left_truncated)) + 
    geom_point(aes(id.ordered, t1, col=left_truncated)) +
    # Plot line (dotted for censored time) representing time from t1 to t2
    #geom_linerange(aes(ymin = t1, ymax = t2, linetype = as.factor(Censored))) +  
    # Plot points representing event
    # Flip coordinates
    coord_flip() + 
    scale_color_manual("Left truncation",values=c("black","indianred"))+
    # Add custom shape scale.  Change the values to get different shapes.
    scale_shape_manual(name = "Event", values = c(1, 4)) +
    # Add main title and axis labels
    xlab("Patient ID") +  ylab("Year of diagnosis") + 
    geom_hline(yintercept = 2011, linetype="dashed", col="dodgerblue") +
    annotate("text", x = 7, y = 2014.5, 
             label = "Study start", color="dodgerblue", 
             size=4 , angle=0, fontface="bold")+
    # I think the bw theme looks better for this graph, 
    # but leave it out if you prefer the default theme
    theme_bw()
}