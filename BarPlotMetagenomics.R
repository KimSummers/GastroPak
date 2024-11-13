# BarPlotMetagenomics
#
# Bar chart plot metagenomic data
#
# file BarPlotMetagenomics
#
# inputs
# 	graphData         - dataframe containing main data
# 	graphData2        - dataframe containing secondary (HuBac/RuBac) data
#   xAxisLabels       - Labels for the x-axis
#   absolute          - TRUE to plot absolute abundances, FALSE for percentage
#   absLab            - y label for absolute date (default is "Abundance per million reads")
#   graphType         - "Abundances" for ARGs or Species, "HuBac" for HuBac/RuBac data only, "Both" to plot both
#   yAxisVals         - Optional parameter to override y-axis limits
#   plotsDir          - directory to store plots in
#   plotTitle         - title for the plot
#
# Version    Author       Date      Affiliation
# 1.00       J K Summers  14/12/23  Wellington Lab - School of Life Sciences - University of Warwick
BarPlotMetagenomics <- function(graphData, graphData2, xAxisLabels, absolute,
                                absLab = "Abundance per million reads",
                                graphType, yAxisVals = NULL, plotsDir, plotTitle) {

  library(pals)

  mgPlot <- ggplot(data = graphData,
                   aes(x = factor(xVals, level = unique(xVals)), y = Abundance))

  if (graphType == "HuBac")
  {
    huBacData <- graphData[graphData$Bac == "HuBac", ]
    ruBacData <- graphData[graphData$Bac == "RuBac", ]

    mgPlot <- mgPlot + geom_point(data = huBacData, 
                                  aes(y = Abundance), colour = "Red", size = 4)
    mgPlot <- mgPlot + geom_point(data = ruBacData, 
                                  aes(y = Abundance), colour = "Blue", size = 4)

    if (!is.null(yAxisVals))
    {
      mgPlot <- mgPlot + ylim(yAxisVals)
    }
    
  }else
  {

    if (absolute)
    {
      mgPlot <- mgPlot + geom_col(aes(fill = ARG))
      mgPlot <- mgPlot + labs(y = absLab)
    }else
    {
      mgPlot <- mgPlot + geom_col(position = "fill",  aes(fill = ARG))
      mgPlot <- mgPlot + labs(y = "Proportion of abundance")
    }

    mgPlot <- mgPlot + guides(color = guide_legend(override.aes = list(size = 10)))
    mgPlot <- mgPlot + scale_fill_manual(values = glasbey(), drop = FALSE)

    if (graphType == "Both")
    {
      huBacData <- graphData2[graphData2$Bac == "HuBac", ]
      huBacData$Abundance <- huBacData$Abundance * 8
      ruBacData <- graphData2[graphData2$Bac == "RuBac", ]
      ruBacData$Abundance <- ruBacData$Abundance * 8
      
      mgPlot <- mgPlot + geom_point(data = huBacData, 
                                    aes(y = Abundance), colour = "Red", size = 4)
      mgPlot <- mgPlot + geom_point(data = ruBacData, 
                                    aes(y = Abundance), colour = "Blue", size = 4)
      
      if (is.null(yAxisVals))
      {
        mgPlot <- mgPlot + 
          scale_y_continuous(sec.axis = sec_axis(~ . / 8, name = "HuBac / RuBac reads per scg"))
      }else
      {
        mgPlot <- mgPlot + 
          scale_y_continuous(limits = yAxisVals, 
                             sec.axis = sec_axis(~ . / 8, name = "HuBac / RuBac reads per scg"))
      }
      
    }else
      
      if (!is.null(yAxisVals))
      {
        mgPlot <- mgPlot + ylim(yAxisVals)
      }
    
  }

  mgPlot <- mgPlot +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11))

  mgPlot <- mgPlot +
    theme(panel.background = element_rect(fill = "white", colour = "grey",
                                          linewidth = 0.5, linetype = "solid"))
  mgPlot <- mgPlot + theme(axis.title.y = element_text(size = 12))
  mgPlot <- mgPlot + theme(axis.title.x = element_blank())
  mgPlot <- mgPlot + theme(legend.title = element_text(size = 10),
                           legend.text = element_text(size = 9))
  mgPlot <- mgPlot + theme(axis.text.y = element_text(size = 10))

  mgPlot <- mgPlot + scale_x_discrete(labels = xAxisLabels)

  # mgPlot <- mgPlot + labs(title = plotTitle)

  mgPlot

  ggsave(paste(plotsDir, plotTitle, ".pdf", sep = ""), mgPlot, width = 10, height = 6)

}
