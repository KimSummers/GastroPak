# BarPlotMetagenomicsSpecies
#
# Bar chart plot metagenomic data
#
# file BarPlotMetagenomics
#
# inputs
# 	graphData         - file containing plate count data
#   xAxisLabels       - Labels for the x-axis
#   absolute          - TRUE to plot absolute abundances, FALSE for percentage
#   absLab            - y label for absolute date (default is "Abundance per million reads")
#   huBac             - TRUE to plot a HuBac graph, FALSE for resistance genes
#   yAxisVals         - Optional parameter to override y-axis limits
#   plotsDir          - directory to store plots in
#   plotTitle         - title for the plot
#
# Version    Author       Date      Affiliation
# 1.00       J K Summers  14/12/23  Wellington Lab - School of Life Sciences - University of Warwick
BarPlotMetagenomicsSpecies <- function(graphData, xAxisLabels, absolute,
                                absLab = "Abundance per million reads",
                                huBac, yAxisVals = NULL, plotsDir, plotTitle) {

  library(pals)

  if (huBac)
  {
    mgPlot <- ggplot(data = graphData,
                     aes(x = factor(xVals, level = unique(xVals)), y = Hubac))

  }else
  {
    mgPlot <- ggplot(data = graphData,
                     aes(x = factor(xVals, level = unique(xVals)), y = Abundance))
  }

  if (huBac)
  {
    mgPlot <- mgPlot + geom_col()
    mgPlot <- mgPlot + labs(y = "Hubac per million reads")
  }else
  {

    if (absolute)
    {
      mgPlot <- mgPlot + geom_col(aes(fill = Family))
      mgPlot <- mgPlot + labs(y = absLab)
    }else
    {
      mgPlot <- mgPlot + geom_col(position = "fill",  aes(fill = Family))
      mgPlot <- mgPlot + labs(y = "Proportion of abundance")
    }

    mgPlot <- mgPlot + guides(color = guide_legend(override.aes = list(size = 10)))
    mgPlot <- mgPlot + scale_fill_manual(values = glasbey(), drop = FALSE)

  }

  mgPlot <- mgPlot +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1, size = 10))

  mgPlot <- mgPlot +
    theme(panel.background = element_rect(fill = "white", colour = "grey",
                                          linewidth = 0.5, linetype = "solid"))
  mgPlot <- mgPlot + theme(axis.title.y = element_text(size = 12))
  mgPlot <- mgPlot + theme(axis.title.x = element_blank())
  mgPlot <- mgPlot + theme(legend.title = element_text(size = 10),
                           legend.text = element_text(size = 8))
  mgPlot <- mgPlot + theme(axis.text.y = element_text(size = 8))

  mgPlot <- mgPlot + scale_x_discrete(labels = xAxisLabels)

  if (!is.null(yAxisVals))
  {
    mgPlot <- mgPlot + ylim(yAxisVals)
  }

  mgPlot <- mgPlot + labs(title = plotTitle)

  mgPlot

  ggsave(paste(plotsDir, plotTitle, ".pdf", sep = ""), mgPlot, width = 10, height = 6)

}
