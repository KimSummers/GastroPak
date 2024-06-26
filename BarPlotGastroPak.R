# BarPlotGastroPak
#
# Bar chart plot GastroPak data
#
# file BarPlotGastroPak
#
# inputs
# 	graphData         - file containing plate count data
#   xAxisLabels       - Labels for the x-axis
#   plotType          - "Sample" to fill by sample, "Species" to fill by Bacteria,
#                       "Media" to fill by media
#   dataType          - "Env" for environmental, "HH" for household
#   mainData          - Type of media if "Sample" or "Species", type of bacteria if "media"
#   subData           - Type of bacteria if by "Sample" or sample if by "Species"
#   fillColours       - Colours to fill by
#   weightOrVol       - by weight or by volume or both
#   rowGraphs         - number of facets in each row
#   plotsDir          - directory to store plots in
#   plotTitle         - title for the plot
#
# Version    Author       Date      Affiliation
# 1.00       J K Summers  20/12/23  Wellington Lab - School of Life Sciences - University of Warwick
BarPlotGastroPak <- function(graphData, xAxisLabels = NULL, plotType, dataType, mainData,
                             subData, fillColours, weightOrVol, rowGraphs,
                             plotsDir, plotTitle) {

  library(scales)
  library(ggplot2)
  library(pals)

  graphData$MeanCfu[graphData$MeanCfu == 0] <- NaN

  if (plotType == "Sample")
  {

    if (dataType == "HH")
    {
      gpPlot <- ggplot(data = graphData, aes(x = SampleID, y = MeanCfu,
                                             fill = SampleType))
    }else
    {
      gpPlot <- ggplot(data = graphData, aes(x = SamplingSite, y = MeanCfu,
                                           fill = SampleType))
    }

  }else
  {

    if (plotType == "Species")
    {

      if (dataType == "HH")
      {
        gpPlot <- ggplot(data = graphData, aes(x = SampleID, y = MeanCfu,
                                               fill = Bacteria))
      }else
      {
        gpPlot <- ggplot(data = graphData, aes(x = SamplingSite, y = MeanCfu,
                                               fill = Bacteria))
      }

    }else
    {

      if (dataType == "HH")
      {
        gpPlot <- ggplot(data = graphData, aes(x = SampleID, y = MeanCfu,
                                               fill = Media))
      }else
      {
        gpPlot <- ggplot(data = graphData, aes(x = SamplingSite, y = MeanCfu,
                                               fill = Media))
      }

    }

  }

  gpPlot <- gpPlot + geom_bar(stat = "identity", position = "dodge")

  gpPlot <- gpPlot + geom_errorbar(aes(ymin = MeanCfu - StdDev,
                                       ymax = MeanCfu + StdDev),
                                   width = 0.2, position = position_dodge(.9))
  gpPlot <- gpPlot +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8))

  if (plotType == "Species")
  {
    fillLabel <- "Bacteria"
  }else
  {

    if (plotType == "Sample")
    {
      fillLabel <- "Sample Type"
    }else
    {
      fillLabel <- "Media"
    }

  }

  if (weightOrVol == "W")
  {
    gpPlot <- gpPlot + labs(y = bquote('Log mean cfu 0.1'~g^-1),
                            title = paste(mainData, subData), fill = fillLabel)
  }else
  {

    if (weightOrVol == "V")
    {
      gpPlot <- gpPlot + labs(y = bquote('Log mean cfu'~ml^-1),
                              title = paste(mainData, subData), fill = fillLabel)
    }else
    {

      if (dataType == "Env")
      {
        gpPlot <- gpPlot + labs(y = "Log mean cfu / ml (water) or cfu / 0.1 g (sediment)",
                                title = paste(mainData, subData), fill = fillLabel)
      }else
      {
        gpPlot <- gpPlot + labs(y = "Log mean cfu / ml (water) or cfu / 0.1 g (faeces)",
                                title = paste(mainData, subData), fill = fillLabel)
      }

    }

  }

  if (dataType == "Env")
  {
    gpPlot <- gpPlot + facet_wrap( ~ River, ncol = rowGraphs, scales = "free_x")
  }else
  {
    gpPlot <- gpPlot + facet_wrap( ~ Household, ncol = rowGraphs, scales = "free_x")
  }

  gpPlot <- gpPlot +
    theme(strip.background = element_rect(fill = "white", colour = "grey"),
          strip.text = element_text(size = 8, colour = "black"))

  gpPlot <- gpPlot +
    theme(panel.background = element_rect(fill = "white", colour = "grey",
                                          linewidth = 0.5, linetype = "solid"))

  gpPlot <- gpPlot + theme(axis.title.y = element_text(size = 12))
  gpPlot <- gpPlot + theme(axis.title.x = element_blank())
  gpPlot <- gpPlot + theme(axis.text = element_text(size = 8))

  gpPlot <- gpPlot + theme(plot.title = element_text(size = 14, hjust = 0.5))
  gpPlot <- gpPlot + theme(legend.title = element_text(size = 10),
                           legend.text = element_text(size = 8))

  gpPlot <- gpPlot + scale_fill_manual(values = fillColours)

  gpPlot <- gpPlot + scale_y_continuous(trans = 'log10')

  if (!is.null(xAxisLabels))
  {
    gpPlot <- gpPlot + scale_x_discrete(labels = xAxisLabels)
  }

  gpPlot

  if (rowGraphs == 3)
  {
    graphWidth <- 6
  }else

    if (rowGraphs == 4)
    {
      graphWidth <- 8
    }else
    {
      graphWidth <- 12
    }

  ggsave(paste(plotsDir, plotTitle, ".pdf", sep = ""), gpPlot,
         width = graphWidth, height = 6)

}
