#### Functions from hmftools repository
purity_ploidy_range_plot <- function(bestFit, range) {
  
  bestPurity = bestFit[1, "purity"]
  bestPloidy = bestFit[1, "ploidy"]
  bestScore = bestFit[1, "score"]
  
  range =  range %>%
    arrange(purity, ploidy) %>%
    group_by(purity) %>%
    mutate(
      absScore = pmin(4, score),
      score = pmin(1, abs(score - bestScore) / score),
      leftPloidy = lag(ploidy),
      rightPloidy = lead(ploidy),
      xmin = ploidy - (ploidy - leftPloidy) / 2,
      xmax = ploidy + (rightPloidy - ploidy) / 2,
      ymin = purity - 0.005,
      ymax = purity + 0.005,
      xmin = ifelse(is.na(xmin), ploidy, xmin),
      xmax = ifelse(is.na(xmax), ploidy, xmax))
  
  maxPloidy = min(range %>% arrange(purity, -ploidy) %>% group_by(purity)  %>% filter(row_number() == 1) %>% dplyr::select(purity, ploidy = xmax) %>% ungroup() %>% dplyr::select(ploidy))
  minPloidy = max(range %>% arrange(purity, ploidy) %>% group_by(purity)  %>% filter(row_number() == 1) %>% dplyr::select(purity, maxPloidy = xmin) %>% ungroup() %>% dplyr::select(maxPloidy))
  
  maxPloidy = max(maxPloidy, bestPloidy)
  minPloidy = min(minPloidy, bestPloidy)
  
  range = range %>%
    filter(xmin <= maxPloidy, xmax >= minPloidy) %>%
    mutate(xmax = pmin(xmax, maxPloidy), xmin = pmax(xmin, minPloidy))
  
  result = ggplot(range) +
    geom_rect(aes(fill=score, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
    scale_fill_gradientn(colours=c("blue","blue", "green", "yellow","orange", "red", "red2"), limits = c(0, 1), values=c(0, 0.0999, 0.1, 0.5, 0.8, 0.9, 1), breaks = c(0.1, 0.5, 1), labels = c("10%", "50%", "100%"), name = "Relative\nScore") +
    geom_segment(aes(y = 0.085, yend = 1.05, x=bestPloidy, xend = bestPloidy), linetype = "dashed", size = 0.1) +
    geom_label(data = data.frame(), aes(x = bestPloidy, y = 1.05, label = round(bestPloidy, 2)), size = 2.5) +
    geom_segment(aes(y = bestPurity, yend = bestPurity, x=minPloidy, xend = maxPloidy + 0.4), linetype = "dashed", size = 0.1) +
    geom_label(data = data.frame(), aes(y = bestPurity, x = maxPloidy + 0.4, label = paste0(bestPurity*100,"%" )), size = 2.5, hjust = 0.7) +
    theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title=element_text(size=6), legend.text=element_text(size=6)) +
    scale_y_continuous(labels = c("25%", "50%", "75%", "100%"), breaks = c(0.25, 0.5, 0.75, 1)) +
    xlab("Ploidy") + ylab("Purity") +  ggtitle("Purity/Ploidy Scores")
  
  return (result)
}

clonality_plot <- function(somaticBuckets, clonalityModel) {
  clonalityVariants = somaticBuckets %>% group_by(variantCopyNumberBucket) %>% summarise(count = sum(count))
  
  subclonalPercentage = clonalityModel %>% 
    group_by(bucket) %>% 
    mutate(totalWeight = sum(bucketWeight, na.rm = T)) %>%
    filter(isSubclonal) %>%
    summarise(
      isSubclonal = T,
      bucketWeight = sum(bucketWeight, na.rm = T), 
      subclonalLikelihood = ifelse(bucketWeight == 0 & is.nan(bucketWeight), 0, bucketWeight / max(totalWeight)))
  
  nonResidualModel = clonalityModel %>% filter(peak != 0)
  
  nonResidualSubclonalPercentage = nonResidualModel %>%
    group_by(bucket) %>%
    mutate(totalWeight = sum(bucketWeight)) %>%
    filter(isSubclonal) %>%
    summarise(
      isSubclonal = T,
      bucketWeight = sum(bucketWeight),
      subclonalLikelihood = ifelse(bucketWeight == 0, 0, bucketWeight / max(totalWeight)))
  
  combinedModel = nonResidualModel %>%
    group_by(bucket) %>% 
    summarise(bucketWeight = sum(bucketWeight))
  
  singleBlue = "#6baed6"
  singleRed = "#d94701"
  
  pTop = ggplot() +
    geom_bar(data=clonalityVariants, aes(x = variantCopyNumberBucket, weight = count), fill=singleBlue, col=singleBlue,  alpha = .4, size = 0.07, width = 0.05) +
    geom_line(data=combinedModel , aes(x = bucket, y = bucketWeight), position = "identity", alpha = 0.8) +
    geom_line(data=nonResidualModel, aes(x = bucket, y = bucketWeight, color = peak), position = "identity") +
    geom_area(data=nonResidualSubclonalPercentage %>% filter(isSubclonal), aes(x = bucket, y = bucketWeight), position = "identity",  alpha = 0.3, fill = singleRed, color = singleRed) +
    ggtitle("") + xlab("Variant Copy Number") + ylab("") +
    scale_y_continuous(expand=c(0.02, 0.02)) +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position="none") +
    scale_x_continuous( expand=c(0.01, 0.01), limits = c(0, 3.5)) 
  
  pBottom = ggplot(data = subclonalPercentage) +
    geom_bar(width = 0.05, aes(x = bucket, y = subclonalLikelihood), stat = "identity", fill=singleRed, col=singleRed,  alpha = 0.3) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank()) +
    xlab("") + ylab("") +
    #scale_y_continuous(labels = c("0%", "25%","50%","75%","100%"), breaks = c(0, 0.25, 0.5, 0.75, 1), expand=c(0.02, 0.02), limits = c(0, 1))# +
    scale_x_continuous( expand=c(0.01, 0.01), limits = c(0, 3.5)) 
  
  pFinal = cowplot::plot_grid(pTop, pBottom, ncol = 1, rel_heights = c(5, 1), align = "v")
  return(pFinal)
}
