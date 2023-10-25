pur_pl <-read_tsv(paste0(HOME, patient, "/", patient, "_purple/", sample, ".purple.purity.range.tsv"), show_col_types = FALSE)
pur_pl <- pur_pl[order(pur_pl$ploidy, pur_pl$purity),]
ploidies <- as.data.frame(unique(pur_pl$ploidy))
colnames(ploidies) <- "ploidy"
ploidies <- ploidies %>% 
  mutate(diff = abs(ploidy - lead(ploidy)))
ploidies[is.na(ploidies$diff), 2] <- max(ploidies$diff, na.rm=TRUE)

for (i in 1:nrow(pur_pl)) {
  diff <- ploidies[ploidies$ploidy == pur_pl[[i, 5]], 2]
  pur_pl[i, 7] <- diff
}

p <- ggplot(pur_pl, aes(ploidy, purity)) +
  geom_tile(aes(fill = score, width = diff),
            data = pur_pl,
            show.legend = FALSE) +
  scale_fill_gradientn(trans = "log",
                       colours = rev(RColorBrewer::brewer.pal(10, "RdBu"))) +
  scale_x_continuous(breaks = 1:8) +
  xlab("psi (Biased ploidy)") +
  ylab("rho (Aberrant cell fraction)")

ggsave(paste0(sample, "_sunrise.png"), path = publishDir, plot = p, width = 40, height = 25, units = "cm")