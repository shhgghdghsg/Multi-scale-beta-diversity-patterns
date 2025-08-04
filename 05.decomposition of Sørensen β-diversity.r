library(tidyverse)
library(betapart)
library(vegan)
library(cowplot)

load("clean.data/westsib.rda")

# The function to compute Sorensen dissimilarity and its components on three nested levels

sor.decompose <- function(cdm, meta) {
  
  meta$veg <- factor(meta$veg)

  agg <- matrix(nrow = nlevels(meta$veg), ncol = ncol(cdm))
  for (ii in 1:nrow(agg)) agg[ii, ] <- colSums(cdm[meta$veg == levels(meta$veg)[ii], ])
  
  meta_agg <- data.frame(veg = levels(meta$veg))
  
  beta.e <- beta.pair((agg > 0) + 0) |> sapply(\(x) mean(x))

  meta$id <- paste0(meta$veg,"/", meta$top)
  meta$id <- factor(meta$id, levels = unique(meta$id))
  agg <- matrix(nrow = nlevels(meta$id), ncol = ncol(cdm))
  for (ii in 1:nrow(agg)) agg[ii, ] <- colSums(cdm[meta$id == levels(meta$id)[ii], ])
  meta_agg <- data.frame(id = levels(meta$id), 
                         sys = sapply(strsplit(levels(meta$id), "/"), function(x) x[1]) |> factor(),
                         bio = sapply(strsplit(levels(meta$id), "/"), function(x) x[2]))
  
  beta.s <- matrix(nrow = nlevels(meta$id), ncol = 3)
  for (ii in 1:nlevels(meta$id)) {
    beta.s[ii, ] <- unlist(beta.pair((cdm[meta$id == levels(meta$id)[ii], ] > 0) + 0)|> sapply(\(x) mean(x)))
  }
  beta.s <- colMeans(beta.s)
  
  beta.b <- matrix(nrow = nlevels(meta_agg$sys), ncol = 3)
  for (ii in 1:nlevels(meta_agg$sys)) {
    beta.b[ii, ] <- unlist(beta.pair((agg[meta_agg$sys == levels(meta_agg$sys)[ii], ] > 0) + 0)|> sapply(\(x) mean(x)))
  }
  beta.b <- colMeans(beta.b)
  
  beta <- c(beta.s, beta.b, beta.e)
  names(beta) <- c("bt.1", "bn.1", "bs.1", "bt.2", "bn.2", "bs.2", "bt.3", "bn.3", "bs.3")
  beta
}

beta_nftu <- oecosimu(cdm[meta$zone == "ftu-n", ], sor.decompose, "quasiswap_count", nsimul = 9999, meta = meta[meta$zone == "ftu-n", ])
beta_ftu <- oecosimu(cdm[meta$zone == "ftu", ],    sor.decompose, "quasiswap_count", nsimul = 9999, meta = meta[meta$zone == "ftu", ])
beta_tai <- oecosimu(cdm[meta$zone == "tai", ],    sor.decompose, "quasiswap_count", nsimul = 9999, meta = meta[meta$zone == "tai", ])
beta_stai <- oecosimu(cdm[meta$zone == "stai", ],  sor.decompose, "quasiswap_count", nsimul = 9999, meta = meta[meta$zone == "stai", ])


# Plot format
plot.theme = theme(plot.title = element_text(size = 20, color = "black", family = "sans", face = "bold", vjust = 0.5, hjust = 0.5),
                   axis.line = element_line(size = 0.5, colour = "black"),
                   axis.ticks = element_line(color = "black"),
                   axis.text = element_text(size = 18, color = "black", family = "sans", vjust = 0.5, hjust = 0.5),
                   axis.title = element_text(size = 20, color = "black", family = "sans", face = "bold", vjust = 0.5, hjust = 0.5),
                   legend.position = c(0.8, 0.9),
                   legend.background = element_blank(),
                   legend.key = element_blank(),
                   legend.text = element_text(colour = 'black', size = 20,  family = "sans",face = 'bold'),
                   legend.title = element_blank(),
                   panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank() )


# Fig 5. Decomposition of Sørensen β-diversity of testate amoeba assemblages into turnover and nestedness-resultant components 
# at ecosystem (a), microhabitat (b) and sample (c) levels of a hierarchical sampling design.
beta <- rbind(beta_nftu$statistic, beta_ftu$statistic,
              beta_tai$statistic, beta_stai$statistic)

# Fig 5a
eco <- data.frame(index = c(beta[, 7], beta[, 8])) |> mutate(
  component = rep(c("turnover", "nestedness"), each = 4),
  zone = factor(rep(c("northern forest-tundra", "forest-tundra", "taiga", "sub-taiga"), 2), levels = c("sub-taiga","taiga","forest-tundra","northern forest-tundra"))
)

p5a <- eco |> 
  ggplot(aes(x = zone, y = index, fill = component)) +
  geom_col(width = 0.5, alpha = 0.8)+
  coord_flip() +
  ylim(0, 0.7) +
  plot.theme +
  scale_fill_manual(values = c("#A9D179", "#AF46B4")) +
  labs(x = "", y = "")

# Fig 5b
bio <- data.frame(index = c(beta[, 4], beta[, 5])) |> mutate(
  component = rep(c("turnover", "nestedness"), each = 4),
  zone = factor(rep(c("northern forest-tundra", "forest-tundra", "taiga", "sub-taiga"), 2), levels = c("sub-taiga","taiga","forest-tundra","northern forest-tundra"))
)

p5b <- bio |>
  ggplot(aes(x = zone, y = index, fill = component)) +
  geom_col(width = 0.5, alpha = 0.8) +
  coord_flip() +
  ylim(0, 0.7) +
  plot.theme +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#A9D179", "#AF46B4")) +
  labs(x = "", y = "")

# Fig 5c
sam <- data.frame(index = c(beta[, 1], beta[, 2])) |> mutate(
  component = rep(c("turnover", "nestedness"), each = 4),
  zone = factor(rep(c("northern forest-tundra", "forest-tundra", "taiga", "sub-taiga"), 2), levels = c("sub-taiga","taiga","forest-tundra","northern forest-tundra"))
)

p5c <- sam |>
  ggplot(aes(x = zone, y = index, fill = component)) +
  geom_col(width = 0.5, alpha = 0.8) +
  coord_flip() +
  ylim(0, 0.7) +
  plot.theme +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#A9D179", "#AF46B4")) +
  labs(x = "", y = "")

# Fig 5
plot_grid(p5a, p5b, p5c, nrow = 3)



# Table 3. Null-model analysis of turnover and nestedness components of Sorensen dissimilarity (in the form of proportions)
# Null-model analysis of turnover and nestedness components for three levels: among sample replicates, among microhabitats, among ecosystems. 
# The observed partitions were compared with the expected values in the null model based on 9999 randomizations: (***) - P < 0.001, (**) - P < 0.01, (*) - P < 0.05, (ns) - P > 0.05.

# The function to compute Sorensen dissimilarity and proportions of its components on three nested levels

sor.decompose.prop <- function(cdm, meta) {
  
  meta$veg <- factor(meta$veg)
  
  agg <- matrix(nrow = nlevels(meta$veg), ncol = ncol(cdm))
  for (ii in 1:nrow(agg)) agg[ii, ] <- colSums(cdm[meta$veg == levels(meta$veg)[ii], ])
  
  meta_agg <- data.frame(veg = levels(meta$veg))
  
  beta.e <- beta.pair((agg > 0) + 0) |> sapply(\(x) mean(x))
  
  meta$id <- paste0(meta$veg,"/", meta$top)
  meta$id <- factor(meta$id, levels = unique(meta$id))
  agg <- matrix(nrow = nlevels(meta$id), ncol = ncol(cdm))
  for (ii in 1:nrow(agg)) agg[ii, ] <- colSums(cdm[meta$id == levels(meta$id)[ii], ])
  meta_agg <- data.frame(id = levels(meta$id), 
                         sys = sapply(strsplit(levels(meta$id), "/"), function(x) x[1]) |> factor(),
                         bio = sapply(strsplit(levels(meta$id), "/"), function(x) x[2]))
  
  beta.s <- matrix(nrow = nlevels(meta$id), ncol = 3)
  for (ii in 1:nlevels(meta$id)) {
    beta.s[ii, ] <- unlist(beta.pair((cdm[meta$id == levels(meta$id)[ii], ] > 0) + 0)|> sapply(\(x) mean(x)))
  }
  beta.s <- colMeans(beta.s)
  
  beta.b <- matrix(nrow = nlevels(meta_agg$sys), ncol = 3)
  for (ii in 1:nlevels(meta_agg$sys)) {
    beta.b[ii, ] <- unlist(beta.pair((agg[meta_agg$sys == levels(meta_agg$sys)[ii], ] > 0) + 0)|> sapply(\(x) mean(x)))
  }
  beta.b <- colMeans(beta.b)
  
  beta <- c(beta.s, beta.b, beta.e)

  beta.p <- c(beta[1] / beta[3], beta[2] / beta[3], beta[3], beta[4] / beta[6], beta[5] / beta[6], beta[6], beta[7] / beta[9], beta[8] / beta[9], beta[9])
  names(beta.p) <- c("btp.1", "bnp.1", "bs1", "btp.2", "bnp.2", "bs2", "btp.3", "bnp.3", "bs3")
  beta.p
}


beta_nftu <- oecosimu(cdm[meta$zone == "ftu-n", ], sor.decompose.prop, "quasiswap_count", nsimul = 9999, meta = meta[meta$zone == "ftu-n", ])
beta_ftu <- oecosimu(cdm[meta$zone == "ftu", ],    sor.decompose.prop, "quasiswap_count", nsimul = 9999, meta = meta[meta$zone == "ftu", ])
beta_tai <- oecosimu(cdm[meta$zone == "tai", ],    sor.decompose.prop, "quasiswap_count", nsimul = 9999, meta = meta[meta$zone == "tai", ])
beta_stai <- oecosimu(cdm[meta$zone == "stai", ],  sor.decompose.prop, "quasiswap_count", nsimul = 9999, meta = meta[meta$zone == "stai", ])

# The contents of Table 3 are extracted from the above data frames