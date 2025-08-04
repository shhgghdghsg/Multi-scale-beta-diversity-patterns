library(tidyverse)
library(vegan)
library(cowplot)

load("clean.data/westsib.rda")

# The function to compute beta diversity (1-CqN) among three nested levels
# q takes values from 0 to 3 by 0.1

hier.CqN <- function(cdm, meta) {
  
  #---among ecosystems inside subzones---#
  meta$veg <- factor(meta$veg)
  
  agg <- matrix(nrow = nlevels(meta$veg), ncol = ncol(cdm))
  for (ii in 1:nrow(agg)) agg[ii, ] <- colSums(cdm[meta$veg == levels(meta$veg)[ii], ])
  
  meta_agg <- data.frame(veg = levels(meta$veg))
  
  q <- seq(0, 3, by = 0.1)
  N <- nrow(agg)
  da <- renyi(agg, scales = q, hill = T) |> colMeans()
  dg <- renyi(colSums(agg), scales = q, hill = T)
  db <- dg/da
  ha <- diversity(agg) |> exp() |> mean() |> log()
  hg <- diversity(colSums(agg))
  
  CqN.e <- 1 - (db^(1-q) - N^(1-q)) / (1 - N^(1-q)) # Chao et al. 2012: formula 12c
  CqN.e[q == 1] <- (hg - ha)/log(N)         # Chao et al. 2012: formula 13b

  
  meta$id <- paste0(meta$veg,"/", meta$top)
  meta$id <- factor(meta$id, levels = unique(meta$id))
  
  agg <- matrix(nrow = nlevels(meta$id), ncol = ncol(cdm))
  for (ii in 1:nrow(agg)) agg[ii, ] <- colSums(cdm[meta$id == levels(meta$id)[ii], ])
  
  meta_agg <- data.frame(id = levels(meta$id), 
                         sys = sapply(strsplit(levels(meta$id), "/"), function(x) x[1]) |> factor(),
                         bio = sapply(strsplit(levels(meta$id), "/"), function(x) x[2]))
  
  #---among samples inside microhabitats---#
  CqN.s <- matrix(nrow = nlevels(meta$id), ncol = length(q))
  
  for (ii in 1:nlevels(meta$id)) {
    N <- nrow(cdm[meta$id == levels(meta$id)[ii], ])
    da <- renyi(cdm[meta$id == levels(meta$id)[ii], ], scales = q, hill = T) |> colMeans()
    dg <- renyi(colSums(cdm[meta$id == levels(meta$id)[ii], ]), scales = q, hill = T)
    db <- dg/da
    ha <- diversity(cdm[meta$id == levels(meta$id)[ii], ]) |> exp() |> mean() |> log()
    hg <- diversity(colSums(cdm[meta$id == levels(meta$id)[ii], ]))
    
    CqN.s[ii,] <- 1 - (db^(1-q) - N^(1-q)) / (1 - N^(1-q)) # Chao et al. 2012: formula 12c
    CqN.s[ii,][q == 1] <- (hg - ha)/log(N)         # Chao et al. 2012: formula 13b
  }
  
  CqN.s <- colMeans(CqN.s)
  
  #---among microhabitats inside ecosystems---#
  CqN.b <- matrix(nrow = nlevels(meta_agg$sys), ncol = length(q))
  
  for (ii in 1:nlevels(meta_agg$sys)) {
    N <- nrow(agg[meta_agg$sys == levels(meta_agg$sys)[ii], ])
    da <- renyi(agg[meta_agg$sys == levels(meta_agg$sys)[ii], ], scales = q, hill = T) |> colMeans()
    dg <- renyi(colSums(agg[meta_agg$sys == levels(meta_agg$sys)[ii], ]), scales = q, hill = T)
    db <- dg/da
    ha <- diversity(agg[meta_agg$sys == levels(meta_agg$sys)[ii], ]) |> exp() |> mean() |> log()
    hg <- diversity(colSums(agg[meta_agg$sys == levels(meta_agg$sys)[ii], ]))
    
    CqN.b[ii,] <- 1 - (db^(1-q) - N^(1-q)) / (1 - N^(1-q)) # Chao et al. 2012: formula 12c
    CqN.b[ii,][q == 1] <- (hg - ha)/log(N)         # Chao et al. 2012: formula 13b
  }
  
  CqN.b <- colMeans(CqN.b)
  
  CqN <- c(CqN.s, CqN.b, CqN.e)
  names(CqN) <- c(paste("sam", seq(0, 3, by = 0.1)), paste("top", seq(0, 3, by = 0.1)), paste("veg", seq(0, 3, by = 0.1)))
  CqN
}

# Calculate beta diversity in 4 sub-zones
# nftu: northern forest-tundra, ftu: forest-tundra", tai: taiga, stai: sub-taiga
beta_nftu <- oecosimu(cdm[meta$zone == "ftu-n", ], hier.CqN, "quasiswap_count", nsimul = 9999, meta = meta[meta$zone == "ftu-n", ])
beta_ftu <- oecosimu(cdm[meta$zone == "ftu", ],    hier.CqN, "quasiswap_count", nsimul = 9999, meta = meta[meta$zone == "ftu", ])
beta_tai <- oecosimu(cdm[meta$zone == "tai", ],    hier.CqN, "quasiswap_count", nsimul = 9999, meta = meta[meta$zone == "tai", ])
beta_stai <- oecosimu(cdm[meta$zone == "stai", ],  hier.CqN, "quasiswap_count", nsimul = 9999, meta = meta[meta$zone == "stai", ])

# The contents of Table 2 are extracted from the above data frames


# Plot format
plot.theme = theme(plot.title = element_text(size = 20, color = "black", family = "sans", face = "bold", vjust = 0.5, hjust = 0.5),
                   axis.line = element_line(size = 0.5, colour = "black"),
                   axis.ticks = element_line(color = "black"),
                   axis.text = element_text(size = 18, color = "black", family = "sans", vjust = 0.5, hjust = 0.5),
                   axis.title = element_text(size = 20, color = "black", family = "sans",  face= "bold", vjust = 0.5, hjust = 0.5),
                   legend.position = "none",
                   legend.background = element_blank(),
                   legend.key = element_blank(),
                   legend.text = element_text(colour = 'black', size = 20, family = "sans", face = 'bold'),
                   legend.title = element_blank(),
                   panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank() ) 


#-------------------------------------- Figure 3 -------------------------------------#

#------ northern forest-tundra ------#

# Calculate beta diversity among three nested levels in the northern forest-tundra
CqN_nftu <- data.frame(beta_nftu$statistic,beta_nftu$oecosimu$pval)

# Abbreviation of the column names
names(CqN_nftu) <- c("obs","pval")

# Add the full names of three nested levels for easier plotting
CqN_nftu$level <- c(rep("replicates", 31), rep("microhabitats", 31), rep("ecosystems", 31))
CqN_nftu$level <- factor(CqN_nftu$level, levels = c("ecosystems", "microhabitats", "replicates"))

# Convert to significance marks by p-value
CqN_nftu$sig <- ifelse(CqN_nftu$pval < 0.05, "*", "")
CqN_nftu$sig <- factor(CqN_nftu$sig, levels = c("*", ""))

# Add the order q for easier plotting
CqN_nftu$q <- rep(seq(0, 3, by = 0.1), 3)

# Fig 3a. Multiplicative partitioning of testate amoeba diversity in the northern forest-tundra zone
# Lines are diversity profiles for three scale levels 
# Points marked with an asterisk represent the significant cases compared to expected results in the null model analysis
p3a <- ggplot(CqN_nftu, aes(x = q, y = obs, color = level)) +
  geom_point() +
  geom_line() +
  ylim(0,0.7) +
  labs(x = "q", y = "1 - CqN") +
  geom_text(data = CqN_nftu, aes(x = q, y = obs + 0.005), label = CqN_nftu$sig, color = "black", size = 5) +
  plot.theme +
  theme(legend.position= c(0.22, 0.9))


# The following three sub-zones are plotted in a similar manner to Fig 3a (northern forest-tundra zone)

#------ forest-tundra zone ------#

CqN_ftu <- data.frame(beta_ftu$statistic,beta_ftu$oecosimu$pval)

names(CqN_ftu) <- c("obs", "pval")

CqN_ftu$level <- c(rep("replicates", 31), rep("microhabitats", 31), rep("ecosystems", 31))
CqN_ftu$level <- factor(CqN_ftu$level, levels = c("ecosystems", "microhabitats", "replicates"))

CqN_ftu$sig <- ifelse(CqN_ftu$pval < 0.05, "*", "")
CqN_ftu$sig <- factor(CqN_ftu$sig, levels = c("*", ""))

CqN_ftu$q <- rep(seq(0, 3, by = 0.1), 3)

# Fig 3b. Multiplicative partitioning of testate amoeba diversity in the forest-tundra zone
p3b <- ggplot(CqN_ftu, aes(x = q, y = obs, color = level)) +
  geom_point() +
  geom_line() +
  ylim(0,0.7) +
  labs(x = "q", y = "1 - CqN") +
  geom_text(data = CqN_ftu, aes(x = q, y = obs + 0.005), label = CqN_ftu$sig, color = "black", size = 5) +
  plot.theme 

#---- taiga zone ----#

CqN_tai <- data.frame(beta_tai$statistic, beta_tai$oecosimu$pval)

names(CqN_tai) <- c("obs", "pval")

CqN_tai$level <- c(rep("replicates", 31), rep("microhabitats", 31), rep("ecosystems", 31))
CqN_tai$level <- factor(CqN_tai$level, levels = c("ecosystems", "microhabitats", "replicates"))

CqN_tai$sig <- ifelse(CqN_tai$pval < 0.05, "*", "")
CqN_tai$sig <- factor(CqN_tai$sig, levels = c("*", ""))

CqN_tai$q <- rep(seq(0, 3, by = 0.1), 3)

# Fig 3c. Multiplicative partitioning of testate amoeba diversity in the taiga zone
p3c <- ggplot(CqN_tai, aes(x = q, y = obs, color = level)) +
  geom_point() +
  geom_line() +
  ylim(0,0.7) +
  labs(x = "q", y = "1 - CqN") +
  geom_text(data = CqN_tai, aes(x = q, y = obs + 0.005), label = CqN_tai$sig, color = "black", size = 5) +
  plot.theme 

#---- sub-taiga zone ----#

CqN_stai <- data.frame(beta_stai$statistic, beta_stai$oecosimu$pval)

names(CqN_stai) <- c("obs", "pval")

CqN_stai$level <- c(rep("replicates", 31), rep("microhabitats", 31), rep("ecosystems", 31))
CqN_stai$level <- factor(CqN_stai$level, levels = c("ecosystems", "microhabitats", "replicates"))

CqN_stai$sig <- ifelse(CqN_stai$pval < 0.05, "*", "")
CqN_stai$sig <- factor(CqN_stai$sig, levels = c("*", ""))

CqN_stai$q <- rep(seq(0, 3, by = 0.1), 3)


# Fig 3d. Multiplicative partitioning of testate amoeba diversity in the sub-taiga zone
p3d <- ggplot(CqN_stai, aes(x = q, y = obs, color = level)) +
  geom_point() +
  geom_line() +
  ylim(0,0.7) +
  labs(x = "q", y = "1 - CqN") +
  geom_text(data = CqN_stai, aes(x = q, y = obs + 0.005), label = CqN_stai$sig, color = "black", size = 5) +
  plot.theme 

# Fig 3
plot_grid(p3a, p3b, p3c, p3d, ncol = 2, nrow = 2)



#-------------------------------------- Figure 4 -------------------------------------#

# Extract observed beta diversity among ecosystems for four sub-zones
beta.veg.obs <- c(beta_nftu$statistic[63:93], beta_ftu$statistic[63:93],
                  beta_tai$statistic[63:93], beta_stai$statistic[63:93])

# Extract expected beta diversity among ecosystems for four sub-zones
beta.veg.exp <- c(beta_nftu$oecosimu$means[63:93], beta_ftu$oecosimu$means[63:93],
                  beta_tai$oecosimu$means[63:93], beta_stai$oecosimu$means[63:93])

# Extract p values of beta diversity among ecosystems for four sub-zones
beta.veg.pval <- c(beta_nftu$oecosimu$pval[63:93], beta_ftu$oecosimu$pval[63:93],
                   beta_tai$oecosimu$pval[63:93], beta_stai$oecosimu$pval[63:93])


# Merge all data on beta diversity among ecosystems
beta.veg <- data.frame(beta.veg.obs, beta.veg.exp, beta.veg.pval)

# Convert to significance marks by p-value
beta.veg$sig <- ifelse(beta.veg$beta.veg.pval < 0.05, "*", "")
beta.veg$sig <- factor(beta.veg$sig, levels = c("*", ""))

# Add the full names of the four sub-zones for easier plotting
beta.veg$zone <- c(rep("northern forest-tundra", 31), rep("forest-tundra", 31), rep("taiga", 31), rep("sub-taiga", 31))
beta.veg$zone <- factor(beta.veg$zone, levels = c("northern forest-tundra", "forest-tundra", "taiga", "sub-taiga"))

# Add the order q for easier plotting
beta.veg$q <- rep(seq(0, 3, by = 0.1), 4)


# Fig 4. Among-ecosystem β-diversity of testate amoeba communities in Western Siberia. 
# Lines are diversity profiles for the sub-zones. 
# Points marked with an asterisk represent the significant cases as compared to expected results in the null model analysis.
ggplot(beta.veg, aes(x = q, y = beta.veg.obs, color = zone)) +
  geom_point() +
  geom_line() +
  ylim(0, 0.7) +
  labs(x = "q", y = "1 - CqN") +
  geom_text(data = beta.veg, aes(x = q, y = beta.veg.obs + 0.005), label = beta.veg$sig, color = "black", size = 5) +
  plot.theme +
  theme(legend.position= c(0.2, 0.85))
