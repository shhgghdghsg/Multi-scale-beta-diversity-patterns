library(tidyverse)
library(vegan)
library(pairwiseAdonis)

load("clean.data/westsib.rda")

# PERMANOVA results
panova1 <- adonis2(cdm ~ zone, meta, method = "bray", permutations = 9999) # sub-zones
panova1

panova2 <- adonis2(cdm ~ veg, meta, method = "bray", permutations = 9999) # ecosystems
panova2

panova3 <- adonis2(cdm ~ top, meta, method = "bray", permutations = 9999) # microhabitats
panova3


# Hellinger transformation of species abundances
spe.h <- decostand(cdm, "hellinger")

# PCA
spe.h.pca <- rda(spe.h)

# Variance explained by the first two axes
axis1 <- round(spe.h.pca$CA$eig[1]/sum(spe.h.pca$CA$eig)*100, 2)
axis2 <- round(spe.h.pca$CA$eig[2]/sum(spe.h.pca$CA$eig)*100, 2)
axis <- axis1 + axis2

sp.sc <- scores(spe.h.pca, display = "species") 
st.sc <- scores(spe.h.pca, display = "sites")

st.sc <- cbind(as.data.frame(st.sc), zone = meta$zone)

cent <- aggregate(cbind(PC1, PC2) ~ zone, data = st.sc, FUN = mean)
segs <- merge(st.sc, setNames(cent, c('zone','oPC1','oPC2')),
              by = 'zone', sort = FALSE)

# ten species with the greatest contribution to the results of ordination
len <- sqrt(rowSums(sp.sc^2))
top <- order(len, decreasing = T)[1:10]
sp.sc.top <- data.frame(sp.sc[top, ])


# Create a dataframe to extract the convex hull points
pca_hull <- 
  st.sc %>% 
  group_by(zone) %>% 
  slice(chull(PC1, PC2))


# Rename the variables so they look nicer on the figure
name <- row.names(sp.sc.top)

name <-  recode(name,
                'Corythion.dubium' = 'Co.du',
                'Assulina.muscorum' = 'As.mu',
                'Centropyxis.sylvatica.minor' = 'Ce.mi',
                'Trinema.lineare' = 'Tr.li',
                'Trinema.enchelys'='Tr.en',
                'Cyclopyxis.kahli' = 'Cy.ka',
                'Assulina.seminulum' = 'As.se',
                'Centropyxis.aerophila' = 'Ce.ae',
                'Centropyxis.sylvatica' = 'Ce.sy',
                'Corythion.orbicularis' = 'Co.or')


# Drawing format
plot.theme = theme(plot.title=element_text(size=20, color="black", family  = "sans", face= "bold",vjust=0.5,hjust=0.5),
                   axis.line=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text= element_text(size=18, color="black", family  = "sans", vjust=0.5, hjust=0.5),
                   axis.title = element_text(size=20, color="black", family  = "sans",face= "bold", vjust=0.5, hjust=0.5),
                   legend.position=c(0.12,0.15),
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text = element_text(colour = 'black', size = 20,  family  = "sans",face = 'bold'),
                   legend.title=element_blank(),
                   panel.background=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank()  ) 

# Adjust the position of the label (species name) in the figure
sp.sc.top$vjust <- c(-1, -1, rep(0, 7),-1)
sp.sc.top$hjust <- c(0, 0, -.5, -.5, -.5, -.5, -.6, 1.1, -.5, .5)


# Fig. 2. Ordination biplot for the principal component analysis (PCA) based on Hellinger-transformed relative abundances of testate amoebae. 
# Arrows represent ten species with the greatest contribution to the results of ordination. 
# (As.mu - Assulina muscorum, As.se - Assulina seminulum, Ce.ae - Centropyxis aerophila, Ce.mi - Centropyxis sylvatica minor, Ce.sy - Centropyxis sylvatica, 
# Co.du - Corythion dubium, Co.or - Corythion orbicularis, Cy.ka - Cyclopyxis kahli, Tr.en - Trinema enchelys, Tr.li - Trinema lineare)
ggplot(st.sc, aes(x = PC1, y = PC2, colour = zone)) +
  geom_segment(data = segs,
               mapping = aes(xend = oPC1, yend = oPC2)) + # spiders
  geom_point(data = cent, size = 5) +  # centroids
  geom_point() +
  geom_polygon(data = pca_hull,
               aes(fill = zone,colour = zone),
               alpha = 0.3,
               show.legend = FALSE)+
  geom_hline(yintercept = 0,lty=2,col="black")+ 
  geom_vline(xintercept = 0,lty=2,col="black")+
  geom_segment(data=sp.sc.top,aes(x=0,y=0,xend=PC1,yend=PC2),
               arrow = arrow(angle = 22.5,length = unit(0.2,"cm"),type="closed"),
               linetype=1,size=0.6,color="black")+
  labs(x="PC1(28.41%)",
       y="PC2(9.18%)")+
  annotate('text', x = sp.sc.top$PC1, y = sp.sc.top$PC2,label = name, size = 6, 
           hjust = sp.sc.top$hjust, vjust= sp.sc.top$vjust)+
  plot.theme+
  guides(colour=guide_legend(override.aes=list(shape=15)))+
  scale_color_discrete(labels=c("northern forest-tundra","forest-tundra","taiga","sub-taiga"))




