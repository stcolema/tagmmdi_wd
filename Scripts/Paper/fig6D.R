require(dplyr)
require(ggplot2)

col_pal <- c(
  "#a7bed3",
  "#adf7b6",
  "#c6e2e9",
  "#f1ffc4",
  "#ffcaaf",
  "#dab894",
  "#70d6ff",
  "#ff70a6",
  "#ff9770",
  "#ffd670",
  "#e9ff70",
  "#f08080",
  "#fbc4ab",
  "#dfb2f4",
  "#f5e960",
  "#f5e5f0",
  "#55d6c2",
  "#ffffff",
  "#84dcc6",
  "#a5ffd6",
  "#79addc",
  "#ffc09f",
  "#ffee93",
  "#fcf5c7",
  "#ffa69e",
  "#ff686b",
  "#808080"
)


mdi_tagm_alloc <- read.csv(file = "~/Desktop/TGOndiiMDIandTAGM.csv")


marker_levels <- c(
  "apical 1",
  "apical 2",
  "micronemes",
  "rhoptries 1",
  "rhoptries 2",
  "dense granules",
  "IMC",
  "tubulin cytoskeleton",
  "PM - peripheral 1",
  "endomembrane vesicles",
  "PM - peripheral 2",
  "PM - integral",
  "Golgi",
  "ER",
  "ER 2",
  "apicoplast",
  "mitochondrion - membranes",
  "mitochondrion - soluble",
  "nucleus - chromatin",
  "nucleus - non-chromatin",
  "nucleolus",
  "40S ribosome",
  "60S ribosome",
  "cytosol",
  "19S proteasome",
  "20S proteasome",
  "unknown"
)
marker_labels <- c(
  "apical 1",
  "apical 2",
  "micronemes",
  "rhoptries 1",
  "rhoptries 2",
  "dense granules",
  "IMC",
  "tubulin cytoskeleton",
  "PM - peripheral 1",
  "endomembrane vesicles",
  "PM - peripheral 2",
  "PM - integral",
  "Golgi",
  "ER 1",
  "ER 2",
  "apicoplast",
  "mitochondrion - membranes",
  "mitochondrion - soluble",
  "nucleus - chromatin",
  "nucleus - non-chromatin",
  "nucleolus",
  "40S ribosome",
  "60S ribosome",
  "cytosol",
  "19S proteasome",
  "20S proteasome",
  "all other proteins"
)

names(col_pal) <- marker_labels[-27] # c(names(pals::alphabet()), "grey50")

cohere <- mdi_tagm_alloc %>% filter(mdi.mcmc.allocation %in%
                                      c("dense granules",
                                        "rhoptries 1",
                                        "rhoptries 2",
                                        "apical 1",
                                        "apical 2",
                                        "micronemes")) %>%
  mutate(samealloc = tagm.mcmc.allocation == mdi.mcmc.allocation)

p1 <- cohere %>%
  ggplot(aes(x = tagm.mcmc.probability,
             y = mdi.mcmc.probability,
             fill = mdi.mcmc.allocation,
             shape = !samealloc)
         ) +
  theme_bw() + 
  geom_point(size = 2, alpha = 0.9, color = "#808080") +
  theme(text = element_text(size = 20)) +
  ylim(c(0.3,1)) + 
  xlim(c(0.3, 1)) +
  coord_fixed() + 
  geom_abline(slope=1, intercept = 0, linewidth = 1) +
  labs(shape = "Different\nclassification", 
       x = "TAGM allocation probability",
       y = "MDI allocation probability",
       fill = "MDI predicted\nlocalisation") +
  scale_fill_manual(values = col_pal) +
  scale_shape_manual(values = c(21, 24)) +
  guides(fill = guide_legend(override.aes = list(shape=21))) +
  theme(legend.position="bottom")

ggsave("./Plots/Fig6/CompareTAGMandMDIallocationFig6D.pdf", 
       device = "pdf",
       width = 7,
       height = 8.5
  )
