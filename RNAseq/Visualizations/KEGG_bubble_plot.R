install.packages("ggplot2")

data2 <- read.csv("zyg_kegg.csv")

head(data2)

ggplot(data2, aes(x = Rich_factor, y = reorder(Pathway_term, -Rich_factor))) +
  geom_point(aes(size = Gene_number, color = P_Value), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", guide = "colour", name = "P-Value") +
  scale_size_continuous(range = c(2, 12), guide = "none") +
  labs(
    x = "Rich Factor", y = "Pathway Term",
    title = "Bubble Plot of Pathway Enrichment Analysis (New Data)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
