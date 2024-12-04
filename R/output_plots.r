####################
# Plot grouped examples
####################

## *** Must run inla_glmm_birds first and output the required Rdata files.

library(dplyr)
library(ggplot2)
library(readxl)

setwd(here::here())

load("output/results.Rda")
load("output/results_mgmt.Rda")
load("output/temporalResults.Rda")

ggplot(data = time.results, aes(x = Year, y = mean, colour = print.spp)) + 
	geom_line(show.legend = FALSE) + facet_wrap(~Site) +
	theme_classic() + 
	ylab("Expected Number of Birds") +	xlab("Year")

res <- lapply(names(results), FUN = function(x){cbind(results[[x]]$summary.random$trend_effect, spp = x)})
res <- data.frame(do.call('rbind', res))
res <- res %>% 
  mutate(parameter = case_when(
    ID == "(Intercept)" ~ "Intercept Lake Poteriteri",
    ID == "SiteWaitutu" ~ "Waitutu Difference",
    ID == "SiteLake Poteriteri:year" ~ "Lake Poteriteri",
    ID == "SiteWaitutu:year" ~ "Waitutu"
  ))
res <- res %>% 
  left_join(time.results %>% 
  group_by(spp) %>% 
  slice(1) %>% 
  select(spp, print.spp))

res$colour <- "Includes Zero"
res$colour[(res$X0.975quant < 0)] <- "Strictly Negative"
res$colour[(res$X0.025quant > 0)] <- "Strictly Positive"

## RANK BY TREND  + Split Native vs. Not-Native
## Organizing the spp:
names <- read_excel("data/Waitutu_birds_for_analysis_Apr2022.xlsx")
res <- res %>% 
  mutate(spp = ifelse(spp == "Cuckoo, Long-tailed ", "Cuckoo, Long-tailed", spp))
res <- res %>% 
  mutate(spp = ifelse(spp == "Creeper, Brown ", "Creeper, Brown", spp))
res <- res %>% 
  left_join(names, by = c(spp = "5MBC_species_name"))

spp.ord <- res %>% 
  filter(parameter %in% c("Lake Poteriteri", "Waitutu")) %>%
  group_by(print.spp, Provenance) %>% 
  summarize(mean = mean(mean)) %>% 
  arrange(Provenance, mean)
res$print.spp.f <- factor(res$print.spp, levels = spp.ord$print.spp, ordered = TRUE)

vals <- spp.ord %>% filter(Provenance == "introduced")

res.trend <- res %>% 
  filter(parameter %in% c("Lake Poteriteri", "Waitutu"))
res.intercept <- res %>% 
  filter(parameter %in% c("Intercept Lake Poteriteri", "Waitutu Difference")) %>% 
  group_by(spp) %>% 
  mutate(sizeW = sum(mean), n = n()) %>%
  mutate(size = ifelse(parameter == "Waitutu Difference", exp(sizeW), exp(mean))) %>%
  ungroup() %>%
  mutate(parameter = ifelse(parameter == "Intercept Lake Poteriteri", "Lake Poteriteri", "Waitutu"))
res.trend <- res.trend %>% 
  left_join(res.intercept %>% 
    select(parameter, spp, size), by = c("parameter", "spp"))

res.trend <- res.trend %>% 
  mutate(Site = ifelse(parameter == "Waitutu", "Waitutu Coast", parameter))

g3 <- ggplot(data = res.trend, 
  aes(x = mean, y = print.spp.f, colour = colour)) +
	facet_wrap(~Site) +
  geom_point(size = 0.01) +  ## Need to put a point down but not sure why... Plotting on top of the rectanges is the point.
  geom_rect(aes(ymin = 0.5, ymax = 6.5,  xmin = -Inf, xmax = Inf), 
            inherit.aes = FALSE, size = 1, fill = 'grey92', alpha = .1) + 
  geom_rect(aes(ymin = 6.5, ymax = 22.5,  xmin = -Inf, xmax = Inf), 
            inherit.aes = FALSE, size = 1, fill = 'white', alpha = .1) +   ## Can change the plot colour if you want.
  geom_point(aes(size = size), shape = 16) +   
	geom_errorbar(aes(xmin = X0.025quant, xmax = X0.975quant), width = 0) +
	theme_classic() + theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin()) + 
  ylab("Species") + xlab("Annual Trend") +
	scale_colour_manual("Direction of Trend", 
		values = c("Includes Zero" = 'black', 'Strictly Negative' = 'red', 'Strictly Positive' = 'green3')) + 
  # scale_size("Expected Initial Count/Station") +
  scale_size_binned("Expected Initial Count/Station", range = c(0.25, 3.5), n.breaks = 5) +## PLAY AROUND WITH THIS.
	geom_vline(xintercept=0, linetype = 'dashed') +
  ylab("") +
  guides(size = guide_legend(order = 2, byrow = TRUE), colour = guide_legend(order = 1, byrow=TRUE)) + 
  # guides(size = guide_bins(axis.arrow = arrow(length = unit(1.5, 'mm'), ends = 'last'))) +
  theme(strip.background = element_blank()) + 
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold"))

g3

ggsave(g3, file="figures/SummaryPlot.png", width=200, height=150, units="mm", dpi = 150)

## Now make the management specific plots:
## BEFORE AFTER ANALYSIS:
res.mgmt <- lapply(names(results.mgmt), FUN = function(x){cbind(results.mgmt[[x]]$summary.random$trend_effect, spp = x)})
res.mgmt <- data.frame(do.call('rbind', res.mgmt))
res.mgmt <- res.mgmt %>% 
  mutate(parameter = case_when(
    ID == "(Intercept)" ~ "Intercept Lake Poteriteri",
    ID == "SiteWaitutu" ~ "Waitutu Difference",
    ID == "SiteLake Poteriteri:mgmt" ~ "Lake Poteriteri (Managed)",
    ID == "SiteWaitutu:mgmt" ~ "Waitutu Difference (Managed)"
  ))
res.mgmt <- res.mgmt %>% 
  left_join(time.results %>% 
  group_by(spp) %>% 
  slice(1) %>% 
  select(spp, print.spp)) %>%
  mutate(site = ifelse(grepl("Waitutu", parameter), "Waitutu", "Lake Poteriteri"))

res.mgmt <- res.mgmt %>% 
  mutate(spp = ifelse(spp == "Cuckoo, Long-tailed ", "Cuckoo, Long-tailed", spp)) %>%
  mutate(spp = ifelse(spp == "Creeper, Brown ", "Creeper, Brown", spp)) %>%
  left_join(names, by = c(spp = "5MBC_species_name"))

res.mgmt$colour <- "Includes Zero"
res.mgmt$colour[(res.mgmt$X0.975quant < 0)] <- "Strictly Negative"
res.mgmt$colour[(res.mgmt$X0.025quant > 0)] <- "Strictly Positive"
res.mgmt$print.spp.f <- factor(res.mgmt$print.spp, levels = spp.ord$print.spp, ordered = TRUE)

res.intercept <- res.mgmt %>% 
  filter(parameter %in% c("Intercept Lake Poteriteri", "Waitutu Difference")) %>% 
  group_by(spp, site) %>% 
  mutate(sizeW = sum(mean), n = n()) %>%
  mutate(size = ifelse(parameter == "Waitutu Difference", exp(sizeW), exp(mean))) %>%
  ungroup() %>%
  mutate(parameter = ifelse(parameter == "Intercept Lake Poteriteri", "Lake Poteriteri", "Waitutu"))
res.mgmt <- res.mgmt %>% 
  left_join(res.intercept %>% 
    select(spp, size, site), by = c("spp", "site"))

res.mgmt.trend <- res.mgmt %>% filter(grepl("Managed", parameter))

ggplot(data = res.mgmt.trend, 
  aes(x = mean, y = print.spp.f, colour = colour)) +
	facet_wrap(~parameter) +
  geom_point() +  ## Need to put a point down but not sure why... Plotting on top of the rectanges is the point.
  geom_rect(aes(ymin = 0.5, ymax = 6.5,  xmin = -Inf, xmax = Inf), 
            inherit.aes = FALSE, size = 1, fill = 'grey96', alpha = .1) + 
  geom_rect(aes(ymin = 6.5, ymax = 22.5,  xmin = -Inf, xmax = Inf), 
            inherit.aes = FALSE, size = 1, fill = 'azure1', alpha = .1) +   ## Can change the plot colour if you want.
	geom_errorbar(aes(xmin = X0.025quant, xmax = X0.975quant), width = 0) +
  geom_point(aes(size = size), shape = 16) +   	
  theme_classic() + theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin()) + 
  ylab("Species") + xlab("Annual Trend") +
	scale_colour_manual("Direction of Trend", 
		values = c("Includes Zero" = 'black', 'Strictly Negative' = 'red', 'Strictly Positive' = 'blue')) + 
  scale_size_binned("Expected Initial Count/Station", range = c(0.25, 3.5), n.breaks = 5) +## PLAY AROUND WITH THIS.
	geom_vline(xintercept=0, linetype = 'dashed') +
  guides(size = guide_legend(order = 2, byrow = TRUE), colour = guide_legend(order = 1, byrow=TRUE)) + 
  # guides(size = guide_bins(axis.arrow = arrow(length = unit(1.5, 'mm'), ends = 'last'))) 
ggsave("figures/SummaryPlot_BA.png", width=200, height=150, units="mm", dpi = 150)  ## Change this for plublication quality potentially.
