#############
# INLA BRU spatial analysis at Waitutu
#############

# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library(inlabru)
library(fmesher)
library(INLA)
library(tidyverse)
library(ggplot2)
library(readxl)
library(sf)
library(concaveman)
library(cowplot)

## Use here to get path to proj.r
library(here)
setwd(here())

## Seed data from James
seeds <- read.csv(here("data/Seed_mast_summary_Waitutu_2011_2022.csv"))
seeds <- seeds %>% 
	mutate(seed = factor(Mast, levels = c("None", "Partial", "Full"))) %>% 
	mutate(seed = as.numeric(seed))
seed <- seeds %>% 
  group_by(site, year) %>% 
  summarize(seed = max(seeds_m2), .groups = "drop")

## Load Processed Bird Data
load(here("data/processed_data.Rda"))
bird.count <- dat.all$dat

## Additional variables required:
## mgmt occurs after 1080 operation that was 2010.
bird.count@data <- bird.count@data %>% 
	mutate(site_station = paste0(Site,Station)) %>%
	mutate(mgmt = ifelse(Year > 2010, 1, 0)) %>%
  mutate(mgmtyear = mgmt*(Year - min(Year[mgmt==1]) + 1)) ## Piecewise linear slope.

## Choose which species to analyze, or can do them all in the loop.
## --------------------------------
## Check bird list: 
## print(unique(bird.count$spp))
## Copy and Paste Candidates from Here:
## 1) "Robin, *NZ", 2) Tomtit*, 3) Bellbird* (mainland) 
## 4) Fantail, *Grey / NZ / Black 5) Kaka* 6) Silvereye / Waxeye / White-eye
## 7) Tui* 8) Warbler, *Grey.
#####################################################
spp.consider <- c("Blackbird", "Bellbird* (mainland)", "Robin, *NZ", "Tomtit*", "Fantail, *Grey / NZ / Black",
	"Kaka*", "Silvereye / Waxeye / White-eye", "Tui*", "Warbler, *Grey", "Chaffinch", "Rifleman*", 
	"Parakeet* spp. / Kakariki spp", "Cuckoo, Long-tailed ", "Pigeon, *NZ", "Thrush, Song", "Dunnock / Hedge Sparrow",
		"Cuckoo, Shining", "Greenfinch", "Kea", "Falcon, *NZ", "Redpoll", "Creeper, Brown ")  ## Hanging white in 'Brown ' from database.

## Special characters for Maori names.
# a = \U0101
# i = \U012B
# u = \U016B

spp.print <- c("Eurasian blackbird / manu pango", "Bellbird / korimako", "South Island robin / toutouwai", "Tomtit / miromiro", 
	"New Zealand fantail / p\U012Bwakawaka","K\U0101k\U0101", "Silvereye / tauhou", 
	"T\U016B\U012B", "Grey warbler / riroriro", "Chaffinch / pahrini", "Rifleman / t\U012Btitipounamu", 
	"Parakeet spp. / k\U0101k\U0101riki", "Long-tailed cuckoo / koekoe\U0101", "New Zealand pigeon / kerer\U016B", 
	"Song thrush / manu-kai-rakau", "Dunnock",
	"Shining cuckoo / p\U012Bp\U012Bwharauroa", "European greenfinch", "Kea", "New Zealand falcon / k\U0101rearea", "Common redpoll", "Brown creeper / pipipi")

## Build the spatial component of the model:
## inner and outer bounds and then the mesh:
df.sp <- bird.count@data %>%
  distinct(coords.x1, coords.x2,  Site, Station) %>%
  mutate(x1 = coords.x1 / 1000, x2 = coords.x2 / 1000)  ## convert to km.
df.sp <- st_as_sf(df.sp, coords = c("x1", "x2"))

bounds.inner <- df.sp %>% 
  st_buffer(dist = 250/1000, endCapStyle = "SQUARE") %>%
  group_by(Site) %>%
  summarize(geometry = st_union(geometry))
  
bounds.outer <- bounds.inner %>% 
  st_buffer(2)

## Use fmesher to build a mesh:
mesh <- fm_mesh_2d_inla(
  loc = st_as_sfc(df.sp),
  boundary = list(bounds.inner, bounds.outer),
  max.edge=c(155/1000, 2000/1000),
  min.angle=c(30, 21),
  max.n=c(48000, 16000), ## Safeguard against large meshes.
  max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
  cutoff=155/1000, ## Filter away adjacent points.
  offset=c(270/1000, 930/1000)
)

## Priors for SPDE
## Priors are defined as upper values on sigma and range.
matern  <- inla.spde2.pcmatern(mesh, 
              prior.sigma =c(3,0.01), 
              prior.range = c(500/1000, 0.01))

## Years to plot:
yearStart <- 3
yearEnd <- 16

## Prediction grid for plotting after fitting:
grd.sp <- df.sp %>% 
  st_make_grid(what = "centers", cellsize = 0.1) %>%
  st_as_sf() %>%
  st_intersection(bounds.inner)
grd.sp <- grd.sp %>% mutate(year = yearStart) %>%
  bind_rows(grd.sp %>% mutate(year = yearEnd))

## Cache the results in loop:
results <- list()
results.mgmt <- list()
time.results <- data.frame()

## If you want to produce the residual plots.
run.checks <- FALSE

## If you want to save the spatial and temporal plots
make.plots <- FALSE

for(i in 1:length(spp.consider))
{
  ## Subset and process bird specific data frame:
	spp.id <-  spp.consider[i]
	print.id <- spp.print[i]
  print(paste("Analyzing", print.id, sep = ": "))
	df <- bird.count@data %>% 
    filter(spp == spp.id)
	df <- df %>% 
    mutate(year = as.integer(Year) - min(as.integer(Year)))
	df$seed <- as.factor(df$seed)
	df <- df %>% 
		mutate(fyear = as.factor(Year)) %>%
		mutate(seed2 = factor(seed, levels = c(1,2,3), 
      labels = c("none", "partial", "full")))
  
  df <- df %>%
    mutate(x1 = coords.x1 / 1000, x2 = coords.x2 / 1000) %>%
    st_as_sf(coords = c("x1", "x2"))
	
	df <- df %>% 
    mutate(reSiteYear = factor(paste(Site, year)))

  ## Model Definition: remove intercept as it is in trend_effect
  mod <- count ~ -1 + grf(main = st_coordinates, model = matern) + 
                 trend_effect(main = ~ Site + year:Site, model = "fixed") +
                 re_site_year(main = reSiteYear, model = "iid") ## Site:Year random effect

  ## Model Fitting:
	bru.spp <- bru(
	  mod,
	  family = "poisson",
	  data = df
	)
  # bru.spp$summary.random$trend_effect ## Summary of fixed effects.

  ## Alternative Definition: Before-After Analysis (Average abundance <= 2010 compared with average abundance after 2010)
  mod2 <- count ~ -1 + grf(main = st_coordinates, model = matern) + 
                 trend_effect(main = ~ Site + mgmt:Site, model = "fixed") +
                 re_site_year(main = reSiteYear, model = "iid") ## Site:Year random effect

  ## Model Fitting:
	bru.spp.mgmt <- bru(
	  mod2,
	  family = "poisson",
	  data = df
	)
  # bru.spp.mgmt$summary.random$trend_effect 
  
	results[[spp.id]] <- bru.spp
	results.mgmt[[spp.id]] <- bru.spp.mgmt

  # bru.spp <- results[[spp.id]]  ## For James to skip bru line ****

  ## Prediction:
  pred <- predict(
    bru.spp, newdata = grd.sp,
    ~ exp(grf + trend_effect)
  )
  
  pred$Year <- pred$year + min(bird.count@data$Year)

	# ggplot() +
	  # geom_tile(pred, mapping = aes(geometry = x, fill = mean, colour = mean), size = 4,  stat = "sf_coordinates") +
	  # facet_wrap(~ Site + Year, scale = "free") + 
	  # scale_fill_viridis_c(name = "Expected # \nof birds",  limits=c(0,ceiling(max(pred$mean)))) +   
	  # scale_colour_viridis_c(name = "Expected # \nof birds",  limits=c(0,ceiling(max(pred$mean)))) + 
	  # theme_classic() +   
    # xlab("") + ylab("") + 
	  # theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
			# axis.ticks.x=element_blank(),axis.ticks.y=element_blank()) +
	  # ggtitle(print.id)
	# ggsave(here("figures", paste0(gsub("\\*|/|\\,", "", spp.id),"SpatPlot.png")), width = 9.45, height = 8.29)

  g1 <- ggplot() +
    geom_tile(pred, mapping = aes(geometry = x, fill = mean, colour = mean), size = 4,  stat = "sf_coordinates") +
    facet_wrap(~ Site + Year, scale = "free", nrow=2) +
  #  facet_grid(Site~Year, scale="free") +
    scale_fill_viridis_c(name = "Expected \nno. of birds",  limits=c(0,ceiling(max(pred$mean)))) +   
    scale_colour_viridis_c(name = "Expected \nno. of birds",  limits=c(0,ceiling(max(pred$mean)))) + 
    theme_bw() +
    theme(axis.line = element_blank(),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()) +
  #    theme_classic() +   
    xlab("") + ylab("") + 
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
      axis.ticks.x=element_blank(),axis.ticks.y=element_blank()) +
    theme(strip.background = element_blank()) + 
    theme(strip.text.x = element_text(size = 0, color = "white", face = "bold.italic")) +
    ggtitle(print.id) + theme(plot.title = element_text(hjust = 0.5))

  g1x <- ggdraw(g1) + 
    draw_label("2009", x = 0.23, y = 0.895) +
    draw_label("2022", x = 0.63, y = 0.895) +
    draw_label("Lake Poteriteri", colour = "black", size = 14, angle = 90, x = 0.03, y = 0.68) +
    draw_label("Waitutu Coast", colour = "black", size = 14, angle = 90, x = 0.03, y = 0.25)


  if(make.plots) ggsave(g1x, file=paste0("./figures/",gsub("\\*|/|\\,", "", spp.id),"_SpatialPlot.png"), width=170, height=110, units="mm", dpi = 150)

	# Now predict average per year:
	avgs <- df %>% 
    group_by(Site, Year, year, mgmt) %>% 
		summarize(pa = mean(count > 0), avg = mean(count, na.rm = TRUE), avg_p = mean(count[count > 0]), .groups = "drop")
	pred.df <- df[!duplicated(paste(df$Site, df$Year)),]
	pred.time <- predict(bru.spp, pred.df, ~ exp(trend_effect))

	pred.time <- pred.time %>% 
    as.data.frame() %>%
    select(-geometry) %>%
    left_join(as.data.frame(avgs) %>% select(-geometry), by = c("Year", "Site", "year", "mgmt"))
  
	# ggplot(data = pred.time, aes(x = Year, y = avg)) + 
    # geom_point() + facet_wrap(~Site) + 
		# geom_line(aes(y = mean)) + 
    # geom_ribbon(aes(ymin=q0.025,ymax=q0.975),alpha=0.3) +
		# annotate("label", x = 2010, y = 0, label = "First 1080 Application", hjust = "left") +
    # geom_vline(xintercept=2010, col = "red", linetype = "dashed") +
		# theme_bw() + 
		# ylab("Average Number of Birds") +
		# ggtitle(print.id) 

  g2 <- ggplot(data = pred.time, aes(x = Year, y = avg)) + 
    geom_point() + facet_wrap(~Site, scales='free') + 
    geom_line(aes(y = mean)) + 
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()) +
    geom_ribbon(aes(ymin=q0.025,ymax=q0.975),alpha=0.3) +
    #annotate("label", x = 2010, y = 0, label = "First 1080 Application", hjust = "left") +
    geom_vline(xintercept=2010, col = "grey50", linetype = "dashed") +
    theme(strip.background = element_blank()) + theme(strip.text = element_text(size=14)) +
    theme(panel.background=element_blank()) + theme(panel.border=element_blank()) +
    theme(axis.line.y = element_line(colour="black", linewidth=0.5, linetype="solid")) + 
    theme(axis.line.x = element_line(colour="black", linewidth=0.5, linetype="solid")) +
    scale_x_continuous(limits=c(2006,2022)) + scale_y_continuous(expand=c(0,0),limits=c(0,max(pred.time$q0.975))) +    #max(pred.time$q0.975
    ylab("Mean count") + expand_limits(y=0) + 
    theme(axis.title = element_text(size=14), axis.text = element_text(size=12)) +
    ggtitle(print.id) + theme(plot.title = element_text(hjust = 0.5, face="bold"))

  g2

  if(make.plots) ggsave(g2, file=here("figures", paste0(gsub("\\*|/|\\,", "", spp.id), "TemporalPlot.png")), width=150, height=90, units="mm", dpi = 150)

	time.results <- rbind(time.results, pred.time)

	# Model checking code
  if(run.checks){
    # Check total count prediction for each site for each year
    # Total observed abundance in each year
    n.obs = df %>% 
      group_by(Site, Year) %>% 
      summarise(n = sum(count))

    # sample poisson rate for each obs location each year
    # Incrase n.samples for final figures, probably want more than
    # 100 but this is fast for developing the plots.
    n.pred = generate(bru.spp,
              data = df,
              formula = ~ exp(trend_effect + re_site_year),
              n.samples = 500)

    # for each simulated poisson rate sample a single count
    # from rpois
    blah = cbind(df %>% data.frame() %>% select(Site, Year, -geometry), n.pred)
    blah2 = blah %>% 
      pivot_longer(cols = -c("Site", "Year"),
             names_to = "sample.id",
             values_to = "rate") %>%
      arrange(Site, Year, sample.id)

    blah2$count = rpois(lambda = blah2$rate,
              n = nrow(blah2))
    blah2$sample.id = as.factor(blah2$sample.id)

    blah3 = blah2 %>% 
      group_by(Site, Year, sample.id) %>% 
      summarise(n = sum(count),
          .groups = "keep") %>%
      arrange(Site, Year)

    # boxplot of posterior counts
    # observed count marked in red
    ggplot(blah3) +
      geom_boxplot(aes(y = n, x = Site)) + 
      geom_point(data = n.obs,
           aes(y = n, x = Site),
           colour = "red") +
      ylab("count") +
      ggtitle(spp.id) +
      facet_wrap( ~ Year)

    if (!dir.exists(here("figures", "model_checking"))){
      dir.create(here("figures", "model_checking"))
    }

    ggsave(here("figures", "model_checking", paste0(gsub("\\*|/|\\,", "", spp.id), "_site_year_counts.png")),
         width = 9.45, 
         height = 8.29)

    # spatial prediction diagnostic plot
    # look at mean squared difference
    # of predicted count vs true count at each sampling location

    blah = cbind(df %>% data.frame() %>% select(Site, Station, Year, count, coords.x1, coords.x2, -geometry), n.pred) %>% 
      mutate(x1 = coords.x1 / 1000, x2 = coords.x2 / 1000) %>%
      dplyr::select(-coords.x1, -coords.x2)

    blah2 = blah %>% 
      pivot_longer(cols = -c("Site", "Station", "Year", "count", "x1", "x2"),
             names_to = "sample.id",
             values_to = "rate") %>%
      mutate(true_count = count) %>%
      arrange(Site, Year, sample.id)

    blah2$post_sampled_count = rpois(lambda = blah2$rate,
                     n = nrow(blah2))
    blah2$sample.id = as.factor(blah2$sample.id)
    blah3 = blah2 %>% 
      mutate(sq_diff = (post_sampled_count - true_count)^2) %>%
      group_by(Year, Station) %>%
      summarise(mean_sq_diff = mean(sq_diff),
          .groups = "keep")
    head(blah3)

    blah4 = blah3 %>%
      left_join(unique(blah[, c("Year", "Station", "x1", "x2", "Site")]), by = c("Year", "Station")) %>%
      arrange(Year, Station)

    # Too many for a single plot so
    # do each site separately
    blah4 %>%
      filter(Site == "Waitutu") %>%
      as.data.frame() %>%
      ggplot() +
      geom_point(aes(x = x1, y = x2, colour = mean_sq_diff)) +
      facet_wrap(~ Year) +
      scale_colour_viridis_c() +
      coord_equal()
     
    # ouch, some locations have large values, dominating colour scale
    # try sqrt of mean square diff

    blah4 %>%
      filter(Site == "Waitutu") %>%
      as.data.frame() %>%
      ggplot() +
      geom_point(aes(x = x1, y = x2, colour = sqrt(mean_sq_diff))) +
      facet_wrap(~ Year) +
      scale_colour_viridis_c("mean difference") +
      coord_equal() +
      xlab("Easting") +
      ylab("Northing") + 
      ggtitle("Waitutu")

    # Think of this plot as analogous to spatial residuals
    ggsave(here("figures", "model_checking", paste0(gsub("\\*|/|\\,", "", spp.id), "_Waitutu_spatial_residuals.png")),
         width = 9.45, 
         height = 8.29)

    blah4 %>%
      filter(Site == "Lake Poteriteri") %>%
      as.data.frame() %>%
      ggplot() +
      geom_point(aes(x = x1, y = x2, colour = sqrt(mean_sq_diff))) +
      facet_wrap(~ Year) +
      scale_colour_viridis_c("mean difference") +
      coord_equal() +
      xlab("Easting") +
      ylab("Northing") + 
      ggtitle("Waitutu")

    # Think of this plot as analogous to spatial residuals
    ggsave(here("figures", "model_checking", paste0(gsub("\\*|/|\\,", "", spp.id), "_Poteriteri_spatial_residuals.png")),
         width = 9.45, 
         height = 8.29)
  } 
}

# time.results <- time.results %>% left_join(data.frame(spp = spp.consider, print.spp = spp.print), by = 'spp')

# save(results, file = "output/results.Rda")
# save(results.mgmt, file = "output/results_mgmt.Rda")
# save(time.results, file = "output/temporalResults.Rda")