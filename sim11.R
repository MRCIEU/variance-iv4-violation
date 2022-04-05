library("dplyr")
library("ggplot2")
library("grid")
library("gtable")
library("TwoSampleMR")
library("viridis")
source("funs.R")
set.seed(123)

n_sim <- 5
r2_u <- 0.05 # U-X main effect
x_b <- log(1.2) # X-Y main effect
u_b <- log(1.2) # U-Y main effect
xu_b <- x_b * 0.5 # X*U-Y effect half the size of the main effect

results <- data.frame()
for (n_obs in c(10000, 50000, 100000)){
    for (r2_z in seq(0.01, 0.05, 0.01)){ # variance explained by main effect of Z-X
        for (phi in seq(0, 2.5, 0.5)){ # size of Z-X interaction effect relative to main effect
            r2_zu <- r2_z * phi
            u_b <- sqrt(r2_u)
            z_b <- sqrt(r2_z)
            zu_b <- sqrt(r2_zu)
            for (i in 1:n_sim){
                # simulate variables
                z <- get_simulated_genotypes(0.25, n_obs); z <- scale(z)
                u <- rnorm(n_obs)
                x <- z*z_b + u*u_b + z*u*zu_b + rnorm(n_obs, sd=sqrt(1-(r2_z+r2_u+r2_zu))); x <- scale(x)

                # simulate binary outcome
                p <- x*x_b + u*u_b + x*u*xu_b
                p <- 1/(1+exp(-p))
                y <- rbinom(n_obs, 1, p)

                # MR effect
                iv_exp <- lm(x ~ z)
                iv_out <- glm(y ~ z, family="binomial")
                b_exp <- iv_exp %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(estimate)
                se_exp <- iv_exp %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(std.error)
                b_out <- iv_out %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(estimate)
                se_out <- iv_out %>% tidy %>% dplyr::filter(term=="z") %>% dplyr::pull(std.error)
                wald <- mr_wald_ratio(b_exp, b_out, se_exp, se_out, NULL)

                # IV-exp variance effect
                result <- model(data.frame(z, x), "z", "x")

                # store result
                result$n_obs <- n_obs
                result$r2_z <- r2_z
                result$r2_zu <- r2_zu
                result$xu_b <- xu_b
                result$b_mr <- wald$b
                result$b <- x_b
                result$phi <- phi
                result$bias <- result$b_mr / result$b

                results <- rbind(results, result)
            }
        }
    }
}

# save data
write.table(file="data/sim11.txt", results)

# estimate power
pwr <- results %>% 
    dplyr::group_by(n_obs, r2_z, phi) %>% 
    dplyr::summarize(cbind(binom.test(sum(phi_p < 0.05), n()) %>% tidy, bias=median(bias))) %>%
    as.data.frame

# plot
pwr$n_obs <- as.factor(pwr$n_obs)
p <- ggplot(pwr, aes(x=n_obs, y=estimate, ymin=conf.low, ymax=conf.high, color=bias)) +
    geom_line() + 
    geom_point() + 
    geom_errorbar() +
    theme_classic() + 
    labs(x = "Sample size", y = expression(Power ~ (alpha==0.05)), color = "Fold-change") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0,.25,.5,.75,1)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey") +
    scale_color_viridis(direction = 1) +
    scale_x_discrete(labels = c('10,000','50,000','100,000')) +
    facet_grid(phi~r2_z) +
    ggtitle("Instrument-exposure main effect variance explained") +
    theme(
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        panel.spacing = unit(1, "lines"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5, size=11)
    )

# text, size, colour for added text
text = "Size of instrument-exposure interaction effect relative to main effect"
size = 11
col = "black"

# Convert the plot to a grob
gt <- ggplot2::ggplotGrob(p)

# Get the positions of the right strips in the layout: t = top, l = left, ...
strip <-c(subset(gt$layout, grepl("strip-r", gt$layout$name), select = t:r))

# Text grob
text.grob = grid::textGrob(text, rot = -90, gp = gpar(fontsize = size, col = col))

# New column to the right of current strip
# Adjusts its width to text size
width = unit(2, "grobwidth", text.grob) + unit(1, "lines")
gt <- gtable::gtable_add_cols(gt, width, max(strip$r))  

# Add text grob to new column
gt <- gtable::gtable_add_grob(gt, text.grob, 
        t = min(strip$t), l = max(strip$r) + 1, b = max(strip$b))

pdf("sim11.pdf")
# Draw it
grid.newpage()
grid.draw(gt)
dev.off()