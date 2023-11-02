# This script aims to generate some of the supplementary figures for the paper
self_relative_path <- "src/make_supplementary_figures.R"

########################
cat("\nStart of:", self_relative_path, "\n")
########################

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
})

{
    # paths 
    gitdir <- here::here()
    source(file.path(gitdir, "R/paths.R"))
    source(file.path(gitdir.functions, "plotting_functions.R"))
}

save_figures <- FALSE
if( exists("OUTDIR")){
    save_figures <- TRUE
    outdir.figures <- file.path(OUTDIR, "figures")
    stopifnot(dir.exists(outdir.figures))
}

##################
#      Main      #
##################

naturemed_reqs()

# supplementary figure on smoothed participation rates...
dcens <- fread(path.census.eligible)
dpart <- readRDS(path.participation.rates)
dcens[, FC := fifelse(COMM_INDEX == 1, yes="fishing", no="inland")]

# Horizonal looks better, but not for displaying the census eligible sizes in fishing communities
facet_formula <- formula(ROUND_LAB ~ FC_LAB)
horizontal <- FALSE
if(horizontal){
    facet_formula <- formula( FC_LAB ~ ROUND_LAB )
}

round_labs_copy <- round_labs

round_labs <- gsub(" to ", "-", round_labs) |>
    gsub(pattern="Jan ", replacement="01/", x=_) |>
    gsub(pattern="Feb ", replacement="02/", x=_) |>
    gsub(pattern="May ", replacement="05/", x=_) |>
    gsub(pattern="Jun ", replacement="06/", x=_) |>
    gsub(pattern="Jul ", replacement="07/", x=_) |>
    gsub(pattern="Sep ", replacement="09/", x=_) |>
    gsub(pattern="Oct ", replacement="10/", x=_) |>
    gsub(pattern="20", replacement="", x=_)

p_smoothed_cens <- plot.pyramid.bysexround(
    dcens[ROUND >= 16], 
    NUM = 'ELIGIBLE',
    DEN='ELIGIBLE',
    .ylab = "",
    .facet = facet_formula,
    percent_lab = FALSE
) +
geom_line(aes(y=ELIGIBLE_SMOOTH * (-1 + 2*as.integer(SEX == "F"))), color="purple" ) +
labs( x="Census Eligible individuals", y="Age") +
theme(legend.key.size = unit(0.5, "cm")) 

dplot <- subset(dpart, ROUND >= 16) |>
    prettify_labels()

p_smoothed_part <- ggplot(dplot, aes(x = AGEYRS, y=PARTRATE_RAW, color=SEX_LAB)) +
geom_point() +
geom_line(aes(y=PARTRATE_SMOOTH.25)) +
facet_grid(facet_formula, scales = "free", labeller = labeller(ROUND_LAB = round_labs, FC_LAB = community_dictionary$long)) +
scale_color_manual(labels=plabels$sex ,values = palettes$sex) +
scale_y_continuous(expand = c(0,0.01), limits=c(0,1), labels = scales::percent) +
theme_default() +
my_labs(y = "Proportion of census eligble individuals who participated") +
NULL

# if(!horizontal){
#     p <- (p_smoothed_cens + nm_reqs + theme(strip.text.y = element_blank())| p_smoothed_part + nm_reqs)
# }else{
# add letters a), b)
p <- (p_smoothed_cens + nm_reqs + labs(y=NULL) )/ (p_smoothed_part + nm_reqs + theme(strip.text.x = element_blank())) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = "a")
# }

if( save_figures ){
    filename <- "suppfig_smoothedcens_smoothedpart.pdf"
    cmd <- ggsave_nature(p, filename=filename, LALA=outdir.figures, w=18, h=24.5)
    system(cmd)
    system(zathura2gthumb(cmd))
}
