#####################
# Set default theme #
#####################

geom.text.size = 5/14*7
update_geom_defaults("text", list(size = geom.text.size))

# update_geom_defaults("text", list(size = 3))
GeomLabel$default_aes$family <- 'sans'

# set default theme to bw and with better looking facets.
theme_default <- function(...) {
    theme_bw() +
    theme( 
        legend.position = "bottom",
        strip.background = element_blank(),
        ...
    )
}

theme_debug <- function(...) {
    theme(
        legend.background = element_rect(fill =  "deepskyblue2"),
        plot.background = element_rect(fill = "lemonchiffon", colour = "black", size = 2),
        panel.background = element_rect(fill = "purple", colour = "red", size = 2),
        legend.key = element_rect(fill = "orange"),
        panel.border = element_rect(colour = "black", fill = NA),
        ...
    )
}
t_nomargin <- theme(plot.margin = unit(c(0,0,0,0), "cm"))


theme_set( theme_default() )


#########################
# Requirements + saving #
#########################

require(ggplot2)

p_reqs <- theme(
    plot.title=element_text(size=12, family='sans'),
    text=element_text(size=10, family = 'sans'),
    axis.title=element_text(size=10, family='sans'), 
    axis.text=element_text(size=7, family = 'sans')
)

ggsave2 <- function(p, file, w, h, LALA=vl.out.dir, u='in')
{
    if(! file %like% 'pdf$')
        stop("`ggsave2` file should end in `.pdf`")
    # I used LALA because I think there may be some naming issues
    filename <- file
    filename2 <- gsub('pdf$','png',filename)
    cat('Saving', filename, '...\n')
    ggsave(p, filename=file.path(LALA,filename2), width=w, height=h, units=u)
    ggsave(p, filename=file.path(LALA, filename), width=w, height=h, units=u)	
    zathura_cmd <- paste0('zathura ', file.path(LALA, filename), ' &')
    return(zathura_cmd)
}

naturemed_reqs <- function() 
{
    # call this before doing your plots
    nm_reqs <<- theme(
        axis.text = element_text(size=5, family='sans'),
        text=element_text(size=7,family='sans'),
        legend.text=element_text(size=7, family='sans'),
        strip.background = element_rect(colour="white", fill="white"),
        plot.tag = element_text(size=8, face='bold', color = "black"),
        plot.tag.position = c(0, 1),
        strip.placement = "outside"
    )

    slides_reqs <<- theme( 
        axis.text = element_text(size=6, family='sans'),
        text=element_text(size=9,family='sans'),
        legend.text=element_text(size=9, family='sans'),
        strip.background = element_rect(colour="white", fill="white"),
        plot.tag = element_text(size=10, face='bold', color = "black"),
        plot.tag.position = c(0, 1),
        strip.placement = "outside"
    )
}

ggarrange_nature <- function(
  ...,
  plotlist = NULL,
  ncol = NULL,
  nrow = NULL,
  labels = NULL,
  label.x = 0,
  label.y = 1,
  hjust = -0.5,
  vjust = 1.5,
  align = c("none", "h", "v", "hv"),
  widths = 1,
  heights = 1,
  legend = NULL,
  common.legend = FALSE,
  legend.grob = NULL, 
  add_reqs=TRUE
){
    if(add_reqs)
        reqs <- theme(axis.text = element_text(size=5, family='sans'), text=element_text(size=7,family='sans'), legend.text=element_text(size=7, family='sans'))

    plots <- c(list(...), plotlist)

    if(add_reqs)
        plots <- lapply(plots, function(p){p + reqs})

    out <- ggarrange(plotlist = plots,
                     ncol = ncol,
                     nrow = nrow,
                     labels = labels,
                     label.x = label.x,
                     label.y = label.y,
                     hjust = hjust,
                     vjust = vjust,
                     font.label = list(size = 8, color = "black", face = "bold", family = 'sans'),
                     align = align,
                     widths = widths,
                     heights = heights,
                     legend = legend,
                     common.legend = common.legend,
                     legend.grob = legend.grob)
    return(out)
}

ggsave_nature <- function(filename, p, LALA=vl.out.dir, w=18,h=24, add_reqs=TRUE)
{
    # check size
    tmp <- sort(c(w,h))
    if(tmp[1] > 18 | tmp[2] > 24)
        warning('Plot is bigger than allowed for EDFs. Maximum size is 18cm x 24cm\n')
    if( tmp[1] < 10)
        warning('w and h represent cm units, not inches. Are you sure you want to save such a small plot?\n')

    if(add_reqs)
        p <- p + nm_reqs

    # I used LALA because I think there may be some naming issues
    filename2 <- gsub('pdf$','png',filename)
    cat('Saving', filename, '...\n')
    ggsave(p, filename=file.path(LALA,filename2), width=w, height=h, units='cm')
    ggsave(p, filename=file.path(LALA, filename), width=w, height=h, units='cm')	
    # save
    ggsave(filename=filename, plot=p, width=w, height=h, units='cm', dpi=310)

    # return command to open up the image with zathura
    zathura_cmd <- paste0('zathura ', file.path(LALA, filename), ' &')
    return(zathura_cmd)
}

####################
# ggplot helper shortcuts
####################

scale_y_expand_lower0 <- scale_y_continuous(expand = expansion(mult = c(0, .1)))
scale_y_percentage <- scale_y_continuous(labels=scales::label_percent(), expand=expansion(mult=0)) 
scale_y_percentage2 <- scale_y_continuous(labels=scales::label_percent(), expand=expansion(mult=0.05)) 
scale_y_percentagef <- function(x)
    scale_y_continuous(
        labels=scales::label_percent(),
        expand=expansion(mult=c(0,x)),
        breaks = c(0,0.25,0.5,0.75,1)
)

####################
# specify palettes #
####################

palettes <- list(

    sex=c(
        F="#F4B5BD",
        M="#9ac0cd",
        Female="#F4B5BD",
        Male="#9ac0cd"
        ),

    comm = c(
        fishing="#9ac0cd",
        inland="#9C964A",
        Fishing="#9ac0cd",
        Inland="#9C964A",
        `Fishing communities with high HIV prevalence` ="#9ac0cd",
        `Inland communities with typical HIV prevalence` ="#9C964A"
    ),

    arvmed = c(
        '#E03531',
        '#51B364'),

    cuarvmed = c(
        '#E03531',
        '#51B364'),

    supp_status =  c(
        "#FF585D",
        "#A9A9A9",
        "#3EB595"),

    comm2 = c(
        "#9C964A",
        "#9ac0cd",
        "#FF9933"),

    bool = c(
        "FALSE"='#ff9999',
        "TRUE"='#99ff99',
        "positive"='#ff9999',
        "negative"='#99ff99',
        `NA`='grey'),

    sexagegroup = c(
        `Male 45-49`= "#0e1c31",
        `Male 40-44`= "#1d3862",
        `Male 35-39`= "#2c5493",
        `Male 30-34`= "#3a70c4",
        `Male 25-29`= "#6b93d2",
        `Male 20-24`= "#9cb7e1",
        `Male 15-19`= "#cddbf0",
        `Female 45-49`="#2e0936",
        `Female 40-44`="#5d126c",
        `Female 35-39`="#8c1ba3",
        `Female 30-34`="#bb25d9",
        `Female 25-29`="#cc5be3",
        `Female 20-24`="#dd92ec",
        `Female 15-19`= "#eec8f5"
    ), 

    hivstatus =  c(
        `HIV positive` = "#FF585D",
        `NA` = "#A9A9A9",
        `HIV negative` = "#3EB595"
    ),

    rakailogo = c(
        "#0777bf",
        "#626263",
        "#D81921",
        "#000000",
        "#FFFFFF"
    ),

    uganda = c(
        "#000000",
        "#fcdc04",
        "#d90000",
        "#ffffff",
        "#9ca69c"
    ),


    minimal = c(
        "#0F4392",
        "#FF4949"
    ),

    minimal2 = c(
        "#374E55",
        "#DF8F44"
    ),

    # minimal3 = c("#FAE48BFF","#FB6467FF"),
    minimal3 = c("purple","darkgreen"),


    ftp = c(
        `All participants`="darkred",
        `First-time participants`="orange"
    ),

    NULL
)


###################
# prettify labels #
###################

# prettify_... functions add a variable to the data table. 
# labs_from_dictionaries returns a function that transforms variables to their labels through the use of dictionaries.

remove_pretty <- function(DT)
{
    nms <- names(DT) %which.like% '_LAB$'
    DT[, (nms) := NULL ]
    return(DT)
}

prettify_round <- function(DT)
{
    nms <- names(DT)
    if("ROUND" %in% nms & ! "ROUND_LAB" %in% nms)
    DT[, ROUND_LAB :=  round_dictionary[as.character(ROUND)]]
    return(DT)
}

prettify_sex <- function(DT)
{
    nms <- names(DT)
    if("SEX_LAB" %in% nms)
    return(DT)
    
    if("SEX_LABEL" %in% nms){
        DT[, SEX_LAB :=  sex_dictionary[as.character(SEX_LABEL)]]

    }else if("SEX" %in% nms){
        DT[, SEX_LAB :=  sex_dictionary[as.character(SEX)]]
    }

    return(DT)
}

prettify_participation_status <- function(DT)
{
    nms <- names(DT)
    if('PARTICIPATION_LAB' %in% nms | (! 'PARTICIPATION_STATUS' %in% nms ) )
        return(DT)

    dict <- postproc_dictionaries$PARTICIPATION_STATUS
    DT[, PARTICIPATION_LAB := dict[PARTICIPATION_STATUS] ]
    return(DT)
}

prettify_hivstatus <- function(DT)
{
    nms <- names(DT)
    if('HIV LAB' %in% nms | ! 'HIV' %in% nms)
        return(DT)

    DT[, HIV_LAB := data.table::fcase(
        HIV==1, 'HIV positive',
        HIV==0, 'HIV negative',
        default = NA_character_
    ) ]

    DT
}


prettify_fc <- function(DT)
{
    nms <- names(DT)
    if('FC_LAB' %in% nms | ! "FC" %in% nms)
        return(DT)

    DT[, FC_LAB := data.table::fcase(
        FC == 'inland', "Inland",
        FC == 'fishing', "Fishing",
        FC == 'trading', "Trading",
        FC == 'agrarian', "Agrarian",
        default = NA_character_
    )]

    DT
}

prettify_loc <- function(DT){

    nms <- names(DT)
    if('LOC_LAB' %in% nms | ! ( "LOC_LABEL" %in% nms | "LOC" %in% nms | "TYPE" %in% nms ))
        return(DT)

    if( "LOC_LABEL" %in% nms){
        DT[, LOC_LAB := loc_dictionary[as.character(LOC_LABEL)]]
    }else if('LOC' %in% nms){
        DT[, LOC_LAB := loc_dictionary[as.character(LOC)]]
    }else if('TYPE' %in% nms){
        DT[, LOC_LAB := loc_dictionary[as.character(TYPE)]]
    } 

    DT
}

prettify_labels <- function(DT)
{
    # preforms all of the above
    nms <- names(DT)

    DT |>  
        prettify_sex() |> 
        prettify_round() |> 
        prettify_hivstatus() |> 
        prettify_fc() |>
        prettify_loc()
}

labs_from_dictionaries <- function(dict)
{
    f <- function(keys)
    {
        NA2textNA <- function(x)
        {
            x[is.na(x)] <- "NA"
            return(x)
        }
        keys <- as.character(keys) |> NA2textNA()
        dict[keys]
    }

    return(f)
}

my_labs <- function (..., title = waiver(), subtitle = waiver(), caption = waiver(), tag = waiver(), alt = waiver(), alt_insight = waiver()) 
{
    # requires my_labs_dictionary

    p <- ggplot_build(last_plot())
    
    # copied from labs
    args <- rlang:::dots_list(..., title = title, subtitle = subtitle, 
        caption = caption, tag = tag, alt = alt, alt_insight = alt_insight, 
        .ignore_empty = "all")
    is_waive <- vapply(args, ggplot2:::is.waive, logical(1))
    args <- args[!is_waive]
    args <- args[!duplicated(names(args))]
    args <- ggplot2:::rename_aes(args)

    # only consider aesthetics with only one variable used.
    mapps <- lapply(p$plot$mapping, all.vars)
    idx <- sapply(mapps, length) == 1
    mapps_one_var <- mapps[idx]
    # exclude args defined in ...
    mapps_final <- mapps[! names(mapps) %in% names(args)]
    # transform labels according to my_labs_dictionary
    mapps_final <- lapply(mapps_final, function(x){
        unname(my_labs_dictionary[x])
    })

    args <- c(args, mapps_final)
    structure(args, class='labels')
}

dfacets <- list(
    sex = setNames(c('Male', 'Female'), c('M', 'F')),
    comm = setNames(c('Fishing', 'Inland'), c('fishing', 'inland') )
)

age_breaks <- seq(15, 50, by=5)
