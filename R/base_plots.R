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
theme_set( theme_default() )


#########################
# Requirements + saving #
#########################

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
        strip.background = element_rect(colour="white", fill="white"))
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
# specify palettes #
####################

palettes <- list(

    sex=c(
        F="#F4B5BD",
        M="#85D4E3",
        Female="#F4B5BD",
        Male="#85D4E3"
        ),

    comm = c(
        fishing="#85D4E3",
        inland="#9C964A"),

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
        "#85D4E3",
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

    round = c(
        `16` = ,
        `17` = ,
        `18` = ,
        `19` = ,
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

    if("SEX" %in% nms)
        DT[, SEX_LAB :=  sex_dictionary[SEX]]

    if("SEX_LABEL" %in% nms)
        DT[, SEX_LAB :=  sex_dictionary[SEX_LABEL]]
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


prettify_labels <- function(DT)
{
    # preforms all of the above
    nms <- names(DT)

    DT <- prettify_sex(DT)
    DT <- prettify_round(DT)
    DT
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


