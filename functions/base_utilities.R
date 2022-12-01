theme_set(theme_bw())
geom.text.size = 5/14*7
update_geom_defaults("text", list(size = geom.text.size))
# update_geom_defaults("text", list(size = 3))
GeomLabel$default_aes$family <- 'sans'

p_reqs <- theme(
                plot.title=element_text(size=12, family='sans'),
                text=element_text(size=10, family = 'sans'),
                axis.title=element_text(size=10, family='sans'), 
                axis.text=element_text(size=7, family = 'sans')
)

palettes <- list(
                 sex = wesanderson::wes_palette("Moonrise3")[c(2,1)],
                 comm = wesanderson::wes_palette("Moonrise3")[c(1,3)],
                 arvmed = c('#E03531', '#51B364'),
                 cuarvmed = c('#E03531', '#51B364'),
                 supp_status =  c('#8F190E', '#DBCB00', '#00DB25'),
                 comm2 = c("#9C964A", "#85D4E3","#FF9933")
)

ggsave2 <- function(p, file, w, h, LALA=vl.out.dir, u='in')
{
        # I used LALA because I think there may be some naming issues
        filename <- file
        filename2 <- gsub('pdf$','png',filename)
        cat('Saving', filename, '...\n')
        ggsave(p, filename=file.path(LALA,filename2), width=w, height=h, units=u)
        ggsave(p, filename=file.path(LALA, filename), width=w, height=h, units=u)	
}

naturemed_reqs <- function() 
{
    # call this before doing your plots
    nm_reqs <<- theme(axis.text = element_text(size=5, family='sans'),
                      text=element_text(size=7,family='sans'),
                      legend.text=element_text(size=7, family='sans'))
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
}
