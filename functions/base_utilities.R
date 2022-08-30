theme_set(theme_bw())
geom.text.size = 5/14*10
update_geom_defaults("text", list(size = geom.text.size))
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
                 supp_status =  c('#8F190E', '#DBCB00', '#00DB25')
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
