palettes <- list(
                 sex = wesanderson::wes_palette("Moonrise3")[c(2,1)],
                 comm = wesanderson::wes_palette("Moonrise3")[c(1,3)],
                 arvmed = c('#E03531', '#51B364'),
                 cuarvmed = c('#E03531', '#51B364')
)

saveplot <- function(filename, plot, w, h, ...)
{
        stopifnot(dir.exists(dirname(filename)))
        ggsave(filename=filename, plot=plot, width=w, height=h, ...)
}

