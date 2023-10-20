#################
#      AIMS     #
#################

# Read UNAIDS data from:
# https://www.unaids.org/en/resources/documents/2023/HIV_estimates_with_uncertainty_bounds_1990-present

#################
#   Libraries   #
#################

library(data.table)
library(ggplot2)
library(readxl)


# paths 
gitdir <- here::here()
source(file.path(gitdir, "R/base_plots.R"))
source(file.path(gitdir, "R/base_utilities.R"))

url <- 'https://www.unaids.org/en/resources/documents/2023/HIV_estimates_with_uncertainty_bounds_1990-present'
url_xlsx <- "http://www.unaids.org/sites/default/files/media_asset/HIV_estimates_from_1990-to-present.xlsx"

sprintf(" curl %s", url_xlsx) |> 
    system()

path_xlsx <- "~/Downloads/HIV_estimates_from_1990-to-present.xlsx"

##################
#      Main      #
##################

unaids_data <- readxl::read_xlsx(path_xlsx, sheet="HIV-Test-&-Treat_ByArea")
unaids_data <- unaids_data[ -(1:3), ]

class <- "Percentage of people living with HIV with suppressed viral load"
(as.character(unaids_data[1, ]) == class ) |>
    any(na.rm=TRUE) |>
    stopifnot("class of interest not found"=_)

idx <- which(as.character(unaids_data[1, ]) == class )

unaids_ss <- unaids_data[ , c(1:3, idx:(idx+14)) ]
fill_na_with_latest_char <- function(x){
    filler <- x[1]
    # stopifnot( ! is.na(filler))
    for( i in seq_along(x)){
        if (i == 1) 
            next 
        if( is.na(x[i]) ){
            x[i] <- filler
        }else{
            filler <- x[i]
        }
    }
    return(x)
}
unaids_ss[1,] <- fill_na_with_latest_char(unaids_ss[1, ])
unaids_ss[2,] <- fill_na_with_latest_char(unaids_ss[2, ])
x <- as.character(unaids_ss[1,]) 
fill_na_with_latest_char(x)
names(unaids_ss) <- paste(unaids_ss[2,], unaids_ss[3, ])
names(unaids_ss)[1:3] <- c("year", "code", "country")
unaids_ss <- unaids_ss[-(1:3), ] 
unaids_ss <- reshape2::melt( unaids_ss, id.vars = c("year", "code", "country"))
setDT(unaids_ss)
unaids_ss[,`:=`(
    age = gsub("^.*\\((.*)\\) (.*)$", "\\1", variable ),
    estimate = gsub("^.*\\((.*)\\) (.*)$", "\\2", variable ),
    variable = NULL,
    value = as.numeric(value)
)]
unaids_ss <- dcast(unaids_ss, year + code + country + age ~ estimate, value.var = 'value')
unaids_ss[, lapply(.SD, range,na.rm=TRUE), .SDcols=c('Estimate', 'Low', 'High')]

dplot <- subset(unaids_ss, 
    (! is.na(Estimate) & ! is.na(High) & ! is.na(Low)) &
    ( age %like% "Men|Women")
)

make_gg <- function(DT){
    ggplot(DT[! is.na(code)], aes(x = country, y = Estimate, ymin=Low, ymax=High, fill=age)) +
        geom_col(position = position_dodge(1)) +
        # geom_linerange(position = position_dodge(1), alpha=.2) +
        geom_hline( yintercept = 86 , linetype = 'dashed', color='grey80') +
        theme_bw() +
        theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1), legend.position = "bottom" ) +
        scale_fill_manual(values=c("lightblue", "pink")) +
        scale_y_continuous(expand=expansion(mult=c(0, .1))) +
        labs( x=NULL, y=paste(class, "\nas reported in UNAIDS 2023 data"), fill=NULL)
}

make_gg(dplot)
make_gg(dplot[High > 80])
     
ssa_countries <- c(
    "Angola",
    "Botswana",
    "Burundi",
    "Gabon",
    "Democratic Republic of Congo",
    "Republic of Congo",
    "Equatorial Guinea",
    "Eswatini",
    "Malawi",
    "Namibia",
    "Kenya",
    "Lesotho",
    "Rwanda", 
    "Uganda",
    "South Africa",
    "Tanzania",
    "Zambia",
    "Zimbabwe" 
)

p1 <- make_gg(dplot[ country %in% ssa_countries])
cmd <- ggsave2(p=p1, w=15, h=10,
    file = "UNAIDS_suppression_2023.pdf",
    LALA="/home/andrea/HPC/ab1820/home/projects/2022/longvl/figures"
)

system(cmd) 
system(zathura2gthumb(cmd))

