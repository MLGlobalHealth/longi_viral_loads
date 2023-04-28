###################
# Generic helpers #
###################


# saving files
# ____________

savefile <- function(data, filename, overwrite=FALSE)
{
    # checks whether filename exists,
    # saves file, and returns a string to open that 
    # file from the cmd line

    nvim_cmd_to_open_file <- paste('! xdg-open', filename, '&')
    

    # check if file exists and whether to overwrite it
    if(file.exists(filename) & !overwrite)
    {
        cat("File", filename, "already exists")
        if( ! overwrite)
        {
            cat(" Not overwriting it.\n")
            return(nvim_cmd_to_open_file)
        }
    }

    # Proceed to save based on file type. 
    success_msg <- paste0("File",filename," successfully saved\n")
    if( grepl( '.csv$', filename ))
    {
        fwrite(x=data, file=filename)
        cat(success_msg)
    }else if(grepl( '.rds$', filename ))
    {
        saveRDS(object=data,file =filename)
        cat(success_msg)
    }

    # 
    return(nvim_cmd_to_open_file)
}


list.files.from.output.directory <- function(pattern, indir=args$indir, vl=args$viremic.viral.load, jobname=args$jobname, rounds=NULL )
{
    # extracts files from the output directory specified by args$vl and args$jobname
    dir <- file.path(indir, paste0('vl_', vl))
    files <- list.files( dir, pattern=pattern, full.names = TRUE, recursive = TRUE)
    
    # subset to file containing roundRR in their basename 
    if(! is.null(rounds) )
    {
        labels <- paste0('round', rounds) |> paste(collapse = '|')
        files <- files[ basename(files) %like% labels ]
    }

    return(files)
}

catn <- function(x) cat('\n----', x, '----\n')

empty2NA <- function(DT)
{
    toNA <- function(vec) 
    {
        vec[vec==""] <- NA
        return(vec)
    }
    DT[, lapply(.SD, toNA), ]
}

.empty2na <- function(x)
{
    y <- gsub(' ', '', x)
    if(any(is.na(x)) & length(y=='')>0) warning('Vector x already contains NA values')
    x[y==''] <- NA
    x
}

na2zero <- function(x){
    x[is.na(x)] <-0; x
}

############################
# Project-specific helpers #
############################

round2numeric <- function(x) 
{
    x <- gsub('S$','.5',x)
    as.numeric(gsub( '[A-z]', '', x))
}

round2factor <- function(x)
{
    stopifnot(is.numeric(x))
    x <- paste0('R', x)
    x <- gsub('\\.5', 'S', x)
    as.factor(x)
}


commcode2commtype <- function(x)
{
    require(data.table)
    out <- fcase(
        x %like% '^f', 'fishing',
        x %like% '^a', 'agrarian',
        x %like% '^t', 'trading')
    return(out)
}
