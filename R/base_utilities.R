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


list.files.from.output.directory <- function(pattern, args=args, indir=args$indir, vl=args$viremic.viral.load, jobname=args$jobname, rounds=NULL )
{
    # extracts files from the output directory specified by args$vl and args$jobname

    dir <- file.path(indir, make.suffix(args) )
    # dir <- file.path(indir, paste0('vl_', vl))
    files <- list.files( dir, pattern=pattern, full.names = TRUE, recursive = TRUE)
    
    # subset to file containing roundnames in their basename 
    if(! is.null(rounds) )
    {
        labels <- paste0('round', rounds) |> paste(collapse = '|')
        files <- files[ basename(files) %like% labels ]
    }

    return(files)
}

catn <- function(..., file="", sep=" ", fill=FALSE, labels = NULL, append = FALSE){
    # wrapper around cat to highlight the output message
    cat('\n---', ..., '---\n' , file = file, sep = sep, fill = fill, labels = labels, append = append)
}

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

empty2naDT <- function(DT)
{
    idx <- DT[, sapply(.SD, is.character)]
    cols <- names(idx)[idx]
    DT[, (cols) := lapply(.SD, .empty2na), .SDcols = cols]
    DT
}

remove.columns.with.unique.entry <- function(DT)
{
    idx <- DT[, which(lapply(.SD, uniqueN) == 1)]
    cols <- names(idx)
    cat('Removing constant columns, as uninformative:\n')
    DT[, lapply(.SD,unique) , .SDcols=cols]
    DT[, (cols):=NULL]
    return(DT)
}

na2zero <- function(x){
    x[is.na(x)] <-0; x
}

`%which.like%` <- function(x, rgx)
    x[x %like% rgx]

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

get.dcomm <- function(split.inland=TRUE)
{
    dcomm <- fread(path.community.types)

    .f <- function(x)
        fifelse(x %like% '^f', yes='fishing', no='inland')

    if(split.inland){
        .f <- function(x)
            fcase(x %like% '^f', 'fishing', x %like% '^t', 'trading', x %like% '^a', 'agrarian')
    }

    dcomm[, TYPE:=.f(COMM_NUM_A)]
    dcomm[, COMM_NUM:=as.integer(COMM_NUM_RAW)]
    dcomm
}

get.dcomm.idx <- function(comm_num, split.inland=TRUE)
{
    dcomm2 <- fread(path.community.idx) |> setkey('comm_num')
    dcomm2[data.table(comm_num=comm_num)]$comm_id
}

split.agegroup <- function(x)
{
    breaks <- c(15, 20, 25, 30, 35, 40, 45, 50)
    n_breaks <- length(breaks)
    labels <- paste(breaks[-n_breaks], breaks[-1] - 1, sep='-' )

    cut(x, breaks=breaks, labels=labels, include.lowest = TRUE, right=FALSE)
}

make.suffix <- function(args)
{
    suffix <- paste0("vl_", VIREMIC_VIRAL_LOAD)
    suffix <- fifelse(args$only.firstparticipants==TRUE, 
        yes=paste0( suffix,'_firstpart'),
        no=suffix)
    suffix <- fifelse(is.na(args$jobname), 
        yes=suffix, 
        no=paste0( suffix,'_',args$jobname))
    return(suffix)
}

fetch.args.from.suffix <- function(suffix, asDT=FALSE)
{
    stopifnot("fetch.args.from.suffix can only process one suffix"=length(suffix) == 1)
    outargs <- list(
        VL=gsub('vl_([0-9]+).*$', '\\1', suffix) |> as.integer(),
        FTP=grepl('firstpart',suffix),
        JOB = NA_character_
    )

    jobname <- gsub("vl_[0-9]+|_firstpart",'',suffix)
    if(jobname != "")
        outargs$JOB <- jobname

    if(asDT)
        outargs <- as.data.table(outargs)

    return(outargs)
}
