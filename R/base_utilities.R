###################
# Generic helpers #
###################

# non-interactive debugging (https://adv-r.hadley.nz/debugging.html#non-interactive-debugging)
dump_and_quit <- function(){ # works in concomitance with 
    # options(error=dump_and_quit)

    # save debugging info to file last.dump.rda
    dump.frames(to.file=TRUE)
    # quit R with error status

    paste(
        "Load dump file in an interactive session as follows:",
        "load('last.dump.rda')",
        "debugger()",
        sep='\n'
    ) |> sprintf("")

    q(status=1)
}

# not sure whether I should keep this here or in the scripts
if( ! interactive() ){
    options(error=dump_and_quit)
}

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


list.files.from.output.directory <- function(pattern, args=args, dir=NA_character_, rounds=NULL )
{
    # extracts files from the output directory specified by args$vl and args$jobname

    if(is.na(dir)){
        dir <- file.path(indir, make.suffix(args) )
    }

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

na2string <- function(x, str){
    if(is.numeric(x)){return(x)}
    if(is.factor(x)){ return(as.factor(str))}
    x[is.na(x)] <- str; x
}

`%which.like%` <- function(x, rgx)
    x[x %like% rgx]

`%which.not.like%` <- function(x, rgx)
    x[ ! x %like% rgx]

quantile2 <- function(x, ps = c(CL = .025, IL = .25, M = .5, IU = .75, CU = .975)) {
    # the posterior package has the same function name!
    out <- quantile(x = x, probs = ps, names = TRUE) |> t()
    if (!is.null(names(ps))) {
        colnames(out) <- names(ps)
    }
    as.data.table(out)
}

setnamesdict <- function(DT, dict) {
    setnames(DT, names(dict), unname(dict), skip_absent = TRUE)
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
    with(fread(path.community.idx), {
        dict <- as.character(comm_id)
        names(dict) <- as.character(comm_num)
        dict <<- dict
    } )
    dict[as.character(comm_num)]
}

split.agegroup <- function(x, breaks=c(15, 20, 25, 30, 35, 40, 45, 50))
{
    n_breaks <- length(breaks)
    labels <- paste(breaks[-n_breaks], breaks[-1] - 1, sep='-' )
    cut(x, breaks=breaks, labels=labels, include.lowest = TRUE, right=FALSE)
}

make.suffix <- function(args, cmdstan=FALSE)
{
    .alpha <- sprintf("alpha%.2f", args$stan.alpha)
    .alpha <- as.character(.alpha) |> gsub(pattern="\\.", replacement="", x=_)
    if(args$shared.hyper == TRUE) .alpha <- paste0(.alpha, 'sharedhyper')
    suffix <- c("", "cmdstan")[as.integer(cmdstan) + 1]
    suffix <- paste(suffix, .alpha, "vl", VIREMIC_VIRAL_LOAD, sep="_")
    suffix <- fifelse(args$only.firstparticipants==TRUE, 
        yes=paste0( suffix,'_firstpart'),
        no=suffix)
    jobname <- gsub('vl_[0-9]+|cmdstan|stan|_firstpart', '', args$jobname) |>
        gsub(pattern="__", replacement="_", x=_) |>
        gsub(pattern="_$|^|", replacement="", x=_)
    suffix <- fifelse(is.na(args$jobname), 
        yes=suffix, 
        no=paste0( jobname, '_',suffix))
    return(suffix)
}

fetch.args.from.suffix <- function(suffix, asDT=FALSE)
{
    stopifnot("fetch.args.from.suffix can only process one suffix"=length(suffix) == 1)
    outargs <- list(
        VL=gsub('^.*vl_([0-9]+).*$', '\\1', suffix) |> as.integer(),
        FTP=grepl('firstpart',suffix),
        JOB = NA_character_
    )

    jobname <- gsub("vl_[0-9]+|_firstpart",'',suffix)
    if(jobname != "")
        outargs$JOB <- jobname
    outargs$JOB<- gsub('_$', '', outargs$JOB)

    if(asDT)
        outargs <- as.data.table(outargs)

    return(outargs)
}

zathura2gthumb <- function(cmd){
    cmd <- sub("zathura", "gthumb", cmd)
    cmd <- sub(".pdf", ".png", cmd)
    return(cmd)
}

tex2png <- function(filename){
    gsub(".tex$", ".png", filename)
}

tex2pdf <- function(filename){
    gsub(".tex$", ".pdf", filename)
}

pdf2png <- function(filename){
    gsub(".pdf$", ".png", filename)
}

# copies data.frame to clipboard for easy copy-paste to google docs.
clipboard <- function(x, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE){
     con <- pipe("xclip -selection clipboard -i", open="w")
     write.table(x, con, sep=sep, row.names=row.names, col.names=col.names, quote=quote)
     close(con)
}

negate.percent.quantiles <- function(DT){
    DT[, `:=` (M = 1-M, CL=1-CU, CU=1-CL)]
    if("IL" %in% names(DT)){
        DT[, `:=` (IL = 1 - IU, IU=1-IL)]
    }
    return(DT)
}
