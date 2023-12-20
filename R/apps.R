#' @importFrom Rcpp sourceCpp
#' @importFrom rhdf5 h5read
#' @useDynLib GxG
NULL

## ================== Custom instantiators and converters for gMatrix ================== ##

#' @name homeology
#' Given a (named) character vector or XStringSet computes a matrix of string distances between
#' all tiles of width + pad.  Returns a gMatrix whose seqnames are the names(seq) and coordinates
#' refer to the provided sequences.  Every GRanges pair ij of the gMatrix represents the edit distance
#' (default levenshtein) between the sequences in the provided sequence tiles (+/- some padding).
#'
#' If flag rc=TRUE, then the distance ij will be the edit distance between the sequence in GRanges i and the reverse
#' complement of sequence in GRanges j. 
#'
#' @param seqs a character vector or XStringSet (e.g. DNAStringSet)
#' @param gr GRanges
#' @param stride specifies both the width and the stride of the sequence bins / tiles around which to measure homeology and which are returned in the output gMatrix
#' @param pad the padding around each tile with which to measure sequence homeology
#' @param rc logical flag specifying whether to provide reverse commplement distance
#' @author M arcin Imielinski
#' @export
#' @return new gPair containing the difference between x and y
homeology = function(seqs, gr = NULL, stride = 1, pad = 10, rc = FALSE, verbose = FALSE, method = "levenshtein", ignoreCase = TRUE, ...)
{
  ## short hand for supplying GRanges
  if (is(seqs, 'GRanges'))
  {
    if (is.null(seqs$seq))
      stop('If GRanges provided to homeology, they must have character or DNAStringSet field $seq that matches the width of the corresponding GRanges')

    gr = seqs
    seqs = seqs$seq
    names(seqs) = seqnames(gr)
  }
  
  if (is.null(names(seqs)))
    names(seqs) = 1:length(seqs)

  if (is.character(seqs))
  {
    seqs = DNAStringSet(seqs)
  }

  if (is.null(gr))
    gr = GRanges(names(seqs), IRanges(1, width(seqs)))
  else if (!identical(width(gr), width(seqs)) || !identical(as.character(seqnames(gr)), names(seqs)))
    stop('If gr is provided, gr must be same length as seqs and each gr[i] must have the same width and seqnames as the nchar and name of the corresponding seqs[i]')

  gr = gr.fix(gr)
  if (any(ix <- strand(gr)=='*'))
    strand(gr)[ix] = '+'

  width = stride
  starts = start(gr)
  names(starts) = seqnames(gr)

  grs = gr.fix(GRanges(seqnames(gr), IRanges(1, width(gr))))
  grs$sn.new = dedup(as.character(seqnames(grs)))
  tiles.og  = suppressWarnings(gr.tile(grs, width))
  strand(tiles.og) = strand(gr)[tiles.og$query.id]
  tiles = suppressWarnings(tiles.og + pad)
  end(tiles) =  suppressWarnings(pmin(seqlengths(tiles)[as.character(seqnames(tiles))], end(tiles)))
  start(tiles) = suppressWarnings(pmax(1, start(tiles)))

  ## need to "dedup" seqnames in case the input granges are on the same chromosome
  tiles = GRanges(grs$sn.new[tiles$query.id], ranges(tiles), strand = strand(tiles))
  names(seqs) = dedup(names(seqs))

  if (rc)
  {
    tiles = c(tiles, gr.flipstrand(tiles))
  }

  if (verbose)
    {
      message('Populating ', length(tiles), ' bins with width ', width(tiles)[1], ' and stride ', width)
    }

  tiles$seq = BSgenome::getSeq(seqs, tiles)

  if (verbose)
  {
    message('Computing ', length(tiles), ' by ', length(tiles), ' distance matrix')
  }

  D = as.matrix(Biostrings::stringDist(tiles$seq, method = method, ignoreCase = ignoreCase, upper = TRUE, ...))

  if (rc) ## only keep +- ij pairs
  {
    D = D[which(as.logical(strand(tiles)=='+')), which(as.logical(strand(tiles)=='-'))]
  }

  ## case gMatrix to original (nonoverlapping) tiles (ie without the padding)
  out = GRanges(as.character(seqnames(tiles.og)),
                IRanges(start(tiles.og) + starts[tiles.og$query.id]-1,
                        end(tiles.og) + starts[tiles.og$query.id]-1),
                strand = strand(tiles.og))
  out$seq = tiles.og$seq
  return(gM(out, D))
}


#' @name straw
#' @title straw
#' @description
#' Instantiates gMatrix from .hic object using straw API https://github.com/theaidenlab/straw/tree/master/R
#' used to extract all of the data in length n query gr (note will do n choose 2 queries since
#' straw only supports pairwise queries)
#'
#' @keywords straw
#' @param hic path to .hic file
#' @param gr granges to query
#' @param norm string specifying normalization to apply ("KR" (default), "NONE", VC")
#' @param res resolution to query (default 10kb)
#' @param type "BP" vs "FRAG"
#' @param mc.cores parallelization to apply for query (useful when length(gr)>2)
#' @return gMatrix 
#' @export
#' @author Marcin Imielinski
straw = function(hic, gr = NULL, norm = "NONE", type = 'BP', res = 1e4, mc.cores = 1, colormap = c('white', 'red', 'black'), ...)
{
  hic = normalizePath(hic)
  if (is.null(gr))
    gr = hg_seqlengths()

  if (is.integer(gr) | is.numeric(gr))
    {
      if (is.null(names(gr))) ## assume this is numeric form of chromosomes
        {
          gr = as.character(gr)
        }
      else ## assume it's a seqlengths
        {
          gr = si2gr(gr)
        }
    }
  
  if (is.character(gr)) {
      gr_character = TRUE #added this to fix differences in seqlengths and ultimately different data from cooler-other line at end to fix it
      gr = parse.gr(gr, seqlengths = hg_seqlengths())
  } else {gr_character=FALSE}
  if (!inherits(gr, 'GRanges'))
      gr = si2gr(gr)
  gr = reduce(gr.stripstrand(gr[, c()]))
  gr.dt = as.data.table(gr)
  levels.seqnames = c(1:22,"X","Y","M")
    gr.dt$seqnames = factor(as.character(gr.dt$seqnames),levels=levels.seqnames) #seqnames is a factor so need to do as.character
  gr.dt = gr.dt[order(seqnames),]
#  grs = as.data.table(gr)[, paste(seqnames, start, end, sep = ":")]
  grs = gr.dt[, paste(seqnames, start-1, end, sep = ":")] #added -1 to start because .hic is 0 based and it doesn't return first bin
  n = length(gr)
  grs.pairs = chr2.pairs = chr1.pairs = NULL
  grs.singles = paste(grs, grs)
  ## new new new
  r1 = grs
  r2 = grs
  chr1.singles = chr2.singles = as.character(seqnames(gr))
  if (n>1)
    {
      pairs = t(combn(n, 2)) ##pairs
      grs.pairs = paste(grs[pairs[,1]], grs[pairs[,2]])   ## combinations
      #' Tuesday, Nov 07, 2017 04:43:52 PM - Julie: adding as.character()
      r1 = c(r1, grs[pairs[,1]])
      r2 = c(r2, grs[pairs[,2]])
      chr1.pairs = as.character(seqnames(gr)[pairs[,1]])
      chr2.pairs = as.character(seqnames(gr)[pairs[,2]])
    }
  str = paste(norm, hic, c(grs.singles, grs.pairs), type, as.integer(res))
  chr1 = as.character(c(chr1.singles, chr1.pairs))
  chr2 = as.character(c(chr2.singles, chr2.pairs))
  ## use strawR
  out = rbindlist(mcmapply(r1 = r1,
                           r2 = r2,
                           FUN = function(r1, r2)
                           {
                             dt = as.data.table(
                                       strawr::straw(norm = norm,
                                                     unit = type,
                                                     binsize = res,
                                                     fname = hic,
                                                     chr1loc = r1,
                                                     chr2loc = r2))
                             setnames(dt, c("start1", "start2", "counts"))
                             dt[, start1 := start1 +1] #adding 1 back since .hic files are zero based
                             dt[, start2 := start2 +1] #adding 1 back since .hic files are zero based
                             dt[, chr1 := gsub("^([0-9XY]+):.*", "\\1", r1)]
                             dt[, chr2 := gsub("^([0-9XY]+):.*", "\\1", r2)]
                             dt[, end1 := start1 + res - 1]
                             dt[, end2 := start2 + res - 1]
                             return(dt)
                           }, SIMPLIFY = FALSE,  mc.cores = mc.cores), fill = TRUE)  
  ## out = rbindlist(mcmapply(str = str, chr1 = chr1, chr2 = chr2,
  ##                          FUN = function(str, chr1, chr2)
  ##                          {
  ##                            dt = as.data.table(
  ##                              tryCatch(GxG:::straw_R(str),
  ##                                       error = function (e)
  ##                                         data.table()))
  ##                            return(dt)
  ##                          }, SIMPLIFY = FALSE,  mc.cores = mc.cores), fill = TRUE)
  out = out[!is.na(counts), ]
  if (!nrow(out))
    stop('Query resulted in no output, please check .hic file or input coordinates')

  out[, end1 := as.integer(start1+res-1)] #added as integer because was reading as 1e* which did not work with mapping below
  out[, end2 := as.integer(start2+res-1)]
  out[, str1 := paste0(chr1, ':', start1, '-', end1)]
  out[, str2 := paste0(chr2, ':', start2, '-', end2)]
  gr.out = unique(dt2gr(rbind(out[, .(seqnames = chr1, start = start1, end = end1)],
                              out[, .(seqnames = chr2, start = start2, end = end2)])))
  gr.map = data.table(grs = gr.string(gr.out), ix = 1:length(gr.out), key = 'grs')
  out$i1 = gr.map[.(out$str1), ix]
  out$j1 = gr.map[.(out$str2), ix]  
  out[, i := pmin(i1, j1)]
  out[, j := pmax(i1, j1)]
  out=out[i > 0 & i <= length(gr.out) & j > 0 & j <= length(gr.out)]
  out = out[order(i,j),]
  if(gr_character==TRUE)
      gr.out = suppressWarnings(GRanges(as.data.table(gr.out),seqlengths= hg_seqlengths())) %>% trim()
  gm = gM(gr.out, out[, .(i, j, value = counts)])
  return(gm) 
}


#' @name hicpro
#' @title hicpro
#' @description
#' Instantiates gMatrix from .matrix and .bed Hi-C pro file 
#'
#' @keywords hicpro
#' @param mat path to .matrix file
#' @param bed path to .bed file
#' @return gMatrix
#' @export
#' @author Marcin Imielinski
hicpro = function(mat, bed = NULL)
{
  if (is.null(bed))
    bed = paste0(gsub('.matrix$', '', mat), '_ord.bed')

  dat = fread(mat)
  setnames(dat, c('i', 'j', 'value'))
  
  if (!file.exists(bed))
    stop('bed file must be provided or one must exist with the same prefix as the .matrix file suffixed with _ord.bed.')
 
  gr = gr.fix(rtracklayer::import(bed))

  if (max(dat$i, dat$j)>length(gr))
    stop('provided .matrix has out of bounds indices relative to provided .bed file, please check to make sure these are compatible')
  return(gM(gr, dat))
}


#' @name cooler
#' @title cooler
#' @description
#' Instantiates gMatrix from .cool and/or .mcool file, +/- at specific locations. 
#'
#' To find viable resolutions (e.g. for mcool files) use info = TRUE and you will get a list
#' of info, including seqlengths. 
#'
#' @keywords cooler
#' @param file path to .cool or m.cool file
#' @param gr optional GRanges of intervals to query
#' @param res optional resolution to pull
#' @param info logical flag (FALSE) specifying whether to return a vector of file info instead of pulling down data, this will include available resolutions (for .mcool format) and fields
#' @return gMatrix or list
#' @export
#' @author Marcin Imielinski, Stanley Clarke

cooler = function(file, res = NULL, gr = NULL, info = FALSE, mc.cores = 1) {
    if (info == TRUE) {
        return(cooler_resolutions(file))
    } else {
                                        # Try to get bins; if it fails, run cooler_resolutions
        tryCatch({
            hic = normalizePath(file)
            if (is.null(gr))
                gr = hg_seqlengths()
            if (is.integer(gr) | is.numeric(gr)) {
                if (is.null(names(gr))) { ## assume this is numeric form of chromosomes
                    gr = as.character(gr) } else { ## assume it's a seqlengths
                    gr = si2gr(gr) }
            }
            if (is.character(gr))
                gr = parse.gr(gr, seqlengths = hg_seqlengths())
            if (!inherits(gr, 'GRanges'))
                gr = si2gr(gr)
            gr = reduce(gr.stripstrand(gr[, c()]))
            anchor = getBins(file, as.integer(res))
            anchor$i = 1:length(anchor)
            #intersect to see whether the whole genome is in the cooler file
            inter_gr = anchor %&% gr
            if ((length(inter_gr) > 0)) {
                gr_sub = gr %Q% (seqnames %in% seqnames(anchor)@values)
                if(length(gr_sub) < length(gr))
                    warning("There are seqlengths specified that are not in the cooler file ")
                sliced = getSlice(anchors = anchor, file = file, res = as.integer(res), gr = gr_sub, mc.cores = 1)
                dat = sliced[[1]]
                dat[,value := as.numeric(value)] #to stay consistent with straw function
                dat = dat[order(i,j),]
                gmat = gM(dat = dat, gr = sliced[[2]])
                return(gmat)
            } else {
                stop("The specified gr does not have sequences in common with the cooler file")
            }
        },
        error = function() {
                                        # If getBins fails, run cooler_resolutions and return its output
            cooler_resolution = cooler_resolutions(file)
            message("The specified resolution is not present in the cooler file. The resolutions present are:")
                                        #stop(print(cooler_resolution))
            print(cooler_resolution)
        })
    }
}


## cooler = function(file, gr = NULL, res = NULL, info = FALSE, field = NULL, nochr = FALSE, mc.cores = 1)
## {
##   contents = subfile = NULL

##   ## if mcool need to pick a resolution
##   if (grepl('mcool$', file))

##   {
##     contents = as.data.table(rhdf5::h5ls(file))[grepl('\\/resolutions\\/\\d+/bins$', group), ][!(name %in% c('chrom', 'end', 'start')), ]
##     contents[, subfile := gsub('\\/bins$', '', group)]
##     contents[, resolution := gsub('(\\/resolutions\\/)|(\\/bins)', '', group) %>% as.integer]
##     contents = contents[order(resolution), .(field = name, dim = dim, resolution)]

##     if (is.null(res))
##       res = contents$resolution[1]

##     if (!(res %in% contents$resolution))
##       stop(sprintf('resolution not present in file %s: choose one of these available resolutions: %s or use default (%s)', file, paste(contents$resolution %>% unique, collapse = ', '), contents$resolution[1]))

##     contents = contents[resolution %in% res, ]

##     subfile = paste('resolutions', contents[, resolution[1]], sep = '/')

##     ## not sure how to get any field other than count .. eg VC, VC_SQRT
##     ## (my guess is these are supposed to get computed on the fly from the marginal
##     ## correction factors rather than stored in the file)
##     ##
##     ## if (!(field %in% contents$field))
##     ##   stop(sprintf('resolution not present in file %s: choose one of these available resolutions: %s or use default (%s)', file, paste(contents$resolution %>% unique, collapse = ', '), contents$resolution[1]))

##     ## subfile = paste('resolutions', contents[field == field, resolution[1]], sep = '/')
##   }

##   ## install cooler if not installed on reticulate python virtualenv
##   reticulate::py_run_string("try:
##     import cooler
##     cooler_installed = True
## except ImportError:
##     cooler_installed = False
## ")
##   if (!reticulate::py$cooler_installed)
##     reticulate::py_install("cooler")

##   ## construct subfile (ie if we have mcool)
##   if (!is.null(subfile))
##   {
##     uri = paste(file, subfile, sep = '::')
##   } else
##   {
##     uri = file
##   }

##   ## load data via cooler
##   reticulate::py_run_string(sprintf('data = cooler.Cooler("%s")', uri))

##   ## get seqlengths from object
##   tmp = reticulate::py_run_string('chromsizes = data.chromsizes')$chromsizes
##   sl = structure(as.vector(tmp), names = names(tmp))

##   if (info)
##   {
##     info = reticulate::py_run_string(sprintf('info = data.info'))$info
##     if (!is.null(contents))
##       {
##         info$resolutions = contents$resolution %>% unique
##         info$default.resolution = contents$resolution[1]
##       }
##     info$fields = colnames(reticulate::py_run_string(sprintf('pix = data.pixels()[:1]'))$pix)[-c(1:2)]
##     info$seqlengths = sl
##     return(info)
##   }

##   ## create our bin GRanges
##   ##  bins = reticulate::py_run_string('bins = data.bins()[:]; print(bins)')$bins

##   ## takes care of weird reticulate conversion error
##   bins = reticulate::py_run_string('bins = data.bins()[:]; bins.chrom= bins.chrom.astype("str")')$bins %>% as.data.table
##   bins[, start := start+1]
##   setnames(bins, 'chrom', 'seqnames')
##   bins = bins %>% dt2gr(seqlengths = sl)

##   if (nochr)
##     bins = gr.sub(bins, 'chr', '')

##   if (!is.null(gr))
##   {
##     gr.old = gr
##     gr = gr.fix(gr, sl, drop = TRUE) %>% gr.stripstrand
##     ## more fixing
##     end(gr) = pmax(end(gr), 1)
##     start(gr) = pmax(start(gr), 1)
##     if (length(gr)<length(gr.old))
##       warning("Some of the provided ranges had to be dropped / clipped to make compatible with the seqlengths of the .cool / .mcool file")

##     if (!length(gr))
##       return(gM(GRanges(seqlengths = sl)))

##     pix = expand.grid(i = 1:length(gr), j = 1:length(gr))

##     res = mclapply(1:nrow(pix), function(k)
##       reticulate::py_run_string(sprintf('res = data.matrix(balance = False, join = False, as_pixels = True).fetch("%s", "%s")',
##                                         gr[pix$i[k]] %>% gr.string,
##                                         gr[pix$j[k]] %>% gr.string))$res %>% as.data.table, mc.cores = mc.cores) %>% rbindlist
##                                         #  res = unique(res, by = c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'))
##   }
##   else
##   {
##     res = reticulate::py_run_string(sprintf('res = data.matrix(balance = False, join = False, as_pixels = True)[:]'))$res %>% as.data.table
##   }

##   if (is.null(field))
##     field = names(res)[-c(1:2)][1]

##   dat = data.table(i = res[[1]]+1, j = res[[2]]+1, value = res[[field]])[, .(i = pmin(i, j), j = pmax(i, j), value = value)]

##   gm = gM(bins, dat = dat)
##   return(gm)
## }


#' @name cooler_resolutions
#' @title cooler_resolutions
#' @description
#' get resolutions from a mcool file
#'
#' 
#' @keywords cooler_resolutions
#' @param file path to .cool or m.cool file
#' @return list of resolutions in cooler
#' @export
#' @author Marcin Imielinski, Stanley Clarke
cooler_resolutions = function(mcool) {
                                        # List the contents of the MCOOL file
    cooler_details.dt = rhdf5::h5ls(mcool) %>% as.data.table()
    suppressWarnings(cooler_details.dt[,name := as.integer(name)])
    resolutions = cooler_details.dt[!is.na(name),]$name %>% sort()
    print(resolutions)
}


#' @name getBins
#' @title getBins
#' @description
#' get the bins of the mcool file for a specific mcool file at a given resolution
#'
#' 
#' @keywords getBins
#' @param file path to .cool or m.cool file
#' @param res resolution to use that is in the mcool file
#' @return 
#' @author Stanley Clarke
getBins = function(file,res = NULL) {

    bins.group = ifelse(is.null(res),'/bins',paste('resolutions',res,'bins',sep='/'))
    
    #s.i = seqinfo.cool(file,res)
    chr.group <- ifelse(is.null(res),'/chroms',paste('resolutions',res,'chroms',sep='/'))
    s.i <- Seqinfo(seqnames = as.vector(rhdf5::h5read(file,name=paste(chr.group,'name',sep="/"))),
                   seqlengths = as.vector(rhdf5::h5read(file,name=paste(chr.group,'length',sep="/"))))

    a.chr = rhdf5::h5read(file,name=paste(bins.group,'chrom',sep="/"))
    ## Cooler use open ended starts, we will use close non-overlapping bins
    a.start = rhdf5::h5read(file,name=paste(bins.group,'start',sep="/"))+1
    a.end = rhdf5::h5read(file,name=paste(bins.group,'end',sep="/"))
    
    anchors = GRanges(a.chr,
                       IRanges(a.start,a.end),
                       seqinfo=s.i)

    return(anchors)
}

#' @name getSlice
#' @title getSlice
#' @description
#' get resolutions from a mcool file
#'
#' 
#' @keywords getSlice
#' @param anchors result from getBins
#' @param res resolution to use that is in the mcool file
#' @param file path to mcool
#' @param gr optional granges to get from mcool file
#' @param mc.cores optional increase the cores when extracting for multiple granges
#' @return list of resolutions in cooler
#' @author Stanley Clarke

getSlice = function(anchors,file,res,gr=NULL,mc.cores=1) {
################################################################################
    ## Get the chromosme bins as a GRanges
################################################################################
    indexes.group = ifelse(is.null(res),"/indexes",paste('resolutions',res,'indexes',sep='/'))
    chr.group = ifelse(is.null(res),"/chroms",paste('resolutions',res,'chroms',sep='/'))
    pixels.group = ifelse(is.null(res),"/pixels",paste('resolutions',res,'pixels',sep='/'))
    
    indexes = list(chr=as.vector(rhdf5::h5read(file,name=paste(chr.group,'name',sep="/"))),
                    chr_idx=as.vector(rhdf5::h5read(file,name=paste(indexes.group,'chrom_offset',sep="/"))),
                    bin1_idx=as.vector(rhdf5::h5read(file,name=paste(indexes.group,'bin1_offset',sep="/")))
                    )
    if(is.null(gr)) {
        seq.leng = hg_seqlengths() %>% as.data.table() %>% setNames("end")
        seq.leng[,start := 1][,seqnames := names(hg_seqlengths())]
        gr = GRanges(seq.leng)
        gr = gr %Q% (seqnames %in% c(1:22,"X","Y"))
        }
################################################################################
    ## Reading chromosome chunk If chr1 is null, return the full cool file
################################################################################
    if (!is.null(gr)) {
            chunks.lst = mclapply(1:length(gr), function(x) {
                chr.start.idx = min((anchors %&% gr[x])$i)
                chr.end.idx = max((anchors %&% gr[x])$i)
#                idx.chunk = seq(chr.start.idx,chr.end.idx)
                idx.chunk = seq(chr.start.idx,chr.end.idx+1)
                bin1.idx = as.vector(rhdf5::h5read(file,
                                                    name=paste(indexes.group,'bin1_offset',sep="/"),
                                                    index=list(idx.chunk)))
                slice1 = sum(bin1.idx[-1] - bin1.idx[-length(bin1.idx)])-1
                chunk  = seq(bin1.idx[1]+1,bin1.idx[1]+1+slice1)
                dt = data.table(bin1_id = as.vector(rhdf5::h5read(file,
                                                                   name=paste(pixels.group,'bin1_id',sep="/"),
                                                                   index=list(chunk)))+1, #had to add +1 to match straw and simulated data
                                 bin2_id = as.vector(rhdf5::h5read(file,
                                                                   name=paste(pixels.group,'bin2_id',sep="/"),
                                                                   index=list(chunk)))+1, #had to add +1 to match straw and simulated data
                                 count = as.vector(rhdf5::h5read(file,
                                                                 name=paste(pixels.group,'count',sep="/"),
                                                                 index=list(chunk)))
                                 )
                return(dt)
        },mc.cores=mc.cores)
            dt = rbindlist(chunks.lst)
        #clean up indices to work with gM since the datatable command adds all indices corresponding to that region
        names(dt) = c("i","j","value")
        anchors.sub = anchors %&% gr
        dt.sub = dt[i %in% anchors.sub$i & j %in% anchors.sub$i,]
#        anchors.sub$i = 1:length(anchors.sub)
        anchors.sub = anchors.sub %Q% (i %in% dt.sub$i | i %in% dt.sub$j)
        anchors.sub$newi = 1:length(anchors.sub)
        anchors.sub$newj = anchors.sub$newi
        index.dt = as.data.table(anchors.sub)[,.(i,newi,newj)]
        anchors.sub$i = anchors.sub$newi
        anchors.sub$newi = NULL
        anchors.sub$newj = NULL
        dt2 = merge.data.table(dt.sub,index.dt[,.(i,newi)],by="i")
        dt = merge.data.table(dt2,index.dt[,.(i,newj)],by.x="j",by.y="i")[,.(newi,newj,value)]
        names(dt) = c("i","j","value")
        dt = dt[order(i,j),]
        anchors=anchors.sub
        }
    ################################################################################
    ## Reading the chunks from the cool file
################################################################################
    if(is.null(gr)) {
        dt = data.table(bin1_id = as.vector(rhdf5::h5read(file,
                                                           name=paste(pixels.group,'bin1_id',sep="/"),
                                                           index=list(chunk)))+1,
                         bin2_id = as.vector(rhdf5::h5read(file,
                                                           name=paste(pixels.group,'bin2_id',sep="/"),
                                                           index=list(chunk)))+1,
                         count = as.vector(rhdf5::h5read(file,
                                                         name=paste(pixels.group,'count',sep="/"),
                                                         index=list(chunk)))
                         )
        names(dt) = c("i","j","value")
    }
    return(list(dt,anchors))
}
