################################################################################
# Utility functions to handle a GDS file
## Create a GDS node
#' @importFrom gdsfmt exist.gdsn addfolder.gdsn index.gdsn add.gdsn
.create_gdsn <- function(root_node,
                         target_node,
                         new_node,
                         val = NULL,
                         storage = "string32",
                         valdim = dim(val),
                         replace = TRUE,
                         attr = NULL,
                         is_folder = FALSE){

  # Check if the target node exists in the root node
  check <- exist.gdsn(node = root_node, path = target_node)
  if(!check){ stop(target_node, " does not exist!") }

  # Get the index of the target node
  target_node_index <- index.gdsn(node = root_node, path = target_node)

  # Create a folder node if specified
  if(is_folder){
    new_node_index <- addfolder.gdsn(node = target_node_index,
                                     name = new_node,
                                     replace = TRUE)
  } else {
    # Create a new node with the specified value, storage type, and dimensions
    new_node_index <- add.gdsn(node = target_node_index, name = new_node,
                               val = val, storage = storage, valdim = valdim,
                               replace = replace)
    # Assign attributes to the new node if provided
    .putAttrGdsn(node = new_node_index, attr = attr)
    # Compress the new node
    .compressNodes(node = new_node_index)
  }
  invisible(new_node_index)
}

# Assign attributes to a GDS node
#' @importFrom gdsfmt put.attr.gdsn
.putAttrGdsn <- function(node, attr){
  # Check if attributes are provided
  if(!is.null(attr)){
    # Loop through each attribute and assign it to the node
    for(i in seq_along(attr)){
      put.attr.gdsn(node = node,
                    name = names(attr)[i], val = attr[[i]])
    }
  }
}

# Compress GDS nodes
#' @importFrom gdsfmt index.gdsn compression.gdsn readmode.gdsn
.compressNodes <- function(node){
  # Loop through each node
  for(i in seq_along(node)){
    # Compress the node using ZIP_RA compression
    compression.gdsn(node = node, compress = "ZIP_RA")
    # Set the node to read mode
    readmode.gdsn(node = node)
  }
}

## Retrieve data from a GDS node
#' @importFrom gdsfmt index.gdsn read.gdsn readex.gdsn
.get_data <- function(object, node, start = NULL, count, sel = NULL){
  # Check if selection criteria (sel) is provided
  if(is.null(sel)){
    # If selection criteria is not provided
    if(is.null(start)){
      # If start is not provided, read the entire node
      out <- read.gdsn(node = index.gdsn(node = object$root, path = node))
    } else {
      # If start is provided, read the node from the specified start position and count
      out <- read.gdsn(node = index.gdsn(node = object$root, path = node),
                       start = start,
                       count = count)
    }
  } else {
    # If selection criteria is provided, read the node with the selection criteria
    out <- readex.gdsn(node = index.gdsn(node = object$root, path = node), sel = sel)
  }

  # Return the retrieved data
  return(out)
}

# Retrieve a scan result for a phenotype
.marker_annotation_table <- function(object) {
  cached <- object@lazydata$marker_table
  if (!is.null(cached)) {
    return(cached)
  }
  cached <- data.frame(
    variant_ID = getMarID(object = object),
    Chr = getChromosome(object = object),
    Pos = getPosition(object = object),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  object@lazydata$marker_table <- cached
  cached
}

.get_scan <- function(object, pheno_name, sel = NULL){
  scan_dat <- .store_read_scan_matrix(object, pheno_name, columns = sel)
  if (is.null(scan_dat)) {
    stop("No scan data for phenotype: ", pheno_name, call. = FALSE)
  }
  mar <- .marker_annotation_table(object = object)
  out <- cbind(mar, scan_dat$matrix)
  colnames(out) <- c("variant_ID", "Chr", "Pos", scan_dat$colnames)
  return(out)
}

.join_peakcall_tables <- function(peaks, blocks, pvalues, recalc = FALSE) {
  peaks$variant_ID <- peaks$peak_variant_ID
  if (recalc) {
    pvalues <- subset(pvalues, select = variant_ID:Pos)
    peaks_join <- dplyr::left_join(x = peaks, y = pvalues, by = "variant_ID")
    names(peaks_join)[-(1:3)] <- paste0("peak_", names(peaks_join)[-(1:3)])
    peaks_join <- subset(peaks_join, select = -variant_ID)
    blocks_join <- dplyr::left_join(x = blocks, y = pvalues, by = "variant_ID")
    hit <- match(peaks_join$peak_variant_ID, blocks_join$variant_ID)
    peaks_join$peak_negLog10P <- blocks_join$negLog10P[hit]
    dplyr::right_join(x = peaks_join, y = blocks_join, by = "peak_ID")
  } else {
    peaks_join <- dplyr::left_join(x = peaks, y = pvalues, by = "variant_ID")
    names(peaks_join)[-(1:3)] <- paste0("peak_", names(peaks_join)[-(1:3)])
    if ("peak_FDR" %in% names(peaks_join)) {
      peaks_join <- subset(peaks_join, select = -c(variant_ID, peak_FDR))
    } else {
      peaks_join <- subset(peaks_join, select = -variant_ID)
    }
    blocks_join <- dplyr::left_join(x = blocks, y = pvalues, by = "variant_ID")
    dplyr::right_join(x = peaks_join, y = blocks_join, by = "peak_ID")
  }
}


# Retrieve a peakcall result for a phenotype
#' @importFrom dplyr right_join left_join
#'
.get_peakcall <- function(object, pheno_name, recalc = FALSE){
  section <- if (recalc) "recalc" else "peakcall"

  if (.store_is_gds(object)) {
    return(.get_peakcall_gds(object, pheno_name, recalc))
  }

  peaks <- .store_read_peaks_df(object, section, pheno_name, "peaks")
  if (is.null(peaks) || nrow(peaks) == 0L) {
    return(NULL)
  }
  blocks <- .store_read_peaks_df(object, section, pheno_name, "blocks")
  meta <- .store_peak_meta(object, section, pheno_name)

  pvalues <- .get_scan(object = object, pheno_name = pheno_name, sel = c("FDR", "negLog10P"))
  out <- .join_peakcall_tables(peaks = peaks, blocks = blocks, pvalues = pvalues, recalc = recalc)

  attributes(out) <- c(
    attributes(out),
    list(signif = meta$signif, threshold = meta$threshold)
  )
  return(out)
}

.get_peakcall_gds <- function(object, pheno_name, recalc = FALSE){
  if (recalc) {
    peaks_gdsn <- "lazygas/recalc/peaks"
    blocks_gdsn <- "lazygas/recalc/blocks"
    peaks_col_names <- gdsfmt::get.attr.gdsn(
      node = gdsfmt::index.gdsn(node = object$root,
                                path = paste0(peaks_gdsn, "/", pheno_name))
    )
    blocks_col_names <- gdsfmt::get.attr.gdsn(
      node = gdsfmt::index.gdsn(node = object$root,
                                path = paste0(blocks_gdsn, "/", pheno_name))
    )
  } else {
    peaks_gdsn <- "lazygas/peakcall/peaks"
    blocks_gdsn <- "lazygas/peakcall/blocks"
    peaks_col_names <- gdsfmt::get.attr.gdsn(
      node = gdsfmt::index.gdsn(node = object$root, path = peaks_gdsn)
    )
    blocks_col_names <- gdsfmt::get.attr.gdsn(
      node = gdsfmt::index.gdsn(node = object$root, path = blocks_gdsn)
    )
  }

  peaks <- .get_data(object = object, node = paste0(peaks_gdsn, "/", pheno_name))
  if (length(peaks) == 0) {
    return(NULL)
  }
  if (!is.matrix(peaks)) {
    peaks <- data.frame(matrix(peaks, nrow = 1))
  } else if (recalc) {
    peaks <- data.frame(peaks)
  } else {
    peaks <- data.frame(t(peaks))
  }
  names(peaks) <- peaks_col_names$col_names

  blocks <- .get_data(object = object, node = paste0(blocks_gdsn, "/", pheno_name))
  if (!is.matrix(blocks)) {
    blocks <- data.frame(matrix(blocks, nrow = 1))
  } else if (recalc) {
    blocks <- data.frame(blocks)
  } else {
    blocks <- data.frame(t(blocks))
  }
  names(blocks) <- blocks_col_names$col_names

  peaks$variant_ID <- peaks$peak_variant_ID
  pvalues <- .get_scan(object = object, pheno_name = pheno_name, sel = c("FDR", "negLog10P"))
  out <- .join_peakcall_tables(peaks = peaks, blocks = blocks, pvalues = pvalues, recalc = recalc)
  attributes(out) <- c(
    attributes(out),
    list(signif = peaks_col_names$signif, threshold = peaks_col_names$threshold)
  )
  out
}


