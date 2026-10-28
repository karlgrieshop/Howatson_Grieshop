# Shared candidate definition for FigS3 and TableS2.
load_grieshop_candidate_ids <- function(path) {
  if (!file.exists(path)) {
    stop("Missing Grieshop et al. (2025) input: ", path, call. = FALSE)
  }
  de <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("FlyBaseID", "Sig") %in% names(de))) {
    stop("Grieshop DE input must contain FlyBaseID and Sig columns.", call. = FALSE)
  }
  ids <- unique(as.character(de$FlyBaseID[
    !is.na(de$Sig) & as.character(de$Sig) == "TRUE" &
      !is.na(de$FlyBaseID) & grepl("^FBgn", de$FlyBaseID)
  ]))
  if (length(ids) != 695L) {
    stop("Expected 695 significant Grieshop genes; found ", length(ids), ".", call. = FALSE)
  }
  ids
}
