#' @import methods rhdf5
#' 
#############################################################################

#' @export
setClass('h5.matrix', representation(did='H5IdComponent'))


#' @title S4 class for H5 matrices
#' @export
h5.matrix <- function(file, name)
{
  fid <- rhdf5::H5Fopen(file)
  did <- rhdf5::H5Dopen(fid, name)
  rhdf5::H5Fclose(fid)
  x <- new("h5.matrix", did=did)
  return(x)
}

#' @title Convert to base R matrix
#' @description Extract values from a \code{h5.matrix} object
#' and convert to a base R matrix object
#' @param x A h5.matrix object
#' @export
setMethod('as.matrix', signature(x='h5.matrix'),
          function(x) return(x[,]))


# #' @template as.h5.matrix_methods_template
# NULL
# 
# setMethod('as.h5.matrix', signature(x='H5IdComponent'),
#           function(x) {
#             if ( rhdf5::H5Iget_type(x) != 'H5I_DATASET') {
#               stop("x is not a H5 Dataset")
#             }
#             y <- h5.matrix(x)
#             return(y)
#           })
# 



#' @rdname h5.matrix
#' @export
setGeneric('is.h5.matrix', function(x) standardGeneric('is.h5.matrix'))

#' @rdname h5.matrix
setMethod('is.h5.matrix', signature(x='h5.matrix'),
          function(x) return(TRUE))

#' @rdname h5.matrix
setMethod('is.h5.matrix', definition=function(x) return(FALSE))

#' @title The Number of Rows/Columns of a h5.matrix
#' @description \code{nrow} and \code{ncol} return the number of
#' rows or columns present in a \code{h5.matrix} object.
#' @param x A h5.matrix object
#' @return An integer of length 1
#' @docType methods
#' @rdname ncol-methods
#' @export
setMethod('ncol', signature(x="h5.matrix"),
          function(x) return(rhdf5::H5Sget_simple_extent_dims(rhdf5::H5Dget_space(x@did))$size[2]))

#' @rdname ncol-methods
#' @export
setMethod('nrow', signature(x="h5.matrix"), 
          function(x) return(rhdf5::H5Sget_simple_extent_dims(rhdf5::H5Dget_space(x@did))$size[1]))

#' @title Dimensions of a h5.matrix object
#' @description Retrieve the dimensions of a \code{h5.matrix} object
#' @param x A \code{h5.matrix} object
#' @export
setMethod('dim', signature(x="h5.matrix"),
          function(x) return(rhdf5::H5Sget_simple_extent_dims(rhdf5::H5Dget_space(x@did))$size))

#' @title Length of a h5.matrix object
#' @description Get the length of a \code{h5.matrix} object
#' @param x A \code{h5.matrix} object
#' @export
setMethod('length', signature(x="h5.matrix"),
          function(x) return(prod(dim(x))))

GetElements.h5 <- function(x, i, j, drop=TRUE)
{
  if ( missing(i) )
    i <- seq_len(nrow(x))
  if ( missing(j) )
    j <- seq_len(ncol(x))
  
  if (!is.numeric(i) & !is.logical(i))
    stop("row indices must be numeric or logical vectors.")
  if (!is.numeric(j) & !is.logical(j))
    stop("column indices must be numeric or logical vectors..")
  if (is.logical(i)) {
    if (length(i) != nrow(x))
      stop("row vector length must match the number of rows of the matrix.")
    i <- which(i)
  }
  if (is.logical(j)) {
    if (length(j) != ncol(x))
      stop(paste("column vector length must match the number of",
                 "columns of the matrix."))
    j <- which(j)
  }
  
  sid <- rhdf5::H5Dget_space(x@did)
  sidmem <- rhdf5::H5Screate_simple(c(length(i), length(j)))
  rhdf5::H5Sselect_index(sid, list(i, j))
  
  mat = rhdf5::H5Dread(x@did, sid, sidmem)[, , drop=drop]

  return(mat)
}



GetAll.h5 <- function(x, drop=TRUE)
{
  return(rhdf5::H5Dread(x@did))
}

#' @title Extract or Replace h5.matrix elements
#' @name Extract,h5.matrix
#' @param x A \code{h5.matrix object}
#' @param i Indices specifying the rows
#' @param j Indices specifying the columns
#' @param drop Logical indication if reduce to minimum dimensions
#' @param value typically an array-like R object of similar class
#' @docType methods
#' @rdname extract-methods
#' @aliases [,h5.matrix,ANY,ANY,missing-method
#' @export
setMethod("[",
          signature(x = "h5.matrix", drop = "missing"),
          function(x, i, j, drop) return(GetElements.h5(x, i, j)))


#' @rdname extract-methods
#' @export
setMethod("[",
          signature(x = "h5.matrix", drop = "logical"),
          function(x, i, j, drop) return(GetElements.h5(x, i, j, drop)))


#' @rdname extract-methods
#' @export
setMethod("[",
          signature(x = "h5.matrix", i="missing", drop = "missing"),
          function(x, i, j, drop) return(GetElements.h5(x, i, j)))


#' @rdname extract-methods
#' @export
setMethod("[",
          signature(x = "h5.matrix", i="missing", drop = "logical"),
          function(x, i, j, drop) return(GetElements.h5(x, i, j)))


#' @rdname extract-methods
#' @export
setMethod("[",
          signature(x = "h5.matrix", j="missing", drop = "missing"),
          function(x, i, j, drop) return(GetElements.h5(x, i, j)))


#' @rdname extract-methods
#' @export
setMethod("[",
          signature(x = "h5.matrix", j="missing", drop = "logical"),
          function(x, i, j, drop) return(GetElements.h5(x, i, j)))


#' @rdname extract-methods
#' @export
setMethod("[",
          signature(x = "h5.matrix", i="missing", j="missing", drop = "missing"),
          function(x, i, j, drop) return(GetElements.h5(x, i, j)))


#' @rdname extract-methods
#' @export
setMethod("[",
          signature(x = "h5.matrix", i="missing", j="missing", drop = "logical"),
          function(x, i, j, drop) return(GetAll.h5(x, drop)))



SetElements.h5 <- function(x, i, j, value)
{
#   checkReadOnly(x)
  
  if ( missing(i) )
    i <- seq_len(nrow(x))
  if ( missing(j) )
    j <- seq_len(ncol(x))
  
  if (!is.numeric(i) & !is.logical(i))
    stop("row indices must be numeric or logical vectors.")
  if (!is.numeric(j) & !is.logical(j))
    stop("column indices must be numeric or logical vectors.")
  
  if (is.logical(i)) {
    if (length(i) != nrow(x))
      stop("row vector length must match the number of rows of the matrix.")
    i <- which(i)
  }
  if (is.logical(j)) {
    if (length(j) != ncol(x))
      stop(paste("column vector length must match the number of",
                 "columns of the matrix."))
    j <- which(j)
  }
  
  
  
  totalts <- max(as.double(length(i)), as.double(length(j)))
  # If we are assigning from a matrix, make sure the dimensions agree.
  if (is.matrix(value))
  {
    if (ncol(value) != length(j) | nrow(value) != length(i)) 
    {
      stop("Matrix dimensions do not agree with h5.matrix instance set size.")
    }
  }
  else if (length(value) != totalts) {
    # Otherwise, make sure we are assigning the correct number of things
    # (rep if necessary)
    numReps <- totalts / length(value)
    if (numReps != round(numReps)) 
    {
      stop(paste("number of items to replace is not a multiple of",
                 "replacement length"))
    }
    value <- rep(value, length.out=totalts)
  }
  
  
  sid <- rhdf5::H5Dget_space(x@did)
  sidmem <- rhdf5::H5Screate_simple(c(length(i), length(j)))
  rhdf5::H5Sselect_index(sid, list(i, j))
  rhdf5::H5Dwrite(x@did, value, sidmem, sid)
  return(x)
}


#' @rdname extract-methods
#' @export
setMethod('[<-',
          signature(x = "h5.matrix"),
          function(x, i, j, value) return(SetElements.h5(x, i, j, value)))

#' @rdname extract-methods
#' @export
setMethod('[<-',
          signature(x = "h5.matrix", i="missing"),
          function(x, i, j, value) return(SetElements.h5(x, i, j, value)))

#' @rdname extract-methods
#' @export
setMethod('[<-',
          signature(x = "h5.matrix", j="missing"),
          function(x, i, j, value) return(SetElements.h5(x, i, j, value)))

#' @rdname extract-methods
#' @export
setMethod('[<-',
          signature(x = "h5.matrix", i="missing", j="missing"),
          function(x, i, j, value) return(SetElements.h5(x, i, j, value)))


#' @title The Type of a h5.matrix Object
#' @description \code{typeof} returns the storage type of a 
#' \code{h5.matrix} object
#' @param x A \code{h5.matrix} object
#' @export
setMethod('typeof', signature(x="h5.matrix"),
          function(x) {
            return(rhdf5::H5Dget_type(x@did))
          }
)
