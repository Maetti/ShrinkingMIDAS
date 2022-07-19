#' Title
#'
#' @param sFolder
#'
#' @return
#' @export
#'
#' @examples
helper_createFolder <- function(sFolder) {
      if (!dir.exists(sFolder)) {
            logger::log_info("Create folder at: {sFolder}")
            dir.create(sFolder)
      }
}
