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


#' Title
#'
#' @param i
#'
#' @return
#' @export
#'
#' @examples
helper_create_number_name <- function(i) {

      ## get simulated data
      if (i < 10) {
            nSimRunSave <- paste0("00", i)
      } else if (i >= 10 & i < 100) {
            nSimRunSave <- paste0("0", i)
      } else {
            nSimRunSave <- as.character(i)
      }
      nSimRunSave
}
