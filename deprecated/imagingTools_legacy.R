# DEPRECATED — these functions used a broken image-path pipeline.
# image_paths() hardcoded an author-specific path (TODO #6).
# imageval() read se$Image_paths[idx] which was always NULL (TODO #7).
# addImgPaths() is the old name for image_paths(); same issues.
# Use addThumbnailPaths() + plotThumbnails() instead.

#' Function that generates all the files paths for the images
#' @param parent_folder Path to where your various experimental data is stored
#' @param idx The well(s) you want to look at
#' @param plate_ID The plate ID(s) in question
#' @param location The location of the plate_ID in the folder names
#' @return A vector of image paths
#' @export
image_paths <- function(parent_folder, idx, plate_ID, location) {
  all_dirs <- list.dirs(parent_folder, full.names = TRUE, recursive = FALSE)
  matched_dir <- NULL
  short_idx <- gsub("^([A-Z])0*", "\\1", idx)

  for (d in all_dirs) {
    dir_name <- basename(d)
    parts <- strsplit(dir_name, "_")[[1]]
    if(length(parts) >= location && parts[location] == plate_ID) {
      particle_path <- file.path(d, "Particle_Analysis")
      if (!dir.exists(particle_path)) next

      subdirs <- list.dirs(particle_path, full.names = TRUE, recursive = FALSE)

      for (sub in subdirs) {
        sub_dir_name <- basename(sub)
        sub_parts <- strsplit(sub_dir_name, "_")[[1]]
        if (length(sub_parts) >= location && sub_parts[location] == plate_ID) {
          matched_dir <- sub
          break
        }
      }
      if (is.null(matched_dir)) break
    }
  }
  if (is.null(matched_dir)) {
    stop(paste("No folder with the plate_ID", plate_ID, "in the", location,"th position found"))
    return(NA)
  }
  img.list <- list.files(matched_dir, pattern = "\\.tif$", recursive = TRUE,
                         full.names = TRUE)
  pattern <- paste0(short_idx, "-\\d+")
  imgs <- img.list[grepl(pattern, img.list)]
  if (length(imgs) == 0) {
    warning(paste("No matching images found for well", idx, "in", matched_dir))
    return(NA)
  }
  return(imgs)
}

#' Helper function for imageval to brighten the images
#' @param img Image array
#' @param factor Factor by which the image should be brightened
#' @return A brightened image
#' @export
brighten_image <- function(img, factor = 1.5) {
  img_normalized <- (img - min(img)) / (max(img) - min(img))
  img_brightened <- img_normalized * factor
  img_brightened <- pmin(img_brightened, 1)
  return(img_brightened)
}

#' Plot the images with the different channels for a given well and plate
#' @param se SummarizedExperiment with Image_paths in colData
#' @param idx Well identifier
#' @param plate_ID Plate ID
#' @param green_slice Which slice is the green channel
#' @param red_slice Which slice is the red channel
#' @return A plot of four images
#' @export
imageval <- function(se, idx, plate_ID, green_slice = 2, red_slice = 3) {
  load_bright <- function(file, factor)
    brighten_image(readImage(file), factor = factor)

  short_idx <- gsub("^([A-Z])0*", "\\1", idx)

  make_grob <- function(img, src_slice = NULL, colour = c("red", "green", "blue")) {
    if (is.null(src_slice)) {
      if (inherits(img, "Image")) {
        x <- normalize(as.array(img))
      } else {
        x <- normalize(img)
      }
    } else {
      colour <- match.arg(colour)
      chan <- normalize(img[,,src_slice])
      rgb_arr <- array(0, dim = c(dim(chan), 3))
      rgb_arr[,, match(colour, c("red", "green", "blue"))] <- chan
      x <- rgb_arr
    }
    rasterGrob(x, interpolate = TRUE)
  }

  imgs <- se$Image_paths[idx]

  all_imgs <- imgs[[1]]
  bf_img <- all_imgs[grepl("BF\\.tif$", all_imgs)]
  img1 <- load_bright(bf_img, factor = 2)
  bf_channel <- normalize(img1)
  img1_grob <- make_grob(EBImage::rotate(img1, 180))

  fluorescent_img <- all_imgs[grepl("nm\\.tif$", all_imgs)]

  if (length(fluorescent_img) != 1) {
    stop(paste("Expected exactly 1 fluorescent image ending in 'nm.tif', but found", length(fluorescent_img)))
  }

  img2 <- load_bright(fluorescent_img, factor = 40)
  img2_grob_green <- make_grob(img2, src_slice = green_slice, colour = "green")
  img2_grob_red   <- make_grob(img2, src_slice = red_slice, colour = "red")

  comp_rgb <- array(0, dim = c(dim(img2[,,3]), 3))
  comp_rgb[,,1] <- normalize(img2[,,red_slice])
  comp_rgb[,,2] <- normalize(img2[,,green_slice])
  comp_rgb[,,3] <- bf_channel
  img2_grob_color <- rasterGrob(comp_rgb, interpolate = TRUE)

  grid.arrange(img1_grob, img2_grob_color, img2_grob_green, img2_grob_red, ncol = 2)
}
