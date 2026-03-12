#' @importFrom magrittr %>%
NULL

# ============================================================
#  localImgBrowser — client-side folder picker for Shiny
# ============================================================
#
#  Works locally AND when deployed on shinyapps.io / Shiny Server
#  because images are read entirely in the browser.
#  No files are uploaded to the server.
#
#  Filename convention (Cluster Analysis thumbnails):
#    ..._<Well>-<site>_<Channel>_<Class>_crop.jpg
#  e.g. Exp_18T39265_A01-1_BF_class1_crop.jpg
#
#  The module extracts the first-level sub-directory name from
#  webkitRelativePath so multi-plate folders work (one sub-dir per
#  plate).  Keys in the returned map:
#    "subdir|well.channel.class"  ->  blob:URL
#  where subdir is "" when files sit directly in the selected root.

# ------------------------------------------------------------------
#  UI
# ------------------------------------------------------------------

#' Client-side local image folder browser
#'
#' Drop-in Shiny module that lets users select a folder on their own
#' computer, even when the app is deployed on a remote server.
#' Files are read in the browser; only filename metadata and blob: URLs
#' reach the Shiny server.  Blob: URLs are valid in the same browser
#' session, so \code{<img src=blob:...>} tags rendered by
#' \code{renderUI} display correctly whether local or deployed.
#'
#' @param id Module namespace id.
#' @param label Button label shown to the user.
#' @export
localImgBrowserUI <- function(id, label = "Select image folder\u2026") {

  ns  <- NS(id)
  jsn <- gsub("[^A-Za-z0-9]", "_", ns("store"))   # safe JS global name

  tagList(
    # Hidden file input — webkitdirectory set via JS after DOM ready
    tags$input(
      id    = ns("folder"),
      type  = "file",
      style = "display:none; position:absolute;"
    ),

    shiny::actionButton(
      ns("browse"), label,
      icon  = shiny::icon("folder-open"),
      style = "width:100%; margin-bottom:4px;"
    ),
    tags$div(
      id    = ns("status"),
      style = "font-size:11px; color:#888; min-height:14px; word-break:break-word;"
    ),

    tags$script(shiny::HTML(paste0("
(function(){
  function byId(x){ return document.getElementById(x); }
  function pad2(n){ return n < 10 ? '0' + n : '' + n; }

  /* Parse thumbnail filename -> {well, channel, cls} or null */
  var RE = /_([A-P][0-9]{1,2}-[0-9]+)_([^_]+)_([^_]+)_crop\\.jpe?g$/i;
  function parseThumb(name){
    var m = name.match(RE);
    if(!m) return null;
    var raw = m[1].split('-')[0];
    var ltr = raw.charAt(0).toUpperCase();
    var num = parseInt(raw.slice(1), 10);
    if(isNaN(num) || num < 1 || num > 24) return null;
    return { well: ltr + pad2(num), channel: m[2], cls: m[3] };
  }

  /* Add webkitdirectory once element exists */
  function initInput(){
    var inp = byId('", ns("folder"), "');
    if(!inp){ setTimeout(initInput, 60); return; }
    inp.setAttribute('webkitdirectory', '');
    inp.setAttribute('directory', '');
  }
  initInput();

  /* Button triggers hidden input */
  $(document).on('click', '#", ns("browse"), "', function(e){
    e.preventDefault();
    var inp = byId('", ns("folder"), "');
    if(inp) inp.click();
  });

  /* Process selected files */
  $(document).on('change', '#", ns("folder"), "', function(e){
    var all   = Array.from(e.target.files);
    var files = all.filter(function(f){ return /\\.(jpg|jpeg)$/i.test(f.name); });
    var status = byId('", ns("status"), "');

    if(files.length === 0){
      status.textContent = 'No JPEG files found.';
      Shiny.setInputValue('", ns("imgdata"), "', null, {priority:'event'});
      return;
    }
    status.textContent = 'Processing ' + files.length + ' images\u2026';

    /* Revoke previous blob: URLs to free browser memory */
    if(window['", jsn, "']){
      Object.values(window['", jsn, "']).forEach(function(u){
        try{ URL.revokeObjectURL(u); } catch(x){}
      });
    }

    /* filename -> blob URL */
    var store = {};
    files.forEach(function(f){ store[f.name] = URL.createObjectURL(f); });
    window['", jsn, "'] = store;

    /* Build key map and collect unique subdirs / channels / classes
       Key: 'subdir|well.channel.class' -> blob URL
       subdir = first-level sub-dir inside the selected root folder,
                '' if files sit directly in root                       */
    var wellMap  = {};
    var subdirSet = {}, chanSet = {}, clsSet = {};

    files.forEach(function(f){
      var p = parseThumb(f.name);
      if(!p) return;

      var relParts = (f.webkitRelativePath || f.name).split('/');
      /* relParts[0] = root folder name (always present)
         relParts[1] = first subdir OR filename (if flat layout)  */
      var subdir = (relParts.length > 2) ? relParts[1] : '';

      var key = subdir + '|' + p.well + '.' + p.channel + '.' + p.cls;
      if(!wellMap[key]) wellMap[key] = store[f.name];  /* first site wins */

      subdirSet[subdir]  = 1;
      chanSet[p.channel] = 1;
      clsSet[p.cls]      = 1;
    });

    var subdirs  = Object.keys(subdirSet).sort();
    var channels = Object.keys(chanSet).sort();
    var classes  = Object.keys(clsSet).sort();

    var nKeys    = Object.keys(wellMap).length;
    var plateTxt = (subdirs.length > 0 && subdirs[0] !== '')
                 ? (' \u2014 ' + subdirs.length + ' plate folder(s)') : '';
    status.textContent = files.length + ' JPEGs' + plateTxt
                       + ' \u2014 ' + nKeys + ' well\u00d7ch\u00d7class entries';

    Shiny.setInputValue(
      '", ns("imgdata"), "',
      { n: files.length, map: wellMap,
        subdirs: subdirs, channels: channels, classes: classes },
      {priority: 'event'}
    );
  });
})();
    ")))
  )
}


# ------------------------------------------------------------------
#  Server
# ------------------------------------------------------------------

#' @rdname localImgBrowserUI
#' @return A reactive returning \code{list(n, map, subdirs, channels, classes)}.
#'   \code{map} is a named character vector with keys
#'   \code{"subdir|well.channel.class"} -> blob: URL.
#'   \code{subdirs} / \code{channels} / \code{classes} are sorted
#'   character vectors of unique values for populating dropdowns.
#' @export
localImgBrowserServer <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::reactive({
      shiny::req(input$imgdata)
      list(
        n        = input$imgdata$n,
        map      = unlist(input$imgdata$map),
        subdirs  = unlist(input$imgdata$subdirs),
        channels = unlist(input$imgdata$channels),
        classes  = unlist(input$imgdata$classes)
      )
    })
  })
}


# ------------------------------------------------------------------
#  filterImgMap — used by tinySEV server to build well -> URL lookup
# ------------------------------------------------------------------

#' Filter the client-side image map to a well -> URL named vector
#'
#' Extracts the \code{well -> url} vector for a specific plate
#' sub-directory, channel, and class from the map returned by
#' \code{localImgBrowserServer}.
#'
#' @param map Named character vector from \code{localImgBrowserServer()$map}.
#' @param subdir First-level sub-directory name (plate folder), or \code{""}
#'   if images sit in the root of the selected folder.
#' @param channel Channel string, e.g. \code{"BF"}.  \code{NULL} = no filter.
#' @param img_class Class string, e.g. \code{"class1"}.  \code{NULL} = no filter.
#' @return Named character vector \code{well -> url}.
#' @export
filterImgMap <- function(map, subdir = "", channel = NULL, img_class = NULL) {
  if (is.null(map) || length(map) == 0)
    return(setNames(character(0), character(0)))

  keys     <- names(map)
  split1   <- strsplit(keys, "|", fixed = TRUE)
  subdirs_k <- vapply(split1, function(x) x[1], character(1))
  wcc_k     <- vapply(split1, function(x) x[2], character(1))
  split2    <- strsplit(wcc_k, ".", fixed = TRUE)
  wells_k   <- vapply(split2, function(x) x[1], character(1))
  chans_k   <- vapply(split2, function(x) x[2], character(1))
  cls_k     <- vapply(split2, function(x) x[3], character(1))

  mask <- subdirs_k == subdir
  if (!is.null(channel))   mask <- mask & (chans_k == channel)
  if (!is.null(img_class)) mask <- mask & (cls_k   == img_class)

  out <- map[mask]
  names(out) <- wells_k[mask]
  out[!duplicated(names(out))]   # keep first entry per well
}


# ------------------------------------------------------------------
#  plateGridUI — renders the 384-well image grid inside renderUI
# ------------------------------------------------------------------

#' Render a 384-well plate grid from a well -> URL map
#'
#' Designed for use inside \code{renderUI}.  URLs may be blob: (client-side)
#' or relative served paths via \code{addResourcePath} (server-side) --
#' the rendering code is identical for both.
#'
#' @param url_map Named character vector \code{well -> URL} as returned by
#'   \code{filterImgMap()} or built from server-side \code{addResourcePath}
#'   paths.
#' @param coldata_colors Optional named vector \code{well -> hex_color} for a
#'   semi-transparent colData overlay on each tile.
#' @param click_input_id Full (namespaced) Shiny input id that receives the
#'   well id on click.  \code{NULL} to disable.
#' @param hover_input_id Full (namespaced) Shiny input id that receives the
#'   well id on mouse-over.  \code{NULL} to disable.
#' @param well_size_px Tile size in pixels (default 36).
#' @return A \code{shiny::tagList} ready for use in \code{renderUI}.
#' @export
plateGridUI <- function(url_map,
                        coldata_colors = NULL,
                        click_input_id = NULL,
                        hover_input_id = NULL,
                        well_size_px   = 36) {

  rows      <- LETTERS[1:16]
  cols      <- sprintf("%02d", 1:24)
  all_wells <- as.vector(t(outer(rows, cols, paste0)))
  sz        <- paste0(well_size_px, "px")

  well_tags <- lapply(all_wells, function(w) {
    url     <- url_map[w]
    has_img <- !is.null(url) && !is.na(url) && nzchar(url)

    overlay <- if (!is.null(coldata_colors)) {
      col <- coldata_colors[w]
      if (!is.null(col) && !is.na(col))
        tags$div(style = paste0(
          "position:absolute;top:0;left:0;",
          "width:100%;height:100%;",
          "background:", col, ";pointer-events:none;"
        ))
    }

    onclick     <- if (!is.null(click_input_id))
      sprintf("Shiny.setInputValue('%s','%s',{priority:'event'})",
              click_input_id, w)
    onmouseover <- if (!is.null(hover_input_id))
      sprintf("Shiny.setInputValue('%s','%s',{priority:'event'})",
              hover_input_id, w)

    tags$div(
      class       = if (has_img) "imgplate-well" else "imgplate-well imgplate-empty",
      style       = paste0("width:", sz, ";height:", sz, ";"),
      onclick     = onclick,
      onmouseover = onmouseover,
      `data-well` = w,
      if (has_img) tags$img(src = url, loading = "lazy"),
      overlay,
      tags$span(class = "imgplate-label", w)
    )
  })

  tagList(
    tags$style(HTML(paste0("
.imgplate-grid{
  display:grid;
  grid-template-columns:repeat(24,", sz, ");
  grid-auto-rows:", sz, ";
  gap:2px;width:max-content;
}
.imgplate-well{
  position:relative;overflow:hidden;border-radius:5px;
  border:1px solid #333;background:#2a2a2a;cursor:pointer;
}
.imgplate-well img{width:100%;height:100%;object-fit:cover;display:block;}
.imgplate-empty{background:#1e1e1e;}
.imgplate-label{
  position:absolute;bottom:1px;right:2px;
  background:rgba(0,0,0,.6);color:#fff;
  font-size:7px;padding:0 2px;border-radius:3px;pointer-events:none;
}
    "))),
    tags$div(class = "imgplate-grid", well_tags)
  )
}
