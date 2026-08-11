#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(igraph)
  library(jsonlite)
})

cmd_args <- commandArgs(trailingOnly = TRUE)
source_network_id <- 128L
id_arg <- grep("^--source-network=", cmd_args, value = TRUE)
if (length(id_arg) == 1L) {
  source_network_id <- as.integer(sub("^--source-network=", "", id_arg))
}

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "reanalysis/scripts/build_kla_ddr_cytoscape_cluster_networks_33groups.R",
    mustWork = TRUE
  )
}
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)

analysis_name <- "kla_ddr_cytoscape_pathway_clusters_33groups_v1"
table_dir <- file.path(project_root, "reanalysis/results/tables", analysis_name)
figure_dir <- file.path(project_root, "reanalysis/results/figures", analysis_name)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

node_import_path <- file.path(table_dir, "cytoscape_node_import_table.csv")
if (!file.exists(node_import_path)) {
  stop(
    "Missing pathway cluster table. Run ",
    "cluster_kla_ddr_pathway_scores_33groups.R first."
  )
}

api_base <- "http://localhost:1234/v1"
identity_header <- add_headers(`Accept-Encoding` = "identity")

assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

get_json <- function(url) {
  response <- GET(url, identity_header)
  stop_for_status(response)
  fromJSON(
    rawToChar(content(response, as = "raw")),
    simplifyVector = FALSE
  )
}

post_json <- function(url, body) {
  response <- POST(
    url,
    identity_header,
    content_type_json(),
    body = toJSON(body, auto_unbox = TRUE, na = "null", null = "null")
  )
  stop_for_status(response)
  fromJSON(
    rawToChar(content(response, as = "raw")),
    simplifyVector = FALSE
  )
}

put_json <- function(url, body, query = list()) {
  response <- PUT(
    url,
    identity_header,
    content_type_json(),
    query = query,
    body = toJSON(body, auto_unbox = TRUE, na = "null", null = "null")
  )
  stop_for_status(response)
  invisible(TRUE)
}

run_command <- function(namespace, command, arguments = list()) {
  url <- paste0(
    api_base,
    "/commands/",
    URLencode(namespace, reserved = TRUE),
    "/",
    URLencode(command, reserved = TRUE)
  )
  response <- GET(url, identity_header, query = arguments)
  stop_for_status(response)
  invisible(rawToChar(content(response, as = "raw")))
}

download_view <- function(network_id, extension, output_path, width = NULL, height = NULL) {
  if (!is.null(width) && !is.null(height)) {
    stop("Cytoscape view export accepts width or height, not both.")
  }
  view_ids <- get_json(
    sprintf("%s/networks/%d/views", api_base, network_id)
  )
  if (length(view_ids) != 1L) {
    stop("Expected exactly one view for Cytoscape network ", network_id, ".")
  }
  view_id <- as.integer(view_ids[[1L]])
  url <- sprintf(
    "%s/networks/%d/views/%d.%s",
    api_base,
    network_id,
    view_id,
    extension
  )
  query <- list()
  if (!is.null(width)) query$w <- width
  if (!is.null(height)) query$h <- height
  response <- GET(url, identity_header, query = query)
  stop_for_status(response)
  writeBin(content(response, as = "raw"), output_path)
}

get_only_view_id <- function(network_id) {
  view_ids <- get_json(
    sprintf("%s/networks/%d/views", api_base, network_id)
  )
  if (length(view_ids) != 1L) {
    stop("Expected exactly one view for Cytoscape network ", network_id, ".")
  }
  as.integer(view_ids[[1L]])
}

get_view_coordinates <- function(network_id) {
  view_id <- get_only_view_id(network_id)
  view <- get_json(
    sprintf("%s/networks/%d/views/%d", api_base, network_id, view_id)
  )
  coordinates <- rbindlist(lapply(view$elements$nodes, function(node) {
    data.table(
      BaseAccession = node$data$BaseAccession,
      NodeSUID = as.integer(node$data$SUID),
      X = as.numeric(node$position$x),
      Y = as.numeric(node$position$y)
    )
  }))
  assert(
    nrow(coordinates) == 507L &&
      uniqueN(coordinates$BaseAccession) == 507L &&
      uniqueN(coordinates$NodeSUID) == 507L,
    paste("Invalid Cytoscape node view for network", network_id)
  )
  coordinates
}

set_view_coordinates <- function(network_id, coordinates) {
  view_id <- get_only_view_id(network_id)
  target <- get_view_coordinates(network_id)[, .(BaseAccession, NodeSUID)]
  values <- merge(
    target,
    coordinates[, .(BaseAccession, X, Y)],
    by = "BaseAccession",
    all.x = TRUE,
    sort = FALSE
  )
  assert(
    nrow(values) == 507L &&
      !anyNA(values[, .(X, Y)]) &&
      uniqueN(values[, .(X, Y)]) == 507L,
    paste("Cannot map all shared coordinates to network", network_id)
  )
  body <- lapply(seq_len(nrow(values)), function(i) {
    list(
      SUID = values$NodeSUID[[i]],
      view = list(
        list(
          visualProperty = "NODE_X_LOCATION",
          value = values$X[[i]]
        ),
        list(
          visualProperty = "NODE_Y_LOCATION",
          value = values$Y[[i]]
        )
      )
    )
  })
  put_json(
    sprintf(
      "%s/networks/%d/views/%d/nodes",
      api_base,
      network_id,
      view_id
    ),
    body,
    query = list(bypass = "false")
  )
}

compact_disconnected_components <- function(coordinates, edges) {
  graph <- graph_from_data_frame(
    edges[, .(Source, Target)],
    directed = FALSE,
    vertices = coordinates$BaseAccession
  )
  membership <- components(graph)$membership
  packed <- copy(coordinates)
  packed[, Component := unname(membership[BaseAccession])]
  component_sizes <- packed[, .N, by = Component][order(-N, Component)]
  main_component <- component_sizes$Component[[1L]]

  main <- packed[Component == main_component]
  main_center_x <- mean(range(main$X))
  main_center_y <- mean(range(main$Y))
  packed[Component == main_component, `:=`(
    X = X - main_center_x,
    Y = Y - main_center_y
  )]

  main <- packed[Component == main_component]
  node_diameter <- 28
  component_gap <- 34
  row_gap <- 34
  shelf_left <- min(main$X)
  shelf_right <- max(main$X)
  cursor_x <- shelf_left
  cursor_y <- max(main$Y) + 90
  row_height <- 0

  for (component_id in component_sizes[Component != main_component, Component]) {
    index <- which(packed$Component == component_id)
    local_x <- packed$X[index] - mean(range(packed$X[index]))
    local_y <- packed$Y[index] - mean(range(packed$Y[index]))

    if (length(index) == 1L) {
      local_x <- 0
      local_y <- 0
    } else {
      span <- max(diff(range(local_x)), diff(range(local_y)))
      if (!is.finite(span) || span == 0) {
        local_x <- seq(-25, 25, length.out = length(index))
        local_y <- 0
      } else {
        scale_factor <- 60 / span
        local_x <- local_x * scale_factor
        local_y <- local_y * scale_factor
      }
    }

    box_width <- max(diff(range(local_x)), node_diameter) + node_diameter
    box_height <- max(diff(range(local_y)), node_diameter) + node_diameter
    if (cursor_x + box_width > shelf_right && cursor_x > shelf_left) {
      cursor_x <- shelf_left
      cursor_y <- cursor_y + row_height + row_gap
      row_height <- 0
    }

    packed$X[index] <- cursor_x + box_width / 2 + local_x
    packed$Y[index] <- cursor_y + box_height / 2 + local_y
    cursor_x <- cursor_x + box_width + component_gap
    row_height <- max(row_height, box_height)
  }

  packed[, Component := NULL]
  assert(
    nrow(packed) == 507L &&
      uniqueN(packed$BaseAccession) == 507L &&
      uniqueN(packed[, .(X, Y)]) == 507L &&
      all(is.finite(packed$X)) &&
      all(is.finite(packed$Y)),
    "Compact shared STRING topology coordinates are invalid."
  )
  packed
}

safe_value <- function(x, name, fallback = NULL) {
  value <- x[[name]]
  if (is.null(value) || length(value) == 0L) fallback else value
}

api_version <- get_json(paste0(api_base, "/version"))
assert(
  identical(api_version$cytoscapeVersion, "3.10.4"),
  paste0(
    "This workflow was validated in Cytoscape 3.10.4; detected ",
    api_version$cytoscapeVersion,
    "."
  )
)

source_network <- get_json(
  sprintf("%s/networks/%d", api_base, source_network_id)
)
source_rows <- get_json(
  sprintf(
    "%s/networks/%d/tables/defaultnode/rows",
    api_base,
    source_network_id
  )
)
row_by_suid <- setNames(
  source_rows,
  vapply(source_rows, function(x) as.character(x$SUID), character(1))
)

node_import <- fread(node_import_path)
node_ids <- node_import$BaseAccession
assert(
  nrow(node_import) == 507L &&
    uniqueN(node_ids) == 507L,
  "Expected exactly 507 unique BaseAccession node rows."
)

source_node_data <- lapply(source_network$elements$nodes, `[[`, "data")
source_id_to_accession <- vapply(source_node_data, function(node) {
  row <- row_by_suid[[as.character(node$SUID)]]
  query_term <- safe_value(row, "query term", NA_character_)
  if (is.na(query_term) && identical(safe_value(row, "name"), "A8MVJ9")) {
    query_term <- "A8MVJ9"
  }
  query_term
}, character(1))
names(source_id_to_accession) <- vapply(
  source_node_data,
  function(node) as.character(node$id),
  character(1)
)

assert(
  length(source_id_to_accession) == 507L &&
    !anyNA(source_id_to_accession) &&
    setequal(source_id_to_accession, node_ids),
  "The live STRING network cannot be mapped one-to-one to the fixed 507 BaseAccession IDs."
)

source_node_rows <- rbindlist(lapply(source_node_data, function(node) {
  row <- row_by_suid[[as.character(node$SUID)]]
  accession <- source_id_to_accession[[as.character(node$id)]]
  mcl_value <- safe_value(row, "MCL_cluster_i3", NA_integer_)
  data.table(
    BaseAccession = accession,
    STRING_NodeSUID = as.integer(node$SUID),
    STRING_Identifier = safe_value(row, "name", accession),
    STRING_Mapped = accession != "A8MVJ9",
    STRING_DisplayName = safe_value(row, "display name", ""),
    MCL_cluster_i3 = as.integer(mcl_value)
  )
}))
source_node_rows[
  is.na(MCL_cluster_i3),
  MCL_cluster_label := "Unclustered"
]
source_node_rows[
  !is.na(MCL_cluster_i3),
  MCL_cluster_label := sprintf("MCL_%02d", MCL_cluster_i3)
]

source_edge_data <- lapply(source_network$elements$edges, `[[`, "data")
string_edges <- rbindlist(lapply(seq_along(source_edge_data), function(i) {
  edge <- source_edge_data[[i]]
  source_accession <- source_id_to_accession[[as.character(edge$source)]]
  target_accession <- source_id_to_accession[[as.character(edge$target)]]
  data.table(
    EdgeID = sprintf("STRING_%04d", i),
    Source = source_accession,
    Target = target_accession,
    Interaction = safe_value(edge, "interaction", "pp"),
    STRING_score = as.numeric(safe_value(edge, "stringdb_score", NA_real_)),
    STRING_experiments = as.numeric(
      safe_value(edge, "stringdb_experiments", NA_real_)
    ),
    STRING_databases = as.numeric(
      safe_value(edge, "stringdb_databases", NA_real_)
    ),
    STRING_coexpression = as.numeric(
      safe_value(edge, "stringdb_coexpression", NA_real_)
    ),
    STRING_textmining = as.numeric(
      safe_value(edge, "stringdb_textmining", NA_real_)
    )
  )
}))
assert(
  nrow(string_edges) == 4458L &&
    all(string_edges$Source %in% node_ids) &&
    all(string_edges$Target %in% node_ids) &&
    all(string_edges$STRING_score >= 0.70),
  "Expected 4,458 high-confidence STRING edges among the fixed proteins."
)

nodes <- merge(
  node_import,
  source_node_rows,
  by = "BaseAccession",
  all.x = TRUE,
  sort = FALSE
)
nodes[, DisplayLabel := fifelse(
  !is.na(GeneSymbol) & nzchar(GeneSymbol),
  GeneSymbol,
  BaseAccession
)]
nodes[, NodeTooltip := paste0(
  DisplayLabel,
  " | ",
  BaseAccession,
  "\n",
  fifelse(is.na(ProteinName), "", ProteinName)
)]
nodes[, K7_cluster_label := sprintf("K7_C%d", K7_Cluster)]
nodes[, K9_cluster_label := sprintf("K9_C%d", K9_Cluster)]
nodes[, K7_assignment_status := fifelse(K7_AllZero, "All-zero", "Assigned")]
nodes[, K9_assignment_status := fifelse(K9_AllZero, "All-zero", "Assigned")]
nodes[, STRING_mapping_status := fifelse(STRING_Mapped, "Mapped", "STRING-unmapped")]

assert(
  nrow(nodes) == 507L &&
    sum(nodes$STRING_Mapped) == 506L &&
    nodes[STRING_Mapped == FALSE, BaseAccession] == "A8MVJ9" &&
    sum(!is.na(nodes$MCL_cluster_i3)) == 446L &&
    uniqueN(na.omit(nodes$MCL_cluster_i3)) == 65L,
  "Live STRING/MCL node counts changed unexpectedly."
)

setorder(nodes, BaseAccession)
fwrite(nodes, file.path(table_dir, "cytoscape_string_nodes_507.csv"))
fwrite(string_edges, file.path(table_dir, "cytoscape_string_edges_4458.csv"))

mcl_summary <- nodes[
  ,
  .(ProteinCount = .N),
  by = .(MCL_cluster_label)
][order(MCL_cluster_label)]
fwrite(mcl_summary, file.path(table_dir, "string_mcl_i3_cluster_sizes.csv"))

mcl_palette <- setNames(
  hcl.colors(65L, palette = "Dynamic"),
  sprintf("MCL_%02d", 1:65)
)
mcl_palette <- c(mcl_palette, Unclustered = "#D9D9D9")
cluster_palette <- c(
  "#0072B2",
  "#D55E00",
  "#009E73",
  "#CC79A7",
  "#E69F00",
  "#56B4E9",
  "#F0E442",
  "#7F3C8D",
  "#999999"
)
k7_palette <- setNames(cluster_palette[1:7], sprintf("K7_C%d", 1:7))
k9_palette <- setNames(cluster_palette, sprintf("K9_C%d", 1:9))
color_key <- rbindlist(list(
  data.table(
    Version = "STRING_MCL_i3",
    Cluster = names(mcl_palette),
    Color = unname(mcl_palette)
  ),
  data.table(
    Version = "PathwayScore_k7_seed25",
    Cluster = names(k7_palette),
    Color = unname(k7_palette)
  ),
  data.table(
    Version = "PathwayScore_k9_seed25",
    Cluster = names(k9_palette),
    Color = unname(k9_palette)
  )
))
fwrite(color_key, file.path(table_dir, "cytoscape_cluster_color_key.csv"))

node_data_columns <- c(
  "BaseAccession",
  "DisplayLabel",
  "GeneSymbol",
  "ProteinName",
  "NodeTooltip",
  "STRING_Mapped",
  "STRING_mapping_status",
  "STRING_Identifier",
  "MCL_cluster_i3",
  "MCL_cluster_label",
  grep("^Score_", names(nodes), value = TRUE),
  "K7_AllZero",
  "K7_Cluster",
  "K7_ClusterLabel",
  "K7_cluster_label",
  "K7_assignment_status",
  "K9_AllZero",
  "K9_Cluster",
  "K9_ClusterLabel",
  "K9_cluster_label",
  "K9_assignment_status"
)

make_network_json <- function(network_name) {
  node_elements <- lapply(seq_len(nrow(nodes)), function(i) {
    values <- as.list(nodes[i, ..node_data_columns])
    values$id <- nodes$BaseAccession[[i]]
    values$name <- nodes$BaseAccession[[i]]
    list(data = values)
  })
  edge_elements <- lapply(seq_len(nrow(string_edges)), function(i) {
    list(data = list(
      id = string_edges$EdgeID[[i]],
      source = string_edges$Source[[i]],
      target = string_edges$Target[[i]],
      interaction = string_edges$Interaction[[i]],
      STRING_score = string_edges$STRING_score[[i]],
      STRING_experiments = string_edges$STRING_experiments[[i]],
      STRING_databases = string_edges$STRING_databases[[i]],
      STRING_coexpression = string_edges$STRING_coexpression[[i]],
      STRING_textmining = string_edges$STRING_textmining[[i]]
    ))
  })
  list(
    data = list(
      name = network_name,
      Scope = "33 paired groups; 507 Kla-intersection-DDR proteins",
      IdentifierKey = "isoform-stripped UniProt BaseAccession",
      STRINGConfidenceCutoff = 0.70,
      STRINGAdditionalInteractors = 0L,
      STRINGMappedProteins = 506L,
      STRINGUnmappedProtein = "A8MVJ9"
    ),
    elements = list(nodes = node_elements, edges = edge_elements)
  )
}

create_network <- function(network_name) {
  existing_ids <- unlist(get_json(paste0(api_base, "/networks")), use.names = FALSE)
  existing_names <- vapply(existing_ids, function(network_id) {
    rows <- get_json(
      sprintf(
        "%s/networks/%d/tables/defaultnetwork/rows",
        api_base,
        network_id
      )
    )
    safe_value(rows[[1L]], "name", "")
  }, character(1))
  matching_ids <- as.integer(existing_ids[existing_names == network_name])
  if (length(matching_ids) == 1L) {
    return(matching_ids[[1L]])
  }
  if (length(matching_ids) > 1L) {
    stop(
      "More than one live Cytoscape network is named '",
      network_name,
      "'. Resolve the duplicate names before rerunning."
    )
  }
  response <- post_json(
    paste0(
      api_base,
      "/networks?format=json&collection=",
      URLencode("Kla_DDR_507_cluster_comparison", reserved = TRUE),
      "&title=",
      URLencode(network_name, reserved = TRUE)
    ),
    make_network_json(network_name)
  )
  as.integer(response$networkSUID)
}

discrete_mapping <- function(column, column_type, visual_property, values) {
  list(
    mappingType = "discrete",
    mappingColumn = column,
    mappingColumnType = column_type,
    visualProperty = visual_property,
    map = lapply(names(values), function(key) {
      list(key = key, value = unname(values[[key]]))
    })
  )
}

make_style <- function(style_name, cluster_column, palette, status_column) {
  list(
    title = style_name,
    defaults = list(
      list(visualProperty = "NETWORK_BACKGROUND_PAINT", value = "#FFFFFF"),
      list(visualProperty = "NODE_SIZE", value = 28),
      list(visualProperty = "NODE_WIDTH", value = 28),
      list(visualProperty = "NODE_HEIGHT", value = 28),
      list(visualProperty = "NODE_SHAPE", value = "ELLIPSE"),
      list(visualProperty = "NODE_FILL_COLOR", value = "#D9D9D9"),
      list(visualProperty = "NODE_BORDER_PAINT", value = "#FFFFFF"),
      list(visualProperty = "NODE_BORDER_WIDTH", value = 0.8),
      list(visualProperty = "NODE_LABEL_FONT_SIZE", value = 8),
      list(visualProperty = "NODE_LABEL_TRANSPARENCY", value = 0),
      list(visualProperty = "EDGE_STROKE_UNSELECTED_PAINT", value = "#73808C"),
      list(visualProperty = "EDGE_TRANSPARENCY", value = 8),
      list(visualProperty = "EDGE_WIDTH", value = 0.2),
      list(visualProperty = "EDGE_CURVED", value = FALSE)
    ),
    mappings = list(
      discrete_mapping(
        cluster_column,
        "String",
        "NODE_FILL_COLOR",
        palette
      ),
      discrete_mapping(
        status_column,
        "String",
        "NODE_SHAPE",
        c(Assigned = "ELLIPSE", `All-zero` = "DIAMOND")
      ),
      list(
        mappingType = "passthrough",
        mappingColumn = "DisplayLabel",
        mappingColumnType = "String",
        visualProperty = "NODE_LABEL"
      ),
      list(
        mappingType = "passthrough",
        mappingColumn = "NodeTooltip",
        mappingColumnType = "String",
        visualProperty = "NODE_TOOLTIP"
      ),
      list(
        mappingType = "continuous",
        mappingColumn = "STRING_score",
        mappingColumnType = "Double",
        visualProperty = "EDGE_WIDTH",
        points = list(
          list(value = 0.70, lesser = "0.10", equal = "0.10", greater = "0.10"),
          list(value = 1.00, lesser = "0.70", equal = "0.70", greater = "0.70")
        )
      ),
      list(
        mappingType = "continuous",
        mappingColumn = "STRING_score",
        mappingColumnType = "Double",
        visualProperty = "EDGE_TRANSPARENCY",
        points = list(
          list(value = 0.70, lesser = "3", equal = "3", greater = "3"),
          list(value = 1.00, lesser = "35", equal = "35", greater = "35")
        )
      )
    )
  )
}

make_mcl_style <- function() {
  style <- make_style(
    "Kla_DDR_STRING_MCL_i3_overview_v4",
    "MCL_cluster_label",
    mcl_palette,
    "STRING_mapping_status"
  )
  style$mappings[[2L]] <- discrete_mapping(
    "STRING_mapping_status",
    "String",
    "NODE_SHAPE",
    c(Mapped = "ELLIPSE", `STRING-unmapped` = "DIAMOND")
  )
  style
}

create_style <- function(style) {
  existing_styles <- unlist(get_json(paste0(api_base, "/styles")), use.names = FALSE)
  if (style$title %in% existing_styles) {
    return(style$title)
  }

  style_file <- tempfile(fileext = ".json")
  response_file <- tempfile(fileext = ".json")
  error_file <- tempfile(fileext = ".log")
  on.exit(unlink(c(style_file, response_file, error_file)), add = TRUE)
  writeLines(
    toJSON(style, auto_unbox = TRUE, na = "null", null = "null"),
    style_file,
    useBytes = TRUE
  )
  response_text <- system2(
    "curl",
    args = c(
      "-sS",
      "--fail-with-body",
      "-X",
      "POST",
      "-H",
      "Content-Type: application/json",
      "--data-binary",
      paste0("@", style_file),
      paste0(api_base, "/styles")
    ),
    stdout = response_file,
    stderr = error_file
  )
  status <- as.integer(response_text)
  if (is.na(status)) status <- 0L
  if (status != 0L) {
    styles_after_request <- unlist(
      get_json(paste0(api_base, "/styles")),
      use.names = FALSE
    )
    if (style$title %in% styles_after_request) {
      return(style$title)
    }
    error_text <- if (file.exists(error_file)) {
      readLines(error_file, warn = FALSE)
    } else {
      character()
    }
    response_body <- if (file.exists(response_file)) {
      readLines(response_file, warn = FALSE)
    } else {
      character()
    }
    stop(
      "Cytoscape style creation failed: ",
      paste(c(error_text, response_body), collapse = "\n")
    )
  }
  response <- fromJSON(
    paste(readLines(response_file, warn = FALSE), collapse = "\n")
  )
  response$title
}

apply_style <- function(style_name, network_id) {
  response <- GET(
    paste0(
      api_base,
      "/apply/styles/",
      URLencode(style_name, reserved = TRUE),
      "/",
      network_id
    ),
    identity_header
  )
  stop_for_status(response)
  invisible(TRUE)
}

network_specs <- list(
  list(
    key = "string_mcl_i3_shared_string_topology",
    name = "Kla_DDR_507_STRING_MCL_i3_shared_STRING_topology",
    style = make_mcl_style(),
    layout = "shared_compact_STRING_topology"
  ),
  list(
    key = "pathway_score_k7_seed25_shared_string_topology",
    name = "Kla_DDR_507_pathway_score_k7_seed25_shared_STRING_topology",
    style = make_style(
      "Kla_DDR_pathway_score_k7_seed25_overview_v4",
      "K7_cluster_label",
      k7_palette,
      "K7_assignment_status"
    ),
    layout = "shared_compact_STRING_topology"
  ),
  list(
    key = "pathway_score_k9_seed25_shared_string_topology",
    name = "Kla_DDR_507_pathway_score_k9_seed25_shared_STRING_topology",
    style = make_style(
      "Kla_DDR_pathway_score_k9_seed25_overview_v4",
      "K9_cluster_label",
      k9_palette,
      "K9_assignment_status"
    ),
    layout = "shared_compact_STRING_topology"
  )
)

created_networks <- lapply(network_specs, function(spec) {
  network_id <- create_network(spec$name)
  style_name <- create_style(spec$style)
  apply_style(style_name, network_id)
  list(spec = spec, network_id = network_id, style_name = style_name)
})

topology_network_id <- created_networks[[1L]]$network_id
topology_network_name <- created_networks[[1L]]$spec$name
run_command(
  "network",
  "set current",
  list(network = topology_network_name)
)
run_command(
  "layout",
  "force-directed",
  list(
    network = topology_network_name,
    defaultSpringLength = 120,
    edgeAttribute = "STRING_score",
    isDeterministic = "true",
    numIterations = 1200,
    singlePartition = "true"
  )
)
raw_topology_coordinates <- get_view_coordinates(topology_network_id)
assert(
  uniqueN(raw_topology_coordinates[, .(X, Y)]) > 400L,
  "Cytoscape did not apply the STRING force-directed layout to the target network."
)
shared_coordinates <- compact_disconnected_components(
  raw_topology_coordinates,
  string_edges
)

for (entry in created_networks) {
  set_view_coordinates(entry$network_id, shared_coordinates)
  run_command(
    "network",
    "set current",
    list(network = entry$spec$name)
  )
  run_command(
    "view",
    "fit content",
    list(view = entry$spec$name)
  )
}

network_manifest <- rbindlist(lapply(created_networks, function(entry) {
  spec <- entry$spec
  network_id <- entry$network_id
  style_name <- entry$style_name

  png_path <- file.path(
    figure_dir,
    paste0("kla_ddr_507_", spec$key, ".png")
  )
  svg_path <- file.path(
    figure_dir,
    paste0("kla_ddr_507_", spec$key, ".svg")
  )
  pdf_path <- file.path(
    figure_dir,
    paste0("kla_ddr_507_", spec$key, ".pdf")
  )
  download_view(network_id, "png", png_path, width = 8000)
  download_view(network_id, "svg", svg_path)
  download_view(network_id, "pdf", pdf_path)

  coordinates <- get_view_coordinates(network_id)[
    ,
    .(BaseAccession, X, Y)
  ]
  fwrite(
    coordinates,
    file.path(table_dir, paste0("cytoscape_layout_", spec$key, ".csv"))
  )

  data.table(
    Version = spec$key,
    NetworkName = spec$name,
    NetworkSUID = network_id,
    StyleName = style_name,
    Layout = spec$layout,
    NodeCount = nrow(nodes),
    EdgeCount = nrow(string_edges),
    PNG = basename(png_path),
    SVG = basename(svg_path),
    PDF = basename(pdf_path)
  )
}))
fwrite(
  network_manifest,
  file.path(table_dir, "cytoscape_network_manifest.csv")
)

session_path <- file.path(
  table_dir,
  "kla_ddr_507_cluster_comparison_shared_string_topology_33groups_v2.cys"
)
run_command("session", "save as", list(file = session_path))

cat(
  "Created three Cytoscape comparison networks and saved the session:\n",
  session_path,
  "\n",
  paste(
    paste0(network_manifest$Version, " (SUID ", network_manifest$NetworkSUID, ")"),
    collapse = "\n"
  ),
  "\n",
  sep = ""
)
