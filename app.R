# ── Multiplex Blockmodeling Explorer ──────────────────────────────────────────
# Interactive Shiny app: explore how weighting M1/M2/M3 shifts the optimal
# task grouping produced by generalised blockmodeling.
#
# Binary blockmodeling only (APPROACH = "bin").
# All helper functions are copied verbatim from sim_unimode.qmd unless noted.
#
# To run: open this file in RStudio and click "Run App", or call
#   shiny::runApp("path/to/blockmodel_app")
# ──────────────────────────────────────────────────────────────────────────────

library(shiny)
library(blockmodeling)

# ── Global constants ───────────────────────────────────────────────────────────
# Reduced reps vs. the thesis run for interactive speed.
# The simulated data is well-structured; 10–20 restarts reliably find the optimum.
APPROACH   <- "bin"
SEED_AOP   <- 101L    # seed for finding all optimal partitions at k*
SEED_KLOOP <- 2025L   # seed for the k-selection loop
REPS_AOP   <- 10L     # restarts in the AOP search at k*
REPS_KLOOP <- 20L     # restarts per k in the k-selection loop

# ── sort_partition (verbatim from sim_unimode.qmd) ────────────────────────────
sort_partition <- function(clu0, node_names) {
  clusters <- sort(unique(clu0))

  split_name <- function(nm) {
    i <- regexpr("-", nm)[[1]]
    if (i > 0) c(prefix = substr(nm, 1, i - 1),
                 suffix = substr(nm, i + 1, nchar(nm)))
    else        c(prefix = nm, suffix = "")
  }
  parsed   <- do.call(rbind, lapply(node_names, split_name))
  prefixes <- parsed[, "prefix"]
  suffixes <- parsed[, "suffix"]
  has_suffixes <- any(nchar(suffixes) > 0)

  clus_prefix <- lapply(clusters, function(g) unique(prefixes[clu0 == g]))
  clus_suffix <- lapply(clusters, function(g) unique(suffixes[clu0 == g]))
  clus_nodes  <- lapply(clusters, function(g) sort(node_names[clu0 == g]))

  is_prefix_pure <- lengths(clus_prefix) == 1
  is_suffix_pure <- has_suffixes & (lengths(clus_suffix) == 1)

  all_pfx      <- sort(unique(prefixes))
  sorted_nodes <- sort(node_names)

  secondary_key <- numeric(length(clusters))
  for (i in seq_along(clusters)) {
    pp  <- is_prefix_pure[i];  sp  <- is_suffix_pure[i]
    pre <- clus_prefix[[i]];   suf <- clus_suffix[[i]]
    fn  <- clus_nodes[[i]][1L]
    if (pp) {
      secondary_key[i] <- match(pre, all_pfx) * 10000 +
                          min(suppressWarnings(as.numeric(suf)), na.rm = TRUE)
    } else if (sp) {
      secondary_key[i] <- 1e6 + suppressWarnings(as.numeric(suf))
    } else {
      secondary_key[i] <- 2e6 + match(fn, sorted_nodes)
    }
  }

  cluster_order <- clusters[order(secondary_key)]
  node_order <- unlist(lapply(cluster_order, function(g) {
    idx <- which(clu0 == g)
    idx[order(node_names[idx])]
  }))
  sort_method <- if (all(is_prefix_pure)) "product"
                 else if (all(is_suffix_pure)) "step"
                 else "mixed"

  list(cluster_order = cluster_order,
       node_order    = node_order,
       sort_method   = sort_method)
}

# ── prespec_block / make_blocks (verbatim) ────────────────────────────────────
prespec_block <- function(prespec, on_diag) {
  switch(prespec,
    comnul = if (on_diag) "com" else "nul",
    dncnul = if (on_diag) "dnc" else "nul",
    comdnc = if (on_diag) "com" else "dnc",
    nuldnc = if (on_diag) "nul" else "dnc",
    stop(paste0('Unknown prespec: "', prespec, '"'))
  )
}

make_blocks <- function(k, prespec = "comnul") {
  arr <- array(NA_character_, dim = c(1L, 1L, k, k))
  for (ri in seq_len(k))
    for (ci in seq_len(k))
      arr[1, 1, ri, ci] <- prespec_block(prespec, ri == ci)
  arr
}

# ── per_layer_err (verbatim, binary branch only) ───────────────────────────────
per_layer_err <- function(M, clu, prespec, approach = "bin") {
  R  <- dim(M)[3]
  k  <- max(clu)
  er <- numeric(R)
  for (r in seq_len(R)) {
    mat <- M[,,r]
    for (ri in seq_len(k)) {
      for (ci in seq_len(k)) {
        ri_idx <- which(clu == ri)
        ci_idx <- which(clu == ci)
        b      <- mat[ri_idx, ci_idx, drop = FALSE]
        btype  <- prespec_block(prespec[r], ri == ci)
        if (btype == "com") {
          if (ri == ci) diag(b) <- 1L   # ignore self-loops
          er[r] <- er[r] + sum(b == 0)
        } else if (btype == "nul") {
          er[r] <- er[r] + sum(b == 1)
        }
        # dnc: contributes 0
      }
    }
  }
  er
}

# ── generate_task_matrices ────────────────────────────────────────────────────
# Binary blockmodeling: all same-step pairs are tied in M1, all same-product
# pairs are tied in M2. Tie intensity variation (0.5 vs 1) from the valued
# blockmodeling version has been dropped.
generate_task_matrices <- function(N_products, N_steps) {
  if (N_products > 26) stop("N_products must be ≤ 26.")

  products  <- LETTERS[seq_len(N_products)]
  steps     <- as.character(seq_len(N_steps))
  last_step <- as.character(N_steps)

  nodes        <- unlist(lapply(products, function(p) paste0(p, "-", steps)))
  n            <- length(nodes)
  node_product <- sub("-.*", "", nodes)
  node_step    <- sub(".*-", "", nodes)

  M1 <- matrix(0, n, n, dimnames = list(nodes, nodes))
  M2 <- matrix(0, n, n, dimnames = list(nodes, nodes))
  M3 <- matrix(0, n, n, dimnames = list(nodes, nodes))

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      pi <- node_product[i];  si <- node_step[i]
      pj <- node_product[j];  sj <- node_step[j]

      # M1: same step, different product → always 1
      if (si == sj && pi != pj)
        M1[i, j] <- 1

      # M2: same product, different step → always 1
      if (pi == pj && si != sj)
        M2[i, j] <- 1

      # M3: goal conflict — last step vs. all other steps, same product
      if (pi == pj && xor(si == last_step, sj == last_step))
        M3[i, j] <- 1
    }
  }
  list(M1 = M1, M2 = M2, M3 = M3)
}

# ── run_blockmodel_with_weights ───────────────────────────────────────────────
run_blockmodel_with_weights <- function(M, max_k, relWeights,
                                        prespec    = NULL,
                                        approach   = "bin",
                                        mingr      = NULL,
                                        reps_kloop = REPS_KLOOP,
                                        reps_aop   = REPS_AOP) {
  R <- dim(M)[3]
  if (is.null(prespec)) prespec <- c(rep("comnul", R - 1L), "nuldnc")

  n <- dim(M)[1]
  if (max_k >= n) max_k <- n - 1L

  approaches      <- rep(approach, R)
  lnames          <- dimnames(M)[[3]]
  if (is.null(lnames)) lnames <- paste0("L", seq_len(R))
  layer_col_names <- paste0("err_", gsub("\\n.*", "", lnames))

  ks           <- seq_len(max_k)
  errs         <- setNames(numeric(max_k), ks)
  layer_errs   <- matrix(NA_real_, nrow = max_k, ncol = R,
                         dimnames = list(NULL, layer_col_names))

  for (i in ks) {
    k          <- i
    blocks_arr <- array(NA_character_, dim = c(1L, R, k, k))
    for (rel in seq_len(R))
      for (ri in seq_len(k))
        for (ci in seq_len(k))
          blocks_arr[1, rel, ri, ci] <- prespec_block(prespec[rel], ri == ci)

    res_k <- suppressWarnings(optRandomParC(
      M          = M,
      k          = k,
      rep        = reps_kloop,
      approaches = approaches,
      blocks     = blocks_arr,
      relWeights = relWeights,
      mingr      = mingr,
      seed       = SEED_KLOOP,
      nCores     = 0,
      printRep   = FALSE
    ))
    errs[i]         <- err(res_k)
    layer_errs[i, ] <- per_layer_err(M, clu(res_k), prespec, approach)
  }

  weighted_layer_errs <- sweep(layer_errs, 2, relWeights, `*`)
  colnames(weighted_layer_errs) <- gsub("^err_", "werr_", colnames(layer_errs))

  interleaved <- as.data.frame(do.call(cbind, lapply(seq_len(R), function(r)
    cbind(layer_errs[, r, drop = FALSE],
          weighted_layer_errs[, r, drop = FALSE])
  )))
  df <- cbind(data.frame(k = ks, err = round(errs, 4)), round(interleaved, 1))

  min_err <- min(df$err)
  k_best  <- min(df$k[abs(df$err - min_err) < 1e-8])

  # Full search at k_best: all distinct optimal partitions via max.iden + switch.names
  blocks_arr_best <- array(NA_character_, dim = c(1L, R, k_best, k_best))
  for (rel in seq_len(R))
    for (ri in seq_len(k_best))
      for (ci in seq_len(k_best))
        blocks_arr_best[1, rel, ri, ci] <- prespec_block(prespec[rel], ri == ci)

  all_partitions <- suppressWarnings(optRandomParC(
    M            = M,
    k            = k_best,
    rep          = reps_aop,
    approaches   = approaches,
    blocks       = blocks_arr_best,
    relWeights   = relWeights,
    mingr        = mingr,
    seed         = SEED_AOP,
    nCores       = 0,
    printRep     = FALSE,
    max.iden     = reps_aop,
    switch.names = TRUE
  ))

  clu0     <- clu(all_partitions$best[[1]])
  clusters <- sort(unique(clu0))
  clu_size <- setNames(sapply(clusters, function(g) sum(clu0 == g)), clusters)
  srt      <- sort_partition(clu0, rownames(M))
  map      <- setNames(seq_along(srt$cluster_order), srt$cluster_order)
  clu_ord  <- as.integer(map[as.character(clu0)])[srt$node_order]

  invisible(list(
    df             = df,
    k_best         = k_best,
    all_partitions = all_partitions,
    clu_size       = clu_size,
    node_order     = srt$node_order,
    clu_ord        = clu_ord,
    prettyM        = M[srt$node_order, srt$node_order, , drop = FALSE]
  ))
}

# ── plot_mat_with_labels ───────────────────────────────────────────────────────
# Replacement for plotMat that adds node names on both axes, cluster block
# labels above the columns, and a two-line header (layer name + block labels).
#
# mat        : square numeric matrix, already reordered by partition
# clu        : integer vector, cluster membership in display order (1-indexed)
# layer_name : character, title for this subplot (may contain \n)
# clu_labels : character vector length k, label per cluster (or NULL to skip)
plot_mat_with_labels <- function(mat, clu, layer_name,
                                 clu_labels     = NULL,
                                 highlight_type = "none") {
  n         <- nrow(mat)
  k         <- max(clu)
  nms       <- rownames(mat)
  clu_sizes <- sapply(seq_len(k), function(g) sum(clu == g))
  cum_sizes <- cumsum(clu_sizes)

  # image() plots y from bottom to top; rows should go top to bottom.
  # z[i, j] = mat value at column i, row j-from-bottom = mat[n+1-j, i]
  # → z = t(mat[n:1, ])
  z <- t(mat[n:1, , drop = FALSE])

  image(seq_len(n), seq_len(n), z,
        col   = c("white", "#111111"),
        breaks = c(-0.5, 0.5, 1.5),
        axes  = FALSE,
        xlab  = "", ylab  = "",
        main  = "")                    # title added via mtext below

  # Column labels at TOP (rotated, las=2) and row labels at LEFT (horizontal, las=1).
  ax_cex <- if (n <= 12) 0.80 else if (n <= 20) 0.70 else 0.60
  for (i in seq_len(n)) {
    mtext(nms[i],       side = 3, at = i, line = 0.3, cex = ax_cex, las = 2)
    mtext(rev(nms)[i],  side = 2, at = i, line = 0.3, cex = ax_cex, las = 1)
  }

  # Layer title above the column labels
  nchar_max  <- max(nchar(nms))
  title_line <- nchar_max * 0.62 + 2.2
  mtext(gsub("\\n", " ", layer_name), side = 3, line = title_line,
        cex = 0.90, font = 2)

  # Cluster boundary lines (blue)
  if (k > 1) {
    for (b in cum_sizes[-k]) {
      abline(v = b + 0.5,     col = "blue", lwd = 1.5)
      abline(h = n - b + 0.5, col = "blue", lwd = 1.5)
    }
  }

  # Coloured rect around matching diagonal blocks
  # M1 (highlight_type="function"): orange around machine/step-pure clusters
  # M2 (highlight_type="product"):  green  around product-pure clusters
  if (highlight_type != "none") {
    starts <- c(0L, cum_sizes[-k]) + 1L
    hl_col <- if (highlight_type == "function") "#e67e22" else "#27ae60"
    for (g in seq_len(k)) {
      members_g <- nms[clu == g]
      prods_g   <- unique(sub("-.*",  "", members_g))
      fns_g     <- unique(sub(".*-", "", members_g))
      do_hl <- (highlight_type == "function" &&
                  length(fns_g) == 1 && length(prods_g) > 1) ||
               (highlight_type == "product" &&
                  length(prods_g) == 1 && length(fns_g) > 1)
      if (do_hl) {
        rect(starts[g] - 0.5,    n - cum_sizes[g] + 0.5,
             cum_sizes[g] + 0.5, n - starts[g] + 1.5,
             border = hl_col, lwd = 8.0, col = NA)
      }
    }
  }

  # Cluster labels at BOTTOM (side 1), centred on each column block.
  if (!is.null(clu_labels) && length(clu_labels) == k) {
    starts_lbl <- c(0L, cum_sizes[-k]) + 1L
    centers    <- (starts_lbl + cum_sizes) / 2
    hl_col_lbl <- if (highlight_type == "function") "#e67e22"
                  else if (highlight_type == "product") "#27ae60"
                  else NA
    for (i in seq_len(k)) {
      members_i <- nms[clu == i]
      prods_i   <- unique(sub("-.*",  "", members_i))
      fns_i     <- unique(sub(".*-", "", members_i))
      use_col <- !is.na(hl_col_lbl) && (
        (highlight_type == "function" && length(fns_i)   == 1 && length(prods_i) > 1) ||
        (highlight_type == "product"  && length(prods_i) == 1 && length(fns_i)   > 1)
      )
      mtext(clu_labels[i], side = 1, at = centers[i],
            line = 1.8, cex = 1.00, font = 2,
            col  = if (use_col) hl_col_lbl else "black")
    }
  }

  box()
}

# ── cluster_label_short ────────────────────────────────────────────────────────
# Returns a short label for a set of node names: "Prod. A", "Mach. 3", or "Mixed".
cluster_label_short <- function(members) {
  prods <- unique(sub("-.*",  "", members))
  fns   <- unique(sub(".*-", "", members))
  if      (length(prods) == 1 && length(fns) > 1)  paste0("Prod. ", prods)
  else if (length(fns)   == 1 && length(prods) > 1) paste0("Mach. ", fns)
  else "Mixed"
}

# ── cluster_label_long ─────────────────────────────────────────────────────────
# Same but long form for the cluster table: "Product A", "Machine 3", "Mixed".
cluster_label_long <- function(members) {
  prods <- unique(sub("-.*",  "", members))
  fns   <- unique(sub(".*-", "", members))
  if      (length(prods) == 1 && length(fns) > 1)  paste0("Product ",  prods)
  else if (length(fns)   == 1 && length(prods) > 1) paste0("Machine ", fns)
  else "Mixed"
}

# ═══════════════════════════════════════════════════════════════════════════════
# UI
# ═══════════════════════════════════════════════════════════════════════════════
ui <- fluidPage(

  tags$head(tags$style(HTML("
    body { font-family: 'Helvetica Neue', Arial, sans-serif; font-size: 14px; }
    .sidebar { background: #f4f6f9; padding: 16px 12px; border-radius: 6px; }
    h4 { color: #2c3e50; margin-top: 18px; margin-bottom: 4px; font-size: 14px;
         font-weight: 600; text-transform: uppercase; letter-spacing: 0.04em; }
    .note { font-size: 12px; color: #777; margin-top: 2px; margin-bottom: 8px; }
    .result-banner { background: #e8f0fe; border-left: 4px solid #3b6bdb;
                     padding: 10px 14px; border-radius: 4px; margin-bottom: 14px; }
    .warn-banner   { background: #fff3cd; border-left: 4px solid #e6a817;
                     padding: 10px 14px; border-radius: 4px; margin-bottom: 14px; }
    table.shiny-table { font-size: 13px; }
    table.shiny-table thead th { background: #3b6bdb; color: white; }
    table.shiny-table tr.k-best td { background: #ddeeff; font-weight: bold; }
    hr { margin: 12px 0; border-color: #ddd; }
  "))),

  titlePanel("Multiplex Blockmodeling Explorer"),
  tags$p("Adjust the network structure and layer weights, then click Run to see which grouping the algorithm selects.",
         style = "color:#555; margin-bottom:16px;"),

  sidebarLayout(

    # ── Sidebar ──────────────────────────────────────────────────────────────
    sidebarPanel(
      width = 3,
      div(class = "sidebar",

        actionButton("run", "▶  Run blockmodel",
                     class = "btn-primary",
                     style = "width:100%; font-size:14px; padding:8px;"),

        hr(),

        h4("Layer weights"),
        p("Relative importance of each rationale. Equal weights = 1, 1, 1.", class = "note"),
        sliderInput("w1", "w₁  Machine-sharing (M1)",  min = 0, max = 3, value = 1, step = 0.1, ticks = FALSE),
        sliderInput("w2", "w₂  Product-sharing (M2)",  min = 0, max = 3, value = 1, step = 0.1, ticks = FALSE),
        sliderInput("w3", "w₃  Goal conflict (M3)",    min = 0, max = 3, value = 1, step = 0.1, ticks = FALSE),

        hr(),

        h4("Computation"),
        p("Higher values = more reliable but slower.", class = "note"),
        numericInput("reps_kloop", "Restarts per k (k-selection loop)",
                     value = 20, min = 5, max = 200, step = 5),
        numericInput("reps_aop",   "Restarts (optimal partition search at k*)",
                     value = 10, min = 5, max = 500, step = 5),
        numericInput("mingr",      "Min. cluster size (0 = auto-derive)",
                     value = 0,  min = 0, max = 20,  step = 1),

        hr(),

        h4("Simulation parameters"),
        sliderInput("N_products", "Products (A, B, …)",
                    min = 2, max = 8, value = 4, step = 1, ticks = FALSE),
        sliderInput("N_steps",   "Machines per product",
                    min = 2, max = 6, value = 4, step = 1, ticks = FALSE)
      )
    ),

    # ── Main panel ────────────────────────────────────────────────────────────
    mainPanel(
      width = 9,

      uiOutput("status_banner"),

      fluidRow(
        column(4,
          h4("Error by k"),
          p("★ = selected k (minimum weighted criterion).", class = "note"),
          tableOutput("err_table")
        ),
        column(8,
          h4("Cluster assignments"),
          p("First distinct optimal partition. If multiple exist, see note below.", class = "note"),
          tableOutput("cluster_table")
        )
      ),

      hr(),
      h4("Block matrices — optimal partition"),
      p("Adjacency matrices reordered by the optimal partition. Lines delimit clusters.",
        class = "note"),
      plotOutput("block_plot", height = "500px")
    )
  )
)

# ═══════════════════════════════════════════════════════════════════════════════
# Server
# ═══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  # ── Main computation — fires only on Run button ─────────────────────────────
  results <- eventReactive(input$run, {

    # Input validation
    validate(
      need(input$w1 + input$w2 + input$w3 > 0,
           "At least one layer weight must be > 0.")
    )

    withProgress(message = "Running blockmodel…", value = 0, {

      setProgress(0.1, detail = "Generating task matrices")
      mats <- generate_task_matrices(
        N_products = input$N_products,
        N_steps    = input$N_steps
      )

      # Stack into n × n × 3 array
      n <- nrow(mats$M1)
      M <- array(NA_real_, dim = c(n, n, 3),
                 dimnames = list(
                   rownames(mats$M1), colnames(mats$M1),
                   c("Machine-sharing\n(functional dependency)",
                     "Product-sharing\n(workflow dependency)",
                     "Goal conflict")
                 ))
      M[,,1] <- mats$M1
      M[,,2] <- mats$M2
      M[,,3] <- mats$M3

      # Normalised weights: each multiplier scales nw = 1 / n(n-1)
      nw         <- 1 / (n * (n - 1L))
      relWeights <- c(input$w1 * nw, input$w2 * nw, input$w3 * nw)

      # Derive max_k and mingr from data dimensions (same logic as sim_unimode.qmd)
      k_M1 <- input$N_steps
      k_M2 <- input$N_products
      mingr_auto <- max(1L, min(input$N_products - 1L, input$N_steps - 1L))
      mingr_used <- if (input$mingr == 0) mingr_auto else as.integer(input$mingr)
      max_k      <- min(n - 1L,
                        k_M1 + k_M2 + 2L,
                        floor(n / mingr_used))
      max_k      <- max(max_k, 1L)

      setProgress(0.25, detail = paste("Fitting k = 1 …", max_k))
      out <- run_blockmodel_with_weights(
        M          = M,
        max_k      = max_k,
        relWeights = relWeights,
        prespec    = c("comnul", "comnul", "nuldnc"),
        approach   = APPROACH,
        mingr      = mingr_used,
        reps_kloop = as.integer(input$reps_kloop),
        reps_aop   = as.integer(input$reps_aop)
      )

      setProgress(1, detail = "Done")
      list(out = out, M = M, nw = nw,
           w          = c(input$w1, input$w2, input$w3),
           mingr_used = mingr_used,
           mingr_auto = mingr_auto,
           reps_kloop = as.integer(input$reps_kloop),
           reps_aop   = as.integer(input$reps_aop))
    })
  })

  # ── Status banner ───────────────────────────────────────────────────────────
  output$status_banner <- renderUI({
    if (input$run == 0) {
      return(div(class = "result-banner",
        "Set parameters and click ", strong("Run blockmodel"), " to see results."))
    }
    req(results())
    out     <- results()$out
    n_parts <- length(out$all_partitions$best)
    w       <- results()$w

    banner_class <- if (n_parts > 1) "warn-banner" else "result-banner"
    part_note    <- if (n_parts > 1)
      paste0(" ⚠ ", n_parts, " distinct optimal partitions — showing the first. ",
             "Consider increasing reps or reviewing all partitions.")
    else " · 1 distinct optimal partition."

    div(class = banner_class,
        strong(paste0("Optimal k = ", out$k_best)),
        paste0(" · Criterion value = ", round(err(out$all_partitions), 4)),
        part_note,
        br(),
        tags$em(paste0("Weights: M1 = ", w[1], " · M2 = ", w[2],
                        " · M3 = ", w[3], "  (× normalised baseline nw)")),
        br(),
        tags$em(paste0(
          "Computation: k-loop ", results()$reps_kloop, " restarts · ",
          "AOP ", results()$reps_aop, " restarts · ",
          "mingr = ", results()$mingr_used,
          if (input$mingr == 0) paste0(" (auto; derived from ", results()$mingr_auto, ")") else ""
        ))
    )
  })

  # ── Error table ─────────────────────────────────────────────────────────────
  output$err_table <- renderTable({
    req(results())
    df     <- results()$out$df
    k_best <- results()$out$k_best

    # Show k, total weighted error, raw per-layer errors
    raw_cols <- grep("^err_", names(df), value = TRUE)
    df_show  <- df[, c("k", "err", raw_cols), drop = FALSE]

    # Clean column names
    layer_labels <- gsub("err_", "", raw_cols)
    layer_labels <- gsub("\\.", " ", layer_labels)
    names(df_show) <- c("k", "Total (wtd)", layer_labels)

    # Mark the selected k
    df_show$k <- ifelse(df_show$k == k_best,
                        paste0("★ ", df_show$k),
                        as.character(df_show$k))
    df_show
  },
  striped  = TRUE,
  hover    = TRUE,
  bordered = TRUE,
  digits   = 2,
  rownames = FALSE)

  # ── Cluster assignment table ─────────────────────────────────────────────────
  output$cluster_table <- renderTable({
    req(results())
    out   <- results()$out
    clu0  <- clu(out$all_partitions$best[[1]])
    nodes <- rownames(results()$M)
    k     <- out$k_best

    do.call(rbind, lapply(seq_len(k), function(g) {
      members <- sort(nodes[clu0 == g])
      data.frame(
        Cluster = paste0("C", g),
        Label   = cluster_label_long(members),
        Size    = length(members),
        Tasks   = paste(members, collapse = ", "),
        stringsAsFactors = FALSE
      )
    }))
  },
  striped  = TRUE,
  hover    = TRUE,
  bordered = TRUE,
  rownames = FALSE)

  # ── Block matrix plot ────────────────────────────────────────────────────────
  output$block_plot <- renderPlot({
    req(results())
    out <- results()$out
    M   <- results()$M
    R   <- dim(M)[3]

    clu0   <- clu(out$all_partitions$best[[1]])
    srt    <- sort_partition(clu0, rownames(M))
    node_t <- srt$node_order
    map_t  <- setNames(seq_along(srt$cluster_order), srt$cluster_order)
    clu_t  <- as.integer(map_t[as.character(clu0)])[node_t]
    mats_t <- M[node_t, node_t, , drop = FALSE]
    nms_t  <- rownames(M)[node_t]
    k      <- max(clu_t)

    # Short cluster labels ("Prod. A", "Mach. 1", "Mixed") used above column blocks
    clu_labels <- sapply(seq_len(k), function(g)
      cluster_label_short(nms_t[clu_t == g])
    )

    max_nchar <- max(nchar(nms_t))
    top_mar   <- max_nchar * 0.62 + 4.0
    lft_mar   <- max_nchar * 0.75 + 1.5
    bot_mar   <- 4.0

    # highlight_type per layer: M1 highlights machine/step clusters,
    # M2 highlights product clusters, M3 none
    hl_types <- c("function", "product", "none")

    par(mfrow = c(1, R), mar = c(bot_mar, lft_mar, top_mar, 0.5), pty = "s")
    lnames <- dimnames(M)[[3]]
    for (r in seq_len(R)) {
      mat_r <- mats_t[,,r]
      rownames(mat_r) <- nms_t
      colnames(mat_r) <- nms_t
      hl <- if (r <= length(hl_types)) hl_types[r] else "none"
      plot_mat_with_labels(mat_r, clu_t, lnames[r], clu_labels,
                           highlight_type = hl)
    }
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
  })
}

# ── Launch ────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
