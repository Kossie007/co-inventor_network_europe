# ===================================================================== --
# 0. General settings ----
# ===================================================================== --

# Reproducibility
set.seed(20251202)
{  req <- c(
    "cluster",
    "data.table",
    "dendextend",
    "dplyr",
    "ergm",
    "factoextra",
    "ggnewscale",
    "ggplot2",
    "ggraph",
    "ggrepel",
    "giscoR",
    "httr",
    "igraph",
    "network",
    "png",
    "purrr",
    "rlang",
    "scales",
    "sf",
    "this.path",
    "tibble",
    "tidyr",
    "viridis"
  )
  
  to_install <- setdiff(req, rownames(installed.packages()))
  if (length(to_install)) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(req, library, character.only = TRUE))
  
  # Set Working Directory
  setwd(this.path::here())
  rm(req, to_install)
  print(getwd())
}

# ===================================================================== --
# 1. Regions: EU subset + lookup ----
# ===================================================================== --

# Load regions data
regions <- fread(
  "REGPAT_REGIONS.txt",
  sep = "|",
  header = TRUE,
  encoding = "UTF-8"
) %>%
  as_tibble() %>%
  rename_with(tolower)
# -> ctry_code, reg_code, reg_label, up_reg_code, up_reg_label

# European countries (EPO member states + UK, NO, IS, CH)
eu_ctry <- c(
  "AT", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "ES",
  "FI", "FR", "GB", "GR", "HR", "HU", "IE", "IS", "IT", "LT",
  "LU", "LV", "MT", "NL", "NO", "PL", "PT", "RO", "SE", "SI", "SK"
)

# Filter for EU regions
regions_eu <- regions %>%
  filter(ctry_code %in% eu_ctry)

eu_nuts3_codes <- unique(regions_eu$reg_code) # 1390 entries
eu_nuts2_codes <- unique(regions_eu$up_reg_code) # 285 entries
eu_ctry_codes  <- unique(regions_eu$ctry_code) # 31 entries



# ===================================================================== --
# 2. EPO ----
# ===================================================================== --

# Load EPO inventor data
inv_epo <- fread(
  "202401_EPO_Inv_reg.txt",
  sep = "|",
  header = TRUE,
  encoding = "UTF-8"
)

# Process inventors: Filter for EU countries and NUTS3 regions
inv_epo_eu <- as_tibble(inv_epo) %>%
  filter(
    ctry_code %in% eu_ctry_codes,
    reg_code %in% eu_nuts3_codes
  )

# EU-relevant patent IDs
epo_eu_appln_ids <- unique(inv_epo_eu$appln_id) # EPO applications with at least one EU inventor

rm(inv_epo, epo_eu_appln_ids)
gc()

# IPC categorization
ipc_epo <- fread(
  "202401_EPO_IPC.txt",
  sep = "|",
  header = TRUE,
  encoding = "UTF-8"
)

ipc_epo_dt <- as_tibble(ipc_epo)

# Join IPC data and create ipc4 column
inv_epo_eu <- inv_epo_eu %>%
  left_join(
    ipc_epo_dt %>% select(appln_id, app_year, IPC),
    by = "appln_id"
  ) %>%
  mutate(
    ipc4 = substr(IPC, 1, 4) # section + class + subclass
  )

# Select common columns
common_cols <- c(
  "appln_id",
  "person_id",
  "reg_code",
  "ctry_code",
  "reg_share",
  "inv_share",
  "app_year",
  "ipc4"
)

inv_epo_eu <- inv_epo_eu %>%
  select(all_of(common_cols))

# Unique inventor–patent links
inv_epo_link <- unique(inv_epo_eu)

# Cleanup

gc()


# ===================================================================== --
# 3. Edge table: co-inventor pairs ----
# ===================================================================== --

# Unique inventor–patent pairs
inv_pat <- inv_epo_link %>%
  filter(!is.na(person_id)) %>%
  distinct(appln_id, person_id)

# Inventors per patent
pat_inv_counts <- inv_pat %>%
  count(appln_id, name = "n_inv")

# Keep only patents with n_inv >= 2
multi_pat <- inv_pat %>%
  inner_join(
    pat_inv_counts %>% filter(n_inv >= 2),
    by = "appln_id"
  )

# Co-inventor pairs within each patent
coinv_el <- multi_pat %>%
  select(appln_id, person_id) %>%
  inner_join(
    .,
    .,
    by = "appln_id",
    suffix = c("_from", "_to")
  ) %>%
  filter(person_id_from < person_id_to) %>% # unordered pairs, no self-pairs
  count(person_id_from, person_id_to, name = "weight") %>% # weight = shared patents
  rename(
    from = person_id_from,
    to   = person_id_to
  )

# Node list: all inventors that ever appear
inventor_nodes <- inv_epo_link %>%
  filter(!is.na(person_id)) %>%
  distinct(person_id) %>%
  rename(name = person_id)

# Create graph object
coinv_net <- graph_from_data_frame(
  d        = coinv_el,
  directed = FALSE,
  vertices = inventor_nodes
)

# Basic sanity checks
# vcount(coinv_net) # number of inventors
# ecount(coinv_net) # number of co-inventor links
# sum(degree(coinv_net) == 0) # isolates

# Component structure
comp <- components(coinv_net)
# table(comp$csize)[1:10]
# max(comp$csize)

# Cleanup

gc()


# ===================================================================== --
# 4. Uzzi type Z scores ----
# ===================================================================== --

tech_df <- inv_epo_eu %>%
  filter(!is.na(ipc4), !is.na(app_year)) %>%
  distinct(
    pat_id = appln_id,
    app_year,
    ipc4
  )

# 1) Time windows
min_year <- 1978L
max_year <- 2023L

starts <- seq(min_year, 2018L, by = 5L)
ends   <- c(starts[-1] - 1L, max_year)

windows <- tibble(
  period_id = seq_along(starts),
  start_y   = starts,
  end_y     = ends
)

# Function to compute windowed Uzzi scores
compute_window_inventor_uzzi <- function(start_y, end_y, period_id,
                                         tech_df, inv_epo_eu) {
  # 1) Patent–IPC for windows
  win_codes <- tech_df %>%
    filter(app_year >= start_y, app_year <= end_y) %>%
    distinct(pat_id, ipc4)
  
  if (nrow(win_codes) == 0L) return(NULL)
  
  # >= 2 IPC code patents
  pat_sizes <- win_codes %>%
    count(pat_id, name = "n_codes")
  
  win_codes <- win_codes %>%
    inner_join(
      pat_sizes %>% filter(n_codes >= 2L),
      by = "pat_id"
    )
  
  N <- dplyr::n_distinct(win_codes$pat_id)
  if (N < 2L) return(NULL)
  
  # 2) n_i: number of patents per IPC
  code_counts <- win_codes %>%
    count(ipc4, name = "n_i")
  
  # 3) As unordered IPC-pair per patent
  pairs_df <- win_codes %>%
    inner_join(
      win_codes,
      by = "pat_id",
      relationship = "many-to-many"
    ) %>%
    filter(ipc4.x < ipc4.y) # i < j, no self-pair
  
  if (nrow(pairs_df) == 0L) return(NULL)
  
  # 4) Observed co-occurence O_ij + n_i, n_j
  O_ij <- pairs_df %>%
    count(ipc4.x, ipc4.y, name = "O_obs") %>%
    rename(
      ipc_i = ipc4.x,
      ipc_j = ipc4.y
    ) %>%
    left_join(
      code_counts %>% rename(ipc_i = ipc4, n_i = n_i),
      by = "ipc_i"
    ) %>%
    left_join(
      code_counts %>% rename(ipc_j = ipc4, n_j = n_i),
      by = "ipc_j"
    ) %>%
    mutate(
      E_ij      = (n_i * n_j) / N,
      sigma2_ij = E_ij * (1 - n_i / N) * ((N - n_j) / (N - 1)),
      sigma_ij  = sqrt(pmax(sigma2_ij, 0)),
      Z_ij      = dplyr::if_else(
        sigma_ij > 0,
        (O_obs - E_ij) / sigma_ij,
        NA_real_
      )
    )
  
  if (all(is.na(O_ij$Z_ij))) return(NULL)
  
  # 5) Uzzi: bottom 10% Z (most negative) = novel pair
  z_cut_10 <- stats::quantile(O_ij$Z_ij, probs = 0.10, na.rm = TRUE)
  
  O_ij <- O_ij %>%
    mutate(
      novel_pair_uzzi = !is.na(Z_ij) & Z_ij <= z_cut_10
    )
  
  # 6) Joining back patent–IPC pairs
  pairs_z <- pairs_df %>%
    left_join(
      O_ij %>% select(ipc_i, ipc_j, Z_ij, novel_pair_uzzi),
      by = c("ipc4.x" = "ipc_i", "ipc4.y" = "ipc_j")
    )
  
  # 7) Patent levels: atypical if there is at least 1 novel IPC-pair
  pat_atyp <- pairs_z %>%
    filter(novel_pair_uzzi) %>%
    distinct(pat_id) %>%
    mutate(ATYPICAL = TRUE)
  
  # 8) Inventor–patent windows (EU inventors)
  inv_pat_win <- inv_epo_eu %>%
    filter(app_year >= start_y, app_year <= end_y) %>%
    transmute(
      pat_id    = appln_id,
      person_id = person_id
    ) %>%
    distinct() %>%
    filter(pat_id %in% win_codes$pat_id) %>%
    left_join(pat_atyp, by = "pat_id") %>%
    mutate(
      ATYPICAL = dplyr::if_else(is.na(ATYPICAL), FALSE, ATYPICAL)
    )
  
  if (nrow(inv_pat_win) == 0L) return(NULL)
  
  # 9) Inventor x period aggregation
  inv_period <- inv_pat_win %>%
    group_by(person_id) %>%
    summarise(
      n_pat      = dplyr::n_distinct(pat_id),
      n_atyp     = sum(ATYPICAL),
      share_atyp = mean(ATYPICAL),
      .groups    = "drop"
    ) %>%
    mutate(
      period_id    = period_id,
      period_start = start_y,
      period_end   = end_y
    )
  
  inv_period
}

# Calculate Uzzi scores
inventor_uzzi_list <- purrr::pmap(
  list(windows$start_y, windows$end_y, windows$period_id),
  ~ compute_window_inventor_uzzi(..1, ..2, ..3, tech_df, inv_epo_eu)
)

inventor_uzzi_all <- bind_rows(inventor_uzzi_list)

# Aggregated overall inventor-level Uzzi
inventor_uzzi_overall <- inventor_uzzi_all %>%
  group_by(person_id) %>%
  summarise(
    n_pat_tot      = sum(n_pat),
    n_atyp_tot     = sum(n_atyp),
    share_atyp_tot = n_atyp_tot / n_pat_tot,
    .groups        = "drop"
  )


gc()

# ==================================================================== --
# 5. Inventor-level centralities & Burt's constraint ----
# ==================================================================== --

inventor_cent <- tibble(
  person_id = as.numeric(V(coinv_net)$name),
  # degree: number of co-inventors
  deg       = degree(coinv_net, mode = "all"),
  # k-core index: how deep in the core
  kcore     = coreness(coinv_net),
  # Burt's constraint: brokerage (lower = more broker-like)
  constraint = constraint(coinv_net)
)

# Keep only meaningful cases for brokerage summaries (deg >= 2)
inventor_cent_nz <- inventor_cent %>%
  filter(deg >= 2, !is.na(constraint))

# Mean and sd of Burt's constraint among non-isolates
mu_con <- mean(inventor_cent_nz$constraint, na.rm = TRUE)
sd_con <- sd(inventor_cent_nz$constraint,   na.rm = TRUE)

# Add z-score of constraint (lower z = more brokerage)
inventor_cent <- inventor_cent %>%
  mutate(
    constraint_z = (constraint - mu_con) / sd_con
  )


# ===================================================================== --
# 6. Broker flags based on mean and 1 SD ----
# ===================================================================== --

inventor_cent <- inventor_cent %>%
  mutate(
    # "Broker below mean": constraint lower than mean (z < 0), deg >= 2
    is_broker_below_mean = if_else(
      deg >= 2 & !is.na(constraint_z) & constraint_z < 0,
      TRUE,
      FALSE
    ),
    # "Strong broker": at least 1 SD below mean (z <= -1), deg >= 2
    is_broker_1sd = if_else(
      deg >= 2 & !is.na(constraint_z) & constraint_z <= -1,
      TRUE,
      FALSE
    )
  )

# Quick sanity checks
table(inventor_cent$is_broker_below_mean, useNA = "ifany")
table(inventor_cent$is_broker_1sd,        useNA = "ifany")

# Top strong brokers (1 SD below mean)
top_brokers_1sd <- inventor_cent %>%
  filter(is_broker_1sd == TRUE) %>%
  arrange(constraint) %>%
  slice_head(n = 50)

# Correlation table
inventor_cent_nz <- inventor_cent %>%
  filter(deg >= 2, !is.na(constraint))

cor(
  inventor_cent_nz %>% select(deg, kcore, constraint),
  use    = "pairwise.complete.obs",
  method = "pearson"
)

# Top hubs
top_hubs <- inventor_cent %>%
  arrange(desc(deg)) %>%
  slice_head(n = 50)
gc()


# ===================================================================== --
# 7. Descriptive statistics: brokers vs non-brokers ----
# ===================================================================== --

# Merge constraint/centrality with Uzzi-type inventor novelty measures
inventor_panel <- inventor_cent %>%
  left_join(inventor_uzzi_overall, by = "person_id")

# For brokerage, we only compare inventors with deg >= 2 and non-missing constraint
inventor_panel_use <- inventor_panel %>%
  filter(deg >= 2, !is.na(constraint))

# Brokers based on "below mean" constraint (z < 0)
desc_below_mean <- inventor_panel_use %>%
  mutate(group_bm = if_else(
    is_broker_below_mean,
    "broker_below_mean",
    "non_broker_below_mean"
  )) %>%
  group_by(group_bm) %>%
  summarise(
    n_inventors       = n(),
    share_of_sample   = n() / nrow(inventor_panel_use),
    
    # network position
    mean_deg          = mean(deg, na.rm = TRUE),
    median_deg        = median(deg, na.rm = TRUE),
    mean_kcore        = mean(kcore, na.rm = TRUE),
    median_kcore      = median(kcore, na.rm = TRUE),
    mean_constraint   = mean(constraint, na.rm = TRUE),
    sd_constraint     = sd(constraint, na.rm = TRUE),
    
    # innovation profile (Uzzi-based atypicality)
    mean_n_pat_tot    = mean(n_pat_tot, na.rm = TRUE),
    median_n_pat_tot  = median(n_pat_tot, na.rm = TRUE),
    mean_share_atyp   = mean(share_atyp_tot, na.rm = TRUE),
    median_share_atyp = median(share_atyp_tot, na.rm = TRUE),
    
    .groups = "drop"
  )

print(desc_below_mean)

gc()

# Strong brokers: constraint <= mean - 1 SD (z <= -1)
desc_1sd <- inventor_panel_use %>%
  mutate(group_1sd = if_else(
    is_broker_1sd,
    "strong_broker_1sd",
    "non_strong_broker"
  )) %>%
  group_by(group_1sd) %>%
  summarise(
    n_inventors       = n(),
    share_of_sample   = n() / nrow(inventor_panel_use),
    
    # network position
    mean_deg          = mean(deg, na.rm = TRUE),
    median_deg        = median(deg, na.rm = TRUE),
    mean_kcore        = mean(kcore, na.rm = TRUE),
    median_kcore      = median(kcore, na.rm = TRUE),
    mean_constraint   = mean(constraint, na.rm = TRUE),
    sd_constraint     = sd(constraint, na.rm = TRUE),
    
    # innovation profile (Uzzi-based atypicality)
    mean_n_pat_tot    = mean(n_pat_tot, na.rm = TRUE),
    median_n_pat_tot  = median(n_pat_tot, na.rm = TRUE),
    mean_share_atyp   = mean(share_atyp_tot, na.rm = TRUE),
    median_share_atyp = median(share_atyp_tot, na.rm = TRUE),
    
    .groups = "drop"
  )

print(desc_1sd)

# T-tests for difference in means
t.test(
  share_atyp_tot ~ is_broker_1sd,
  data = inventor_panel_use
)

lm(
  share_atyp_tot ~ is_broker_1sd + log1p(n_pat_tot) + deg + kcore,
  data = inventor_panel_use
) %>%
  summary()

gc()

# ===================================================================== --
# 8. Geographical distribution of brokers ----
# ===================================================================== --

inv_loc <- inv_epo_link %>%
  filter(
    !is.na(person_id),
    !is.na(ctry_code)
  ) %>%
  distinct(person_id, appln_id, ctry_code)

# Main country = country with most patents for that inventor
inv_main_ctry <- inv_loc %>%
  count(person_id, ctry_code, name = "n_pat_ctry") %>%
  group_by(person_id) %>%
  slice_max(order_by = n_pat_ctry, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(main_ctry = ctry_code)

# Join main country to inventor centrality table
inventor_geo <- inventor_cent %>%
  left_join(inv_main_ctry, by = "person_id") %>%
  filter(!is.na(main_ctry))

# Country-level broker counts
ctry_stats <- inventor_geo %>%
  group_by(main_ctry) %>%
  summarise(
    n_inv               = n(),
    n_broker_below_mean = sum(is_broker_below_mean, na.rm = TRUE),
    n_broker_1sd        = sum(is_broker_1sd,        na.rm = TRUE),
    .groups             = "drop"
  )

# Totals for Lorenz curves
tot_inv        <- sum(ctry_stats$n_inv)
tot_brok_bmean <- sum(ctry_stats$n_broker_below_mean)
tot_brok_1sd   <- sum(ctry_stats$n_broker_1sd)

# Lorenz-style cumulative shares
ctry_lorenz <- ctry_stats %>%
  arrange(n_inv) %>% # order by inventor count
  mutate(
    c_country      = row_number() / n(), # cumulative share of countries
    cum_inv        = cumsum(n_inv) / tot_inv, # cumulative share of inventors
    cum_brok_bmean = if (tot_brok_bmean > 0) {
      cumsum(n_broker_below_mean) / tot_brok_bmean
    } else {
      NA_real_
    },
    cum_brok_1sd   = if (tot_brok_1sd > 0) {
      cumsum(n_broker_1sd) / tot_brok_1sd
    } else {
      NA_real_
    }
  )

lorenz_long <- ctry_lorenz %>%
  select(
    main_ctry,
    c_country,
    cum_inv,
    cum_brok_bmean,
    cum_brok_1sd
  ) %>%
  pivot_longer(
    cols      = c(cum_inv, cum_brok_bmean, cum_brok_1sd),
    names_to  = "series",
    values_to = "cum_share"
  ) %>%
  mutate(
    series = dplyr::recode(
      series,
      "cum_inv"        = "Inventors",
      "cum_brok_bmean" = "Brokers: constraint < mean",
      "cum_brok_1sd"   = "Brokers: constraint \u2264 mean - 1 SD"
    )
  )
gc()

# Lorenz curves plot
ggplot(lorenz_long, aes(x = c_country, y = cum_share, color = series)) +
  geom_abline(
    slope     = 1,
    intercept = 0,
    linetype  = "dashed",
    linewidth = 0.7,
    colour    = "grey60"
  ) +
  geom_line(linewidth = 1.1, alpha = 0.9, na.rm = TRUE) +
  scale_color_manual(
    values = c(
      "Inventors"                         = "grey30",
      "Brokers: constraint < mean"        = "steelblue4",
      "Brokers: constraint \u2264 mean - 1 SD" = "firebrick4"
    )
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 5),
    name   = "Cumulative share of countries (ordered by inventor count)"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 5),
    name   = "Cumulative share of inventors / brokers"
  ) +
  coord_equal() +
  labs(
    title    = "Lorenz curves of inventors and brokers across countries",
    subtitle = "If brokers were proportional to inventors, broker curves would track the inventor curve",
    color    = NULL,
    caption  = "Source: own calculations based on OECD REGPAT data"
  ) +
  theme_minimal(base_family = "serif", base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 8)),
    axis.title    = element_text(size = 11),
    axis.text     = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    legend.text     = element_text(size = 10),
    plot.caption.position = "plot",
    plot.caption = element_text(
      hjust  = 1,
      face   = "italic",
      size   = 9,
      margin = margin(t = 10)
    ),
    plot.margin = margin(10, 15, 10, 10)
  )

# Gini coefficients function
gini_simple <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0) return(NA_real_)
  mu <- mean(x)
  if (mu == 0) return(NA_real_)
  diff_sum <- sum(abs(outer(x, x, "-")))
  diff_sum / (2 * n^2 * mu)
}

gini_inv        <- gini_simple(ctry_stats$n_inv)
gini_brok_bmean <- gini_simple(ctry_stats$n_broker_below_mean)
gini_brok_1sd   <- gini_simple(ctry_stats$n_broker_1sd)

gini_vals <- c(
  Gini_inventors               = gini_inv,
  Gini_brokers_constraint_mean = gini_brok_bmean,
  Gini_brokers_constraint_1sd  = gini_brok_1sd
)

print(gini_vals)

gc()

# ===================================================================== --
# 9. Technological distribution of brokers (by IPC4) ----
# ===================================================================== --

# Inventor–IPC mapping (main IPC4 per inventor)
inv_tech <- inv_epo_link %>%
  filter(
    !is.na(person_id),
    !is.na(ipc4)
  ) %>%
  distinct(person_id, appln_id, ipc4)

# Main IPC4 = IPC where the inventor has most patents
inv_main_ipc <- inv_tech %>%
  count(person_id, ipc4, name = "n_pat_ipc") %>%
  group_by(person_id) %>%
  slice_max(order_by = n_pat_ipc, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(main_ipc4 = ipc4)

# Join main IPC to inventor centrality table
inventor_tech <- inventor_cent %>%
  left_join(inv_main_ipc, by = "person_id") %>%
  filter(!is.na(main_ipc4))

# IPC-level broker counts
ipc_stats <- inventor_tech %>%
  group_by(main_ipc4) %>%
  summarise(
    n_inv               = n(),
    n_broker_below_mean = sum(is_broker_below_mean, na.rm = TRUE),
    n_broker_1sd        = sum(is_broker_1sd,        na.rm = TRUE),
    .groups             = "drop"
  )

# Totals for normalisation
tot_inv_ipc        <- sum(ipc_stats$n_inv)
tot_brok_bmean_ipc <- sum(ipc_stats$n_broker_below_mean)
tot_brok_1sd_ipc   <- sum(ipc_stats$n_broker_1sd)

# Lorenz-style cumulative shares across IPC4 classes
ipc_lorenz <- ipc_stats %>%
  arrange(n_inv) %>%
  mutate(
    c_ipc          = row_number() / n(), # cumulative share of IPC fields
    cum_inv        = cumsum(n_inv) / tot_inv_ipc, # inventors
    cum_brok_bmean = if (tot_brok_bmean_ipc > 0) {
      cumsum(n_broker_below_mean) / tot_brok_bmean_ipc
    } else {
      NA_real_
    },
    cum_brok_1sd   = if (tot_brok_1sd_ipc > 0) {
      cumsum(n_broker_1sd) / tot_brok_1sd_ipc
    } else {
      NA_real_
    }
  )

lorenz_ipc_long <- ipc_lorenz %>%
  select(
    main_ipc4,
    c_ipc,
    cum_inv,
    cum_brok_bmean,
    cum_brok_1sd
  ) %>%
  pivot_longer(
    cols      = c(cum_inv, cum_brok_bmean, cum_brok_1sd),
    names_to  = "series",
    values_to = "cum_share"
  ) %>%
  mutate(
    series = dplyr::recode(
      series,
      "cum_inv"        = "Inventors",
      "cum_brok_bmean" = "Brokers: constraint < mean",
      "cum_brok_1sd"   = "Brokers: constraint \u2264 mean - 1 SD"
    )
  )

# Lorenz curves plot for IPC classes
ggplot(lorenz_ipc_long, aes(x = c_ipc, y = cum_share, color = series)) +
  geom_abline(
    slope     = 1,
    intercept = 0,
    linetype  = "dashed",
    linewidth = 0.7,
    colour    = "grey60"
  ) +
  geom_line(linewidth = 1.1, alpha = 0.9, na.rm = TRUE) +
  scale_color_manual(
    values = c(
      "Inventors"                         = "grey30",
      "Brokers: constraint < mean"        = "steelblue4",
      "Brokers: constraint \u2264 mean - 1 SD" = "firebrick4"
    )
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 5),
    name   = "Cumulative share of IPC4 fields (ordered by inventor count)"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 5),
    name   = "Cumulative share of inventors / brokers"
  ) +
  coord_equal() +
  labs(
    title    = "Lorenz curves of inventors and brokers across IPC4 fields",
    subtitle = "If brokers were proportional to inventors, broker curves would track the inventor curve",
    color    = NULL,
    caption  = "Source: own calculations based on OECD REGPAT data"
  ) +
  theme_minimal(base_family = "serif", base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 8)),
    axis.title    = element_text(size = 11),
    axis.text     = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    legend.text     = element_text(size = 10),
    plot.caption.position = "plot",
    plot.caption = element_text(
      hjust  = 1,
      face   = "italic",
      size   = 9,
      margin = margin(t = 10)
    ),
    plot.margin = margin(10, 15, 10, 10)
  )

# Gini coefficients for concentration across IPC4
ipc_gini_inv        <- gini_simple(ipc_stats$n_inv)
ipc_gini_brok_bmean <- gini_simple(ipc_stats$n_broker_below_mean)
ipc_gini_brok_1sd   <- gini_simple(ipc_stats$n_broker_1sd)

ipc_gini_vals <- c(
  IPC_Gini_inventors               = ipc_gini_inv,
  IPC_Gini_brokers_constraint_mean = ipc_gini_brok_bmean,
  IPC_Gini_brokers_constraint_1sd  = ipc_gini_brok_1sd
)

print(ipc_gini_vals)

gc()

# ===================================================================== --
# 10. Elite clubs ----
# ===================================================================== --

# Broker-only subgraph + centralities in the broker-broker network
inventor_cent_tbl <- as_tibble(inventor_cent)

# Broker lists (two definitions)
brokers_bmean <- inventor_cent_tbl %>%
  filter(is_broker_below_mean, deg >= 2) %>%
  pull(person_id)

brokers_1sd <- inventor_cent_tbl %>%
  filter(is_broker_1sd, deg >= 2) %>%
  pull(person_id)

# Broker-only induced subgraphs
coinv_brok_bmean <- induced_subgraph(
  coinv_net,
  vids = V(coinv_net)[name %in% brokers_bmean]
)

coinv_brok_1sd <- induced_subgraph(
  coinv_net,
  vids = V(coinv_net)[name %in% brokers_1sd]
)

# Centrality measures *within the broker-broker network*
brok_bmean_cent <- tibble(
  person_id        = as.numeric(V(coinv_brok_bmean)$name),
  deg_broker_net   = degree(coinv_brok_bmean, mode = "all"),
  kcore_broker_net = coreness(coinv_brok_bmean),
  clustering_broker_net = transitivity(
    coinv_brok_bmean,
    type     = "local",
    isolates = "zero"
  )
)

brok_1sd_cent <- tibble(
  person_id        = as.numeric(V(coinv_brok_1sd)$name),
  deg_broker_net   = degree(coinv_brok_1sd, mode = "all"),
  kcore_broker_net = coreness(coinv_brok_1sd),
  clustering_broker_net = transitivity(
    coinv_brok_1sd,
    type     = "local",
    isolates = "zero"
  )
)

# Merging back to full inventor centrality table
inventor_cent_tbl <- inventor_cent_tbl %>%
  left_join(
    brok_bmean_cent %>%
      rename(
        deg_broker_net_bmean        = deg_broker_net,
        kcore_broker_net_bmean      = kcore_broker_net,
        clustering_broker_net_bmean = clustering_broker_net
      ),
    by = "person_id"
  ) %>%
  left_join(
    brok_1sd_cent %>%
      rename(
        deg_broker_net_1sd        = deg_broker_net,
        kcore_broker_net_1sd      = kcore_broker_net,
        clustering_broker_net_1sd = clustering_broker_net
      ),
    by = "person_id"
  )

# Edges from igraph to tibble
edges_df <- igraph::as_data_frame(coinv_net, what = "edges") %>%
  as_tibble() %>%
  mutate(
    from = as.character(from),
    to   = as.character(to)
  )

# Brokerage ends with broker flags
broker_flags <- inventor_cent_tbl %>%
  select(person_id, is_broker_below_mean, is_broker_1sd) %>%
  mutate(person_id = as.character(person_id))

edges_with_flags <- edges_df %>%
  left_join(
    broker_flags,
    by = c("from" = "person_id")
  ) %>%
  rename(
    is_broker_below_from = is_broker_below_mean,
    is_broker_1sd_from   = is_broker_1sd
  ) %>%
  left_join(
    broker_flags,
    by = c("to" = "person_id")
  ) %>%
  rename(
    is_broker_below_to = is_broker_below_mean,
    is_broker_1sd_to   = is_broker_1sd
  )

# Neighbour composition for 1 SD brokers
broker_neighbor_1sd <- edges_with_flags %>%
  mutate(
    broker_from = is_broker_1sd_from %in% TRUE,
    broker_to   = is_broker_1sd_to   %in% TRUE
  ) %>%
  # handled as egos
  pivot_longer(
    cols = c(from, to),
    names_to = "end",
    values_to = "person_id"
  ) %>%
  mutate(
    is_broker_self = dplyr::if_else(
      end == "from", broker_from, broker_to
    ),
    is_broker_nei = dplyr::if_else(
      end == "from", broker_to, broker_from
    )
  ) %>%
  group_by(person_id) %>%
  summarise(
    deg_total              = n(), # degree full network
    n_broker_neigh_1sd     = sum(is_broker_nei, na.rm = TRUE),
    share_broker_neigh_1sd = n_broker_neigh_1sd / deg_total,
    .groups = "drop"
  ) %>%
  left_join(
    inventor_cent_tbl %>%
      select(person_id, is_broker_1sd) %>%
      mutate(person_id = as.character(person_id)),
    by = "person_id"
  )

# Mean-based neighbour composition
broker_neighbor_bmean <- edges_with_flags %>%
  mutate(
    broker_from = is_broker_below_from %in% TRUE,
    broker_to   = is_broker_below_to   %in% TRUE
  ) %>%
  pivot_longer(
    cols = c(from, to),
    names_to = "end",
    values_to = "person_id"
  ) %>%
  mutate(
    is_broker_self = dplyr::if_else(
      end == "from", broker_from, broker_to
    ),
    is_broker_nei = dplyr::if_else(
      end == "from", broker_to, broker_from
    )
  ) %>%
  group_by(person_id) %>%
  summarise(
    deg_total                = n(),
    n_broker_neigh_bmean     = sum(is_broker_nei, na.rm = TRUE),
    share_broker_neigh_bmean = n_broker_neigh_bmean / deg_total,
    .groups = "drop"
  ) %>%
  left_join(
    inventor_cent_tbl %>%
      select(person_id, is_broker_below_mean) %>%
      mutate(person_id = as.character(person_id)),
    by = "person_id"
  )

# Only neighbors of 1 SD brokers
broker_neighbor_1sd %>%
  filter(is_broker_1sd) %>%
  summarise(
    mean_share_broker_neigh   = mean(share_broker_neigh_1sd, na.rm = TRUE),
    median_share_broker_neigh = median(share_broker_neigh_1sd, na.rm = TRUE)
  )

# Broker vs non-broker: own network position
broker_panel_1sd <- inventor_cent_tbl %>%
  filter(deg >= 2) %>%
  mutate(
    broker_group_1sd = if_else(is_broker_1sd, "broker_1sd", "non_broker")
  )

broker_panel_1sd %>%
  group_by(broker_group_1sd) %>%
  summarise(
    n_inventors = n(),
    mean_deg    = mean(deg, na.rm = TRUE),
    median_deg  = median(deg, na.rm = TRUE),
    mean_kcore  = mean(kcore, na.rm = TRUE),
    median_kcore = median(kcore, na.rm = TRUE),
    mean_constr = mean(constraint, na.rm = TRUE),
    sd_constr   = sd(constraint, na.rm = TRUE),
    .groups     = "drop"
  )

gc()

# ===================================================================== --
# 11. Internationality of brokers' ties ----
# ===================================================================== --

# Main country per inventor (country with most patents)
inv_loc <- inv_epo_link %>%
  filter(
    !is.na(person_id),
    !is.na(ctry_code)
  ) %>%
  distinct(person_id, appln_id, ctry_code)

inv_main_ctry <- inv_loc %>%
  count(person_id, ctry_code, name = "n_pat_ctry") %>%
  group_by(person_id) %>%
  slice_max(order_by = n_pat_ctry, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(main_ctry = ctry_code) %>%
  mutate(person_id = as.character(person_id))

# Edge list from igraph + countries on both endpoints
edges_df <- igraph::as_data_frame(coinv_net, what = "edges") %>%
  as_tibble() %>%
  mutate(
    from = as.character(from),
    to   = as.character(to)
  )

edges_ctry <- edges_df %>%
  left_join(inv_main_ctry, by = c("from" = "person_id")) %>%
  rename(ctry_from = main_ctry) %>%
  left_join(inv_main_ctry, by = c("to" = "person_id")) %>%
  rename(ctry_to = main_ctry) %>%
  filter(!is.na(ctry_from), !is.na(ctry_to)) %>%
  mutate(international = ctry_from != ctry_to)

# Inventor-level share of international ties
broker_flags <- inventor_cent %>%
  as_tibble() %>%
  mutate(person_id = as.character(person_id)) %>%
  select(person_id, is_broker_below_mean, is_broker_1sd)

edge_long <- edges_ctry %>%
  select(from, to, international) %>%
  tidyr::pivot_longer(
    cols      = c(from, to),
    names_to  = "endpoint",
    values_to = "person_id"
  ) %>%
  mutate(person_id = as.character(person_id))

inv_intl <- edge_long %>%
  group_by(person_id) %>%
  summarise(
    deg_edge   = n(), # degree based on edges
    intl_ties  = sum(international), # nr. of international ties
    share_intl = intl_ties / deg_edge, # share of international ties
    .groups    = "drop"
  ) %>%
  left_join(broker_flags, by = "person_id")

# Only inventors with deg >= 2
inv_intl_use <- inv_intl %>%
  filter(deg_edge >= 2)

# Descriptive stats for brokers vs non-brokers
intl_summary_bmean <- inv_intl_use %>%
  group_by(is_broker_below_mean) %>%
  summarise(
    n_inventors       = n(),
    mean_share_intl   = mean(share_intl, na.rm = TRUE),
    median_share_intl = median(share_intl, na.rm = TRUE),
    sd_share_intl     = sd(share_intl, na.rm = TRUE),
    .groups           = "drop"
  )

intl_summary_1sd <- inv_intl_use %>%
  group_by(is_broker_1sd) %>%
  summarise(
    n_inventors       = n(),
    mean_share_intl   = mean(share_intl, na.rm = TRUE),
    median_share_intl = median(share_intl, na.rm = TRUE),
    sd_share_intl     = sd(share_intl, na.rm = TRUE),
    .groups           = "drop"
  )

print(intl_summary_bmean)
print(intl_summary_1sd)

# t-test for "strong brokers" (1 SD below mean constraint)
t.test(
  share_intl ~ is_broker_1sd,
  data = inv_intl_use
)

# International share by strong broker status
ggplot(inv_intl_use, aes(x = is_broker_1sd, y = share_intl)) +
  geom_boxplot(
    outlier.alpha = 0.15,
    width         = 0.6,
    fill          = "grey80",
    colour        = "grey30"
  ) +
  scale_x_discrete(
    labels = c("FALSE" = "Non-broker", "TRUE" = "Broker (1 SD)")
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    name   = "Share of international co-inventor ties"
  ) +
  labs(
    x        = "Strong broker status (constraint \u2264 mean - 1 SD)",
    title    = "International orientation of brokers vs non-brokers",
    subtitle = "Distribution of the share of international co-inventor ties",
    caption  = "Source: own calculations based on OECD REGPAT data"
  ) +
  theme_minimal(base_family = "serif", base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 8)),
    axis.title    = element_text(size = 11),
    axis.text     = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.caption.position = "plot",
    plot.caption = element_text(
      hjust  = 1,
      face   = "italic",
      size   = 9,
      margin = margin(t = 10)
    ),
    plot.margin = margin(10, 15, 10, 10)
  )

gc()

# ===================================================================== --
# 12. Country–country co-inventor relations ----
# ===================================================================== --

# For each edge ctry_from -> ctry_to; and ctry_to -> ctry_from
ctry_edge_long <- edges_ctry %>%
  filter(!is.na(ctry_from), !is.na(ctry_to)) %>%
  transmute(ctry = ctry_from, partner_ctry = ctry_to) %>%
  bind_rows(
    edges_ctry %>%
      filter(!is.na(ctry_from), !is.na(ctry_to)) %>%
      transmute(ctry = ctry_to, partner_ctry = ctry_from)
  )

# Number of ties and shares per (country, partner country)
ctry_partner_stats <- ctry_edge_long %>%
  count(ctry, partner_ctry, name = "n_ties") %>%
  group_by(ctry) %>%
  mutate(
    tot_ties   = sum(n_ties),
    share_ties = n_ties / tot_ties
  ) %>%
  ungroup()

# For each country, sum = 1
ctry_partner_stats %>%
  group_by(ctry) %>%
  summarise(sum_share = sum(share_ties)) %>%
  arrange(desc(sum_share))

# How international each country is
ctry_intl_summary <- ctry_partner_stats %>%
  group_by(ctry) %>%
  summarise(
    share_domestic      = sum(share_ties[partner_ctry == ctry], na.rm = TRUE),
    share_international = 1 - share_domestic,
    n_partners          = dplyr::n_distinct(partner_ctry),
    .groups             = "drop"
  ) %>%
  arrange(desc(share_international))

print(ctry_intl_summary)

# Top partners per country (only foreign partners)
top_partners <- ctry_partner_stats %>%
  filter(partner_ctry != ctry) %>%
  group_by(ctry) %>%
  slice_max(order_by = share_ties, n = 5, with_ties = FALSE) %>%
  ungroup()

print(top_partners)

# Heatmap of country–partner shares
ggplot(ctry_partner_stats,
       aes(x = ctry, y = partner_ctry, fill = share_ties)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    x        = "Country",
    y        = "Partner country",
    fill     = "Share of ties",
    title    = "Distribution of co-inventor ties across country pairs",
    subtitle = "Row-normalised shares of all co-inventor ties",
    caption  = "Source: own calculations based on OECD REGPAT data"
  ) +
  theme_minimal(base_family = "serif", base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title  = element_text(size = 11),
    panel.grid  = element_blank(),
    legend.position = "right",
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 8)),
    plot.caption.position = "plot",
    plot.caption = element_text(
      hjust  = 1,
      face   = "italic",
      size   = 9,
      margin = margin(t = 10)
    ),
    plot.margin = margin(10, 15, 10, 10)
  )

ctry_partner_int <- ctry_partner_stats %>%
  filter(ctry != partner_ctry)

# Heatmap excluding domestic ties
ggplot(ctry_partner_int,
       aes(x = ctry, y = partner_ctry, fill = share_ties)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", , direction = -1) +
  labs(
    title    = "Distribution of international co-inventor ties",
    subtitle = "Row-normalised shares, domestic ties excluded",
    x        = "Country",
    y        = "Partner country",
    fill     = "Share of ties",
    caption  = "Source: own calculations based on OECD REGPAT data"
  ) +
  theme_minimal(base_family = "serif", base_size = 12) +
  theme(
    # white background
    panel.background = element_rect(colour = "white", fill = "white", ),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title  = element_text(size = 11),
    panel.grid  = element_blank(),
    legend.position = "right",
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 8)),
    plot.caption.position = "plot",
    plot.caption = element_text(
      hjust  = 1,
      face   = "italic",
      size   = 9,
      margin = margin(t = 10)
    ),
    plot.margin = margin(10, 15, 10, 10)
  )

# International share summary again
ctry_intl_summary <- ctry_partner_stats %>%
  group_by(ctry) %>%
  summarise(
    share_domestic      = sum(share_ties[partner_ctry == ctry], na.rm = TRUE),
    share_international = 1 - share_domestic,
    n_partners          = n_distinct(partner_ctry),
    .groups             = "drop"
  )
print(ctry_intl_summary)

gc()

# ===================================================================== --
# 13. Centrality and periphery ----
# ===================================================================== --

# Centralities & Burt's constraint (re-calculating per original structure)
inventor_cent <- tibble(
  person_id = as.numeric(V(coinv_net)$name),
  deg       = degree(coinv_net, mode = "all"),
  kcore     = coreness(coinv_net),
  constraint = constraint(coinv_net)
)

# Meaningful nodes: degree >= 2
inventor_cent_nz <- inventor_cent %>%
  filter(deg >= 2, !is.na(constraint))

# Mean & sd of constraint
mu_con <- mean(inventor_cent_nz$constraint, na.rm = TRUE)
sd_con <- sd(inventor_cent_nz$constraint, na.rm = TRUE)

# Z-score of constraint (lower = more broker-like)
inventor_cent <- inventor_cent %>%
  mutate(
    constraint_z = (constraint - mu_con) / sd_con
  )

# Define broker flags
inventor_cent <- inventor_cent %>%
  mutate(
    is_broker_below_mean = if_else(
      deg >= 2 & !is.na(constraint_z) & constraint_z < 0,
      TRUE,
      FALSE
    ),
    is_broker_1sd = if_else(
      deg >= 2 & !is.na(constraint_z) & constraint_z <= -1,
      TRUE,
      FALSE
    )
  )

# Filtering for usable inventors: deg >= 2, non-missing constraint
# Re-using inventor_panel from section 7 to ensure consistency if variables were added
inventor_panel_use <- inventor_panel %>%
  filter(deg >= 2, !is.na(constraint))

# Sanity checks
table(inventor_panel_use$is_broker_below_mean, useNA = "ifany")
table(inventor_panel_use$is_broker_1sd, useNA = "ifany")

# Top brokers and hubs
top_brokers_1sd <- inventor_panel_use %>%
  filter(is_broker_1sd) %>%
  arrange(constraint) %>%
  slice_head(n = 50)

top_hubs <- inventor_panel_use %>%
  arrange(desc(deg)) %>%
  slice_head(n = 50)

# Correlation table: degree, k-core, constraint
cor(
  inventor_panel_use %>% select(deg, kcore, constraint),
  use = "pairwise.complete.obs",
  method = "pearson"
)

# Core-periphery: k-core values distribution brokers and non-brokers
inventor_core <- inventor_panel_use %>%
  mutate(
    broker_group = case_when(
      is_broker_1sd ~ "Strong broker",
      is_broker_below_mean & !is_broker_1sd ~ "Weak broker",
      TRUE ~ "Non-broker"
    )
  )

# Proportion of strong and weak brokers
sum(inventor_core$broker_group == "Strong broker") / nrow(inventor_core)
sum(inventor_core$broker_group == "Weak broker") / nrow(inventor_core)
sum(inventor_core$broker_group == "Non-broker") / nrow(inventor_core)

# Plot: k-core index distribution
ggplot(inventor_core, aes(x = kcore, fill = broker_group)) +
  geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c(
    "Strong broker" = "red2",
    "Weak broker"   = "yellow",
    "Non-broker"    = "darkgrey"
  )) +
  labs(
    title = "Core periphery among brokers",
    caption = "Source: own calculations based on OECD REGPAT data",
    x = "k-core index",
    y = "Inventors",
    fill = "Broker type"
  ) +
  theme_minimal(base_family = "serif", base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 8)),
    axis.title    = element_text(size = 11),
    axis.text     = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.caption.position = "plot",
    plot.caption = element_text(
      hjust  = 1,
      face   = "italic",
      size   = 9,
      margin = margin(t = 10)
    ),
    plot.margin = margin(10, 15, 10, 10)
  )

# k-core based core-periphery categories
inventor_core <- inventor_core %>%
  mutate(
    core_periphery = case_when(
      kcore >= 5 ~ "Core",
      kcore < 3  ~ "Periphery",
      TRUE       ~ "Intermediate"
    )
  )

core_periphery_table <- inventor_core %>%
  group_by(broker_group, core_periphery) %>%
  summarise(n_inventors = n(), .groups = "drop") %>%
  group_by(broker_group) %>%
  mutate(
    share = n_inventors / sum(n_inventors)
  )

print(core_periphery_table)

# Order factors for plotting
core_periphery_table$broker_group <- factor(
  core_periphery_table$broker_group,
  levels = c("Non-broker", "Weak broker", "Strong broker")
)

# Plot
ggplot(core_periphery_table, aes(x = broker_group, y = share, fill = core_periphery)) +
  geom_col(alpha = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(
    "Intermediate" = "grey70",
    "Core"         = "red2",
    "Periphery"    = "green"
  )) +
  labs(
    title = "Core–Periphery distribution by broker type",
    caption = "Source: own calculations based on OECD REGPAT data",
    x = "Broker type",
    y = "Share of inventors",
    fill = "Network layer"
  ) +
  theme_minimal(base_family = "serif", base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 8)),
    axis.title    = element_text(size = 11),
    axis.text     = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.caption.position = "plot",
    plot.caption = element_text(
      hjust  = 1,
      face   = "italic",
      size   = 9,
      margin = margin(t = 10)
    ),
    plot.margin = margin(10, 15, 10, 10)
  )
gc()

# ===================================================================== --
# 14. Louvain community detection ----
# ===================================================================== --

# Louvain community detection on the co-inventor network
# (weight: number of shared patents)
comm_louvain <- cluster_louvain(
  coinv_net,
  weights = E(coinv_net)$weight
)

# Basic stats
print(modularity(comm_louvain))      # global modularity
print(length(comm_louvain))          # number of communities
print(sizes(comm_louvain)[1:10])     # sizes of the largest communities

# Membership vector (name = vertex name, value = community ID)
comm_vec <- membership(comm_louvain)

# Node-level table: inventor ~ community
inventor_comm <- tibble(
  person_id = as.numeric(names(comm_vec)),
  comm_id   = as.integer(comm_vec)
)

head(inventor_comm)


# ===================================================================== --
# 15. Number of neighbor communities per inventor ----
# ===================================================================== --

# Edge list: from, to (vertex names)
edges_tbl <- igraph::as_data_frame(coinv_net, what = "edges") %>%
  as_tibble()

# Attach community IDs to both endpoints using the named vector comm_vec
edges_tbl <- edges_tbl %>%
  mutate(
    comm_from = as.integer(comm_vec[as.character(from)]),
    comm_to   = as.integer(comm_vec[as.character(to)])
  )

# Node-level community info (needed later for counting "other communities")
node_comm_tbl <- inventor_comm %>%
  select(person_id, own_comm = comm_id)

# EGO-ORIENTED edge representation: each edge becomes two rows
# 1. ego = from, neighbor's community = comm_to
# 2. ego = to,   neighbor's community = comm_from

ego_comm <- bind_rows(
  edges_tbl %>%
    transmute(
      person_id  = as.numeric(from),
      neigh_comm = as.integer(comm_to)
    ),
  edges_tbl %>%
    transmute(
      person_id  = as.numeric(to),
      neigh_comm = as.integer(comm_from)
    )
)

# Join ego's own community
ego_comm <- ego_comm %>%
  left_join(node_comm_tbl, by = "person_id")

# Summary: number of neighbors, number of distinct neighbor communities,
# and among those, how many are "other" (different from ego's own community)
ego_comm_summary <- ego_comm %>%
  group_by(person_id) %>%
  summarise(
    n_neigh          = n(),                         # should match degree
    n_comm_neighbors = n_distinct(neigh_comm),      # all distinct neighbor communities
    n_other_comms    = n_distinct(neigh_comm[neigh_comm != own_comm]), # communities other than ego's own
    .groups          = "drop"
  )

head(ego_comm_summary)

# Note: isolated vertices (deg = 0) are not included here;
# we will assign 0 to them in the merge step.
gc()


# ===================================================================== --
# 16. Merge: number of communities + broker status + Uzzi ----
# ===================================================================== --

# Merge community and neighbor-community info into the main panel
inventor_panel <- inventor_panel %>%
  # Louvain community membership
  left_join(inventor_comm, by = "person_id") %>%
  # Number of neighbor communities
  left_join(ego_comm_summary, by = "person_id") %>%
  # Handle isolates (replace NA with 0)
  mutate(
    n_neigh          = if_else(is.na(n_neigh), 0L, as.integer(n_neigh)),
    n_comm_neighbors = if_else(is.na(n_comm_neighbors), 0L, as.integer(n_comm_neighbors)),
    n_other_comms    = if_else(is.na(n_other_comms), 0L, as.integer(n_other_comms))
  )

# Check: degree and n_neigh should mostly match
inventor_panel %>%
  summarise(
    cor_deg_n_neigh = cor(deg, n_neigh, use = "complete.obs")
  ) %>%
  print()

head(inventor_panel)

gc()

# ===================================================================== --
# 17. Brokers vs. non-brokers: how many communities do they connect to? ----
# ===================================================================== --

# Filter sample: only deg >= 2 and non-missing constraint
inventor_panel_use <- inventor_panel %>%
  filter(deg >= 2, !is.na(constraint))

# Descriptive statistics by strong broker (is_broker_1sd)
desc_comm_1sd <- inventor_panel_use %>%
  mutate(group_1sd = if_else(
    is_broker_1sd,
    "strong_broker_1sd",
    "non_strong_broker"
  )) %>%
  group_by(group_1sd) %>%
  summarise(
    n_inventors           = n(),
    share_of_sample       = n() / nrow(inventor_panel_use),
    
    # Network position
    mean_deg              = mean(deg, na.rm = TRUE),
    median_deg            = median(deg, na.rm = TRUE),
    
    # Community connections
    mean_n_comm_neighbors   = mean(n_comm_neighbors, na.rm = TRUE),
    median_n_comm_neighbors = median(n_comm_neighbors, na.rm = TRUE),
    sd_n_comm_neighbors     = sd(n_comm_neighbors, na.rm = TRUE),
    
    mean_n_other_comms    = mean(n_other_comms, na.rm = TRUE),
    median_n_other_comms  = median(n_other_comms, na.rm = TRUE),
    sd_n_other_comms      = sd(n_other_comms, na.rm = TRUE),
    
    .groups = "drop"
  )

print(desc_comm_1sd)

# T-tests:
# 1. Total neighbor communities
t_comm_total <- t.test(
  n_comm_neighbors ~ is_broker_1sd,
  data = inventor_panel_use
)
print(t_comm_total)

# 2. Other communities (excluding own)
t_comm_other <- t.test(
  n_other_comms ~ is_broker_1sd,
  data = inventor_panel_use
)
print(t_comm_other)

# Simple regressions with controls
# Total number of neighbor communities
lm_comm_total <- lm(
  n_comm_neighbors ~ is_broker_1sd + log1p(deg) + kcore,
  data = inventor_panel_use
)
summary(lm_comm_total)

# Number of other communities
lm_comm_other <- lm(
  n_other_comms ~ is_broker_1sd + log1p(deg) + kcore,
  data = inventor_panel_use
)
summary(lm_comm_other)
gc()


# ===================================================================== --
# 18. Visualizing distributions ----
# ===================================================================== --

# Boxplot: how many other communities an inventor connects to
ggplot(inventor_panel_use, aes(x = is_broker_1sd, y = n_other_comms)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_x_discrete(
    name = "Strong broker (constraint <= mean - 1 SD)"
  ) +
  scale_y_continuous(
    name = "Number of other communities (neighbor communities excluding own)"
  ) +
  theme_minimal(base_family = "serif", base_size = 12) +
  labs(
    title    = "Strong brokers and non-brokers: connections to communities",
    subtitle = "Co-inventor network, based on Louvain communities"
  )


# ===================================================================== --
# 19. Community size distribution ----
# ===================================================================== --

# Community sizes
comm_sizes <- inventor_comm %>%
  count(comm_id, name = "comm_size")

# Summary statistics
summary(comm_sizes$comm_size)
quantile(comm_sizes$comm_size, probs = c(0.5, 0.9, 0.99, 0.999))

# Histogram
ggplot(comm_sizes, aes(x = comm_size)) +
  geom_histogram(
    aes(fill = after_stat(count)),
    bins  = 40,
    color = "white",
    alpha = 0.9
  ) +
  scale_x_log10(
    breaks = c(1, 2, 5, 10, 50),
    labels = comma_format(),
    name   = "Community size (log10 scale, number of inventors)"
  ) +
  coord_cartesian(xlim = c(1, 50)) +
  scale_fill_viridis_c(option = "C", direction = -1, guide = "none") +
  labs(
    y       = "Number of communities",
    title   = "Size distribution of Louvain communities",
    caption = "Source: own calculations based on OECD REGPAT data"
  ) +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(
    plot.background   = element_rect(fill = "grey98", colour = NA),
    panel.background  = element_rect(fill = "grey96", colour = NA),
    panel.grid.major  = element_line(linewidth = 0.3, colour = "grey85"),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title        = element_text(size = 11),
    axis.text         = element_text(size = 10),
    plot.caption.position = "plot",
    plot.caption = element_text(
      hjust  = 1,
      face   = "italic",
      size   = 9,
      margin = margin(t = 10)
    ),
    plot.margin = margin(10, 15, 10, 10)
  )
gc()


# ===================================================================== --
# 20. Community size and share of strong brokers ----
# ===================================================================== --

# Community-level aggregation:
# Only inventors for whom brokerage is defined (deg >= 2 & !is.na(constraint))
comm_broker <- inventor_panel %>%
  filter(
    !is.na(comm_id),
    deg >= 2,
    !is.na(constraint)
  ) %>%
  group_by(comm_id) %>%
  summarise(
    comm_size       = n(),
    n_strong_broker = sum(is_broker_1sd, na.rm = TRUE),
    share_strong    = n_strong_broker / comm_size,
    .groups         = "drop"
  )

# Only sufficiently large communities
comm_broker_filt <- comm_broker %>%
  filter(comm_size >= 10)

# Quick stats
summary(comm_broker_filt$share_strong)

# Scatter plot: community size (log) vs. share of strong brokers
ggplot(comm_broker_filt, aes(x = comm_size, y = share_strong)) +
  geom_point(alpha = 0.15, size = 1.3, color = "grey40") +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se      = FALSE,
    span    = 0.5,
    size    = 1,
    color   = "steelblue4"
  ) +
  scale_x_log10(
    labels = comma_format(),
    name   = "Community size (log10, number of inventors – degree ≥ 2)"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    name   = "Share of strong brokers within the community"
  ) +
  labs(
    title    = "Concentration of strong brokers by community size",
    subtitle = "Louvain communities with at least 10 inventors",
    caption  = "Source: own calculations based on OECD REGPAT data"
  ) +
  theme_minimal(base_family = "serif", base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 8)),
    axis.title    = element_text(size = 11),
    axis.text     = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.caption.position = "plot",
    plot.caption = element_text(
      hjust  = 1,
      face   = "italic",
      size   = 9,
      margin = margin(t = 10)
    ),
    plot.margin = margin(10, 15, 10, 10)
  )

gc()

# ===================================================================== --
# 21. Country-level community analysis ----
# ===================================================================== --

# Inventor -> main country computed DIRECTLY from inv_epo_link
main_ctry_tbl <- inv_epo_link %>%
  filter(
    !is.na(person_id),
    !is.na(ctry_code)
  ) %>%
  distinct(person_id, appln_id, ctry_code) %>%
  count(person_id, ctry_code, name = "n_pat_ctry") %>%
  group_by(person_id) %>%
  slice_max(order_by = n_pat_ctry, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(main_ctry = ctry_code)

# Inventors + community + main country together
inventor_panel_ctry <- inventor_panel %>%
  left_join(main_ctry_tbl, by = "person_id") %>%
  filter(
    deg >= 2,
    !is.na(constraint),
    !is.na(comm_id),
    !is.na(main_ctry)
  )

# Community-level indicators: size, share of strong brokers
comm_broker_ctry <- inventor_panel_ctry %>%
  group_by(comm_id) %>%
  summarise(
    comm_size       = n(),
    n_strong_broker = sum(is_broker_1sd, na.rm = TRUE),
    share_strong    = n_strong_broker / comm_size,
    .groups         = "drop"
  )

# Assign dominant country to communities
comm_dom_ctry <- inventor_panel_ctry %>%
  group_by(comm_id, main_ctry) %>%
  summarise(
    n_inv_ctry = n(),
    .groups    = "drop"
  ) %>%
  group_by(comm_id) %>%
  slice_max(order_by = n_inv_ctry, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(dom_ctry = main_ctry)

# Merge: community size + broker share + dominant country
comm_stats_ctry <- comm_broker_ctry %>%
  left_join(comm_dom_ctry, by = "comm_id")

# Only communities with size >= 10
comm_stats_ctry_filt <- comm_stats_ctry %>%
  filter(comm_size >= 10)

# Select top countries (with the largest number of such communities)
top_ctrys <- comm_stats_ctry_filt %>%
  count(dom_ctry, sort = TRUE) %>%
  slice_head(n = 6) %>%
  pull(dom_ctry)

comm_stats_top <- comm_stats_ctry_filt %>%
  filter(dom_ctry %in% top_ctrys)

# Plot: community size vs. share of strong brokers, by country (facet)
ggplot(comm_stats_top, aes(x = comm_size, y = share_strong)) +
  geom_point(alpha = 0.15) +
  geom_smooth(
    method  = "loess",
    formula = y ~ x,
    se      = FALSE,
    span    = 0.6
  ) +
  scale_x_log10(
    labels = comma_format(),
    name   = "Community size (log10, number of inventors with deg ≥ 2)"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.5),
    name   = "Share of strong brokers within the community"
  ) +
  facet_wrap(~ dom_ctry) +
  theme_minimal(base_family = "serif", base_size = 12) +
  labs(
    title    = "Share of strong brokers in co-inventor communities, by country",
    subtitle = "Louvain communities, only communities with size ≥10; dominant country = country with most inventors"
  )

gc()

# ===================================================================== --
# 22. Average shortest path length (Broker analysis) ----
# ===================================================================== --

# Components, selecting the largest connected component (GCC)
comp <- components(coinv_net)
gcc_id  <- which.max(comp$csize)
gcc_vid <- which(comp$membership == gcc_id)

coinv_gcc <- induced_subgraph(coinv_net, vids = gcc_vid)

# Vertex table: which vertex belongs to which person_id, broker status, etc.
gcc_vertices_df <- tibble(
  vid_gcc   = seq_along(gcc_vid),
  vid_orig  = gcc_vid,
  person_id = as.numeric(V(coinv_net)$name)[gcc_vid]
) %>%
  left_join(
    inventor_panel %>% select(person_id, deg, constraint, is_broker_1sd),
    by = "person_id"
  )

# Brokerage definition: deg >= 2 & non-missing constraint
gcc_vertices_df <- gcc_vertices_df %>%
  mutate(
    broker_flag = case_when(
      deg >= 2 & !is.na(constraint) & is_broker_1sd  ~ "broker",
      deg >= 2 & !is.na(constraint) & !is_broker_1sd ~ "non_broker",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(broker_flag))

# Broker / non-broker vertex indices in the GCC
broker_vids <- gcc_vertices_df %>%
  filter(broker_flag == "broker") %>%
  pull(vid_gcc)

nonbroker_vids <- gcc_vertices_df %>%
  filter(broker_flag == "non_broker") %>%
  pull(vid_gcc)

print(length(broker_vids))
print(length(nonbroker_vids))

# Sampling to keep the number of pairs manageable
n_samp_broker    <- min(500, length(broker_vids))
n_samp_nonbroker <- min(500, length(nonbroker_vids))

samp_broker    <- sample(broker_vids,    size = n_samp_broker)
samp_nonbroker <- sample(nonbroker_vids, size = n_samp_nonbroker)

# Distance matrices (unweighted geodesic distance)

# 1. Broker–Broker
dist_bb_mat <- distances(
  coinv_gcc,
  v       = samp_broker,
  to      = samp_broker,
  mode    = "all",
  weights = NA
)
dist_bb <- dist_bb_mat[upper.tri(dist_bb_mat, diag = FALSE)]
dist_bb <- dist_bb[is.finite(dist_bb)]
avg_bb  <- mean(dist_bb)
med_bb  <- median(dist_bb)

# 2. Non-broker–Non-broker
dist_nn_mat <- distances(
  coinv_gcc,
  v       = samp_nonbroker,
  to      = samp_nonbroker,
  mode    = "all",
  weights = NA
)
dist_nn <- dist_nn_mat[upper.tri(dist_nn_mat, diag = FALSE)]
dist_nn <- dist_nn[is.finite(dist_nn)]
avg_nn  <- mean(dist_nn)
med_nn  <- median(dist_nn)

# 3. Broker–Non-broker
dist_bn_mat <- distances(
  coinv_gcc,
  v       = samp_broker,
  to      = samp_nonbroker,
  mode    = "all",
  weights = NA
)
dist_bn <- as.vector(dist_bn_mat)
dist_bn <- dist_bn[is.finite(dist_bn)]
avg_bn  <- mean(dist_bn)
med_bn  <- median(dist_bn)

# Summarize results
avg_path_summary <- tibble(
  pair_type = c("broker–broker", "broker–non-broker", "non-broker–non-broker"),
  mean_dist = c(avg_bb,           avg_bn,               avg_nn),
  med_dist  = c(med_bb,           med_bn,               med_nn),
  n_pairs   = c(length(dist_bb), length(dist_bn),      length(dist_nn))
)

print(avg_path_summary)

gc()

# ===================================================================== --
# 23. Sample size sensitivity: avg. path length for different n ----
# ===================================================================== --

# Function: for a given sample size, compute average paths
estimate_paths_for_n <- function(n_broker_samp, n_nonbroker_samp, seed = 20251202) {
  set.seed(seed)
  
  # Sample from indices (within GCC)
  n_broker_samp    <- min(n_broker_samp, length(broker_vids))
  n_nonbroker_samp <- min(n_nonbroker_samp, length(nonbroker_vids))
  
  samp_broker    <- sample(broker_vids,    size = n_broker_samp)
  samp_nonbroker <- sample(nonbroker_vids, size = n_nonbroker_samp)
  
  ## broker–broker
  dist_bb_mat <- distances(
    coinv_gcc,
    v       = samp_broker,
    to      = samp_broker,
    mode    = "all",
    weights = NA
  )
  dist_bb <- dist_bb_mat[upper.tri(dist_bb_mat, diag = FALSE)]
  dist_bb <- dist_bb[is.finite(dist_bb)]
  
  ## non-broker–non-broker
  dist_nn_mat <- distances(
    coinv_gcc,
    v       = samp_nonbroker,
    to      = samp_nonbroker,
    mode    = "all",
    weights = NA
  )
  dist_nn <- dist_nn_mat[upper.tri(dist_nn_mat, diag = FALSE)]
  dist_nn <- dist_nn[is.finite(dist_nn)]
  
  ## broker–non-broker
  dist_bn_mat <- distances(
    coinv_gcc,
    v       = samp_broker,
    to      = samp_nonbroker,
    mode    = "all",
    weights = NA
  )
  dist_bn <- as.vector(dist_bn_mat)
  dist_bn <- dist_bn[is.finite(dist_bn)]
  
  tibble(
    n_broker_samp    = n_broker_samp,
    n_nonbroker_samp = n_nonbroker_samp,
    pair_type        = c("broker–broker", "broker–non-broker", "non-broker–non-broker"),
    mean_dist        = c(mean(dist_bb),   mean(dist_bn),        mean(dist_nn)),
    med_dist         = c(median(dist_bb), median(dist_bn),      median(dist_nn)),
    n_pairs          = c(length(dist_bb), length(dist_bn),      length(dist_nn))
  )
}

# Chosen sample sizes (symmetric for brokers/non-brokers)
n_grid <- c(200, 400, 600, 800, 1000)

# Run for each n
avg_path_by_n <- purrr::map_dfr(
  n_grid,
  ~ estimate_paths_for_n(n_broker_samp = .x, n_nonbroker_samp = .x)
)

print(avg_path_by_n)

gc()

# ===================================================================== --
# 24. Targeted vs random broker removal in the GCC ----
# ===================================================================== --

# GCC and broker list
# Components, largest connected component
comp <- components(coinv_net)
gcc_id  <- which.max(comp$csize)
gcc_vid <- which(comp$membership == gcc_id)

coinv_gcc <- induced_subgraph(coinv_net, vids = gcc_vid)

# Vertex table: GCC index + person_id + brokerage info
gcc_vertices_df <- tibble(
  vid_gcc   = seq_along(gcc_vid),
  vid_orig  = gcc_vid,
  person_id = as.numeric(V(coinv_net)$name)[gcc_vid]
) %>%
  left_join(
    inventor_panel %>%
      select(person_id, deg, constraint, is_broker_1sd),
    by = "person_id"
  ) %>%
  mutate(
    broker_flag = case_when(
      deg >= 2 & !is.na(constraint) & is_broker_1sd  ~ "broker",
      deg >= 2 & !is.na(constraint) & !is_broker_1sd ~ "non_broker",
      TRUE ~ NA_character_
    )
  )

# Only vertices with defined brokerage
gcc_vertices_use <- gcc_vertices_df %>%
  filter(!is.na(broker_flag))

# Degree in the GCC
deg_gcc <- degree(coinv_gcc, mode = "all")

gcc_vertices_use <- gcc_vertices_use %>%
  mutate(deg_gcc = deg_gcc[vid_gcc])

# Broker and non-broker indices
broker_vids <- gcc_vertices_use %>%
  filter(broker_flag == "broker") %>%
  pull(vid_gcc)

nonbroker_vids <- gcc_vertices_use %>%
  filter(broker_flag == "non_broker") %>%
  pull(vid_gcc)

print(length(broker_vids))
print(length(nonbroker_vids))

# Selecting the largest brokers (targeted removal list)
top_brokers <- gcc_vertices_use %>%
  filter(broker_flag == "broker") %>%
  arrange(desc(deg_gcc))

top_broker_vids <- top_brokers$vid_gcc  # sorted in decreasing order of degree

# Function: estimated average path length in the (remaining) GCC
compute_apl_gcc <- function(g, n_sample = 500, seed = 1L) {
  set.seed(seed)
  
  # Always work with the current largest component
  comp <- components(g)
  gcc_id  <- which.max(comp$csize)
  gcc_vid <- which(comp$membership == gcc_id)
  g_gcc   <- induced_subgraph(g, vids = gcc_vid)
  
  n <- vcount(g_gcc)
  if (n < 2) return(NA_real_)
  
  n_samp <- min(n_sample, n)
  vids   <- sample(V(g_gcc), size = n_samp)
  
  # Random vertex pairs (n_samp x n_samp matrix, upper triangle)
  dmat <- distances(
    g_gcc,
    v       = vids,
    to      = vids,
    mode    = "all",
    weights = NA
  )
  
  d <- dmat[upper.tri(dmat, diag = FALSE)]
  d <- d[is.finite(d) & d > 0]
  if (!length(d)) return(NA_real_)
  
  mean(d)
}

# Parameters: how many brokers to remove, how many random repetitions
k_vec <- c(100, 500, 1000, 5000, 10000)
# Ensure k does not exceed the broker population
k_vec <- k_vec[k_vec < length(broker_vids)]

n_rep  <- 10   # Number of random experiment repetitions
n_samp <- 500  # Vertex sample size for path length estimation

# Baseline APL (no removal)
apl_baseline <- compute_apl_gcc(coinv_gcc, n_sample = n_samp, seed = 100)

# Running targeted vs random removal
results_list <- list()

for (k in k_vec) {
  message("Running: k = ", k)
  
  # 1. Targeted: top-k brokers by degree
  g_target <- delete_vertices(coinv_gcc, top_broker_vids[1:k])
  apl_target <- compute_apl_gcc(g_target, n_sample = n_samp, seed = 100 + k)
  
  # 2. Random: remove k brokers at random, n_rep repetitions
  apl_random_vec <- numeric(n_rep)
  
  for (r in seq_len(n_rep)) {
    rand_brokers <- sample(broker_vids, size = k, replace = FALSE)
    g_rand       <- delete_vertices(coinv_gcc, rand_brokers)
    apl_random_vec[r] <- compute_apl_gcc(g_rand, n_sample = n_samp, seed = 200 + k * 10 + r)
  }
  
  results_list[[as.character(k)]] <- tibble(
    k_removed        = k,
    removal_type     = c("targeted", "random_mean"),
    apl_estimate     = c(apl_target, mean(apl_random_vec)),
    random_sd        = c(NA_real_,   sd(apl_random_vec)),
    random_rep_count = c(NA_integer_, n_rep)
  )
}

results_df <- bind_rows(
  tibble(
    k_removed        = 0L,
    removal_type     = "baseline",
    apl_estimate     = apl_baseline,
    random_sd        = NA_real_,
    random_rep_count = NA_integer_
  ),
  bind_rows(results_list)
)

print(results_df)

# Visualization: change in APL under broker removal
results_plot <- results_df %>%
  # Only the two types: targeted and random
  filter(removal_type %in% c("targeted", "random_mean")) %>%
  mutate(
    removal_type = factor(
      removal_type,
      levels = c("targeted", "random_mean"),
      labels = c(
        "Targeted: top brokers",
        "Random: random brokers"
      )
    )
  )

ggplot(results_plot, aes(x = k_removed, y = apl_estimate, colour = removal_type)) +
  geom_point(size = 2) +
  geom_line() +
  scale_x_continuous(
    breaks       = c(100, 500, 1000, 5000, 10000),
    minor_breaks = c(200, 300, 400, 2000, 3000, 4000),
    labels       = comma_format(),
    name         = "Number of removed brokers (in GCC)"
  ) +
  scale_y_continuous(
    name = "Average path length"
  ) +
  scale_colour_manual(
    values = c("red", "steelblue"),
    name   = "Removal type"
  ) +
  labs(
    title    = "Targeted vs random removal of brokers in the co-inventor GCC",
    subtitle = "Top-degree brokers vs random brokers – impact on average path length",
    caption  = "Source: own calculations based on OECD REGPAT data"
  ) +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 8)),
    axis.title    = element_text(size = 11),
    axis.text     = element_text(size = 10),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.caption.position = "plot",
    plot.caption = element_text(
      hjust  = 1,
      face   = "italic",
      size   = 9,
      margin = margin(t = 10)
    ),
    plot.margin = margin(10, 15, 10, 10)
  )

gc()

# ===================================================================== --
# 25. Region-level edges (EPO only) & TL3→TL2 mapping ----
# ===================================================================== --

# We REUSE:
#   regions_eu, eu_ctry_codes, eu_nuts3_codes  (defined at the top)
#   inv_epo_eu                                  (EPO EU inventor–region data)
# No re-loading of REGPAT, no PCT here.

regions_map <- regions_eu %>%
  select(
    ctry_code,
    reg_code,
    reg_label,
    up_reg_code,
    up_reg_label
  )

# TL3 -> TL2 mapping
map_tl3_to_tl2 <- regions_map %>%
  distinct(
    reg_code,
    nuts2      = up_reg_code,
    nuts2_name = up_reg_label,
    ctry_code
  )

# TL3 region–region edges from EPO inventor-reg data (EU-only)
dt_unique <- inv_epo_eu %>%
  filter(
    ctry_code %in% eu_ctry_codes,
    !is.na(reg_code),
    reg_code %in% eu_nuts3_codes
  ) %>%
  distinct(appln_id, reg_code)

edges_tl3 <- dt_unique %>%
  inner_join(
    dt_unique,
    by = "appln_id",
    suffix = c("_from", "_to")
  ) %>%
  filter(reg_code_from < reg_code_to) %>%   # unordered, no self-loops
  count(from = reg_code_from, to = reg_code_to, name = "weight")

# Aggregate TL3 edges to TL2 / NUTS2 level (EPO-only)
edges_tl2 <- edges_tl3 %>%
  left_join(map_tl3_to_tl2, by = c("from" = "reg_code")) %>%
  rename(from_tl2 = nuts2) %>%
  left_join(map_tl3_to_tl2, by = c("to" = "reg_code")) %>%
  rename(to_tl2 = nuts2) %>%
  filter(!is.na(from_tl2), !is.na(to_tl2)) %>%
  filter(from_tl2 != to_tl2) %>%           # drop within-region TL2 edges
  mutate(
    new_from = pmin(from_tl2, to_tl2),
    new_to   = pmax(from_tl2, to_tl2)
  ) %>%
  group_by(from = new_from, to = new_to) %>%
  summarise(weight = sum(weight), .groups = "drop")


# ===================================================================== --
# 26. NUTS2 Network Aggregation and Mapping ----
# ===================================================================== --

# TL2 / NUTS2 node table + strength (weighted degree)
nuts2_ids_used <- unique(c(edges_tl2$from, edges_tl2$to))

nodes_tl2 <- regions_eu %>%
  filter(up_reg_code %in% nuts2_ids_used) %>%
  transmute(
    nuts_id   = up_reg_code,
    ctry_code,
    reg_label = up_reg_label
  ) %>%
  distinct()

# Node strength
node_strength <- bind_rows(
  edges_tl2 %>% select(nuts_id = from, weight),
  edges_tl2 %>% select(nuts_id = to, weight)
) %>%
  group_by(nuts_id) %>%
  summarise(strength = sum(weight), .groups = "drop")

nodes_tl2 <- nodes_tl2 %>%
  left_join(node_strength, by = "nuts_id") %>%
  mutate(strength = replace_na(strength, 0))

# Get NUTS2 shapes and centroids
nuts2_sf <- gisco_get_nuts(
  nuts_level   = 2,
  year         = 2010,
  resolution   = "20",
  update_cache = FALSE
) %>%
  filter(CNTR_CODE %in% eu_ctry_codes) %>%
  rename(nuts_id = NUTS_ID)

# Full NUTS2 layer with strength
nuts2_full_sf <- nuts2_sf %>%
  left_join(nodes_tl2 %>% select(nuts_id, strength), by = "nuts_id") %>%
  mutate(strength = replace_na(strength, 0))

centroids_full <- st_point_on_surface(nuts2_full_sf)

# Build LINESTRING edges between NUTS2 centroids
centroids_network <- centroids_full %>%
  filter(nuts_id %in% nuts2_ids_used) %>%
  select(nuts_id, geometry)

edges_geo <- edges_tl2 %>%
  left_join(centroids_network, by = c("from" = "nuts_id")) %>%
  rename(geom_from = geometry) %>%
  left_join(centroids_network, by = c("to" = "nuts_id")) %>%
  rename(geom_to = geometry) %>%
  filter(!st_is_empty(geom_from), !st_is_empty(geom_to))

# Create spatial lines
edge_lines <- map(seq_len(nrow(edges_geo)), function(i) {
  coords_from <- st_coordinates(edges_geo$geom_from[i])
  coords_to   <- st_coordinates(edges_geo$geom_to[i])
  st_linestring(rbind(coords_from, coords_to))
})

edges_sf <- st_sf(
  from     = edges_geo$from,
  to       = edges_geo$to,
  weight   = edges_geo$weight,
  geometry = st_sfc(edge_lines, crs = st_crs(centroids_full))
)

# Visualization preparation (LAEA Europe)
crs_laea       <- 3035
countries_laea <- gisco_get_countries(region = "Europe") %>% st_transform(crs_laea)
nodes_full_laea <- st_transform(centroids_full, crs_laea)
edges_laea      <- st_transform(edges_sf, crs_laea)

# Filter for mainland Europe for plotting (wider to include e.g. Ukraine)
nodes_full_ll <- st_transform(centroids_full, 4326)
coords        <- st_coordinates(nodes_full_ll)
in_europe     <- coords[, "X"] >= -10 & coords[, "X"] <= 40 &
  coords[, "Y"] >= 35  & coords[, "Y"] <= 71

nodes_main_laea <- nodes_full_laea[in_europe, ]
main_ids        <- nodes_main_laea$nuts_id

edges_main_laea <- edges_laea %>%
  filter(from %in% main_ids, to %in% main_ids)

# Top 10% ties
q_strong  <- quantile(edges_main_laea$weight, 0.90, na.rm = TRUE)
edges_vis <- edges_main_laea %>% filter(weight >= q_strong)

# Map (only plot nodes with strength > 0 in colour)
nodes_main_laea_df <- nodes_main_laea %>%
  mutate(
    has_patents = strength > 0
  )

ggplot() +
  geom_sf(
    data  = countries_laea,
    fill  = "#f7f7f7",
    color = "#d0d0d0",
    linewidth = 0.2
  ) +
  geom_sf(
    data  = edges_vis,
    aes(linewidth = weight),
    colour = "grey40",
    alpha   = 0.45
  ) +
  scale_linewidth_continuous(range = c(0.2, 0.9), guide = "none") +
  geom_sf(
    data = nodes_main_laea_df %>% filter(!has_patents),
    color = "grey70",
    alpha = 0.4,
    size  = 1.1
  ) +
  geom_sf(
    data = nodes_main_laea_df %>% filter(has_patents),
    aes(color = log1p(strength)),
    size  = 1.6
  ) +
  scale_color_viridis_c(
    name   = "Node strength\nlog(1 + sum of ties)",
    option = "C"
  ) +
  coord_sf(
    xlim   = c(st_bbox(nodes_main_laea)["xmin"], st_bbox(nodes_main_laea)["xmax"]),
    ylim   = c(st_bbox(nodes_main_laea)["ymin"], st_bbox(nodes_main_laea)["ymax"]),
    expand = FALSE
  ) +
  theme_void(base_family = "serif") +
  theme(
    legend.position  = "right",
    panel.background = element_rect(fill = "#f8f9fa", colour = NA),
    plot.title       = element_text(face = "bold", size = 16, hjust = 0.5)
  ) +
  ggtitle("European Co-Inventor Network (NUTS2 level)")

gc()

# ===================================================================== --
# 27. Broker Structural Roles: Data Prep ----
# ===================================================================== --



# EPO inventor–region + IPC information (REUSED objects)
# inv_epo already loaded earlier; we keep inv_epo_eu as constructed before

inv_epo_eu <- inv_epo_eu %>%
  select(
    appln_id, person_id, reg_code, ctry_code, reg_share,
    inv_share, app_year, ipc4
  )

inv_epo_link <- unique(inv_epo_eu)

# Main country per inventor (EPO-based)
inv_main_ctry <- inv_epo_link %>%
  distinct(person_id, appln_id, ctry_code) %>%
  count(person_id, ctry_code) %>%
  group_by(person_id) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  transmute(person_id = as.numeric(person_id), main_ctry = ctry_code) %>%
  ungroup()

# Load base data and merge
if (!exists("inventor_panel")) stop("`inventor_panel` not found.")

dat <- inventor_panel %>%
  as_tibble() %>%
  left_join(inv_main_ctry, by = "person_id")

# Filter to meaningful brokerage cases
inventor_panel_geo <- dat %>%
  filter(deg >= 2, !is.na(constraint), !is.na(main_ctry))

print(nrow(inventor_panel_geo))

# Helper theme
theme_brokers <- function() {
  theme_minimal(base_family = "serif", base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 8)),
      axis.title    = element_text(size = 11),
      axis.text     = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 1, face = "italic", size = 9, margin = margin(t = 10))
    )
}

# Country-level broker vs non-broker profiles
ctry_broker_profiles <- inventor_panel_geo %>%
  mutate(group = if_else(is_broker_1sd, "strong_broker_1sd", "non_broker")) %>%
  group_by(main_ctry, group) %>%
  summarise(
    n_inventors       = n(),
    mean_deg          = mean(deg, na.rm = TRUE),
    mean_kcore        = mean(kcore, na.rm = TRUE),
    mean_constraint   = mean(constraint, na.rm = TRUE),
    mean_n_pat_tot    = mean(n_pat_tot, na.rm = TRUE),
    mean_share_atyp   = mean(share_atyp_tot, na.rm = TRUE),
    .groups = "drop"
  )

ctry_strong_brokers <- ctry_broker_profiles %>%
  filter(group == "strong_broker_1sd") %>%
  arrange(desc(n_inventors))

gc()

# ===================================================================== --
# 29. Feature-based Role Discovery (K-Means) ----
# ===================================================================== --

# Identify strong brokers and Calculate Local Features
target_brokers <- inventor_panel_geo %>%
  filter(is_broker_1sd == TRUE)

broker_ids_char <- as.character(target_brokers$person_id)

# Helper: local features
get_node_features <- function(g, node_names) {
  vids <- match(node_names, V(g)$name)
  vids <- vids[!is.na(vids)]
  if (length(vids) == 0) stop("No matching vertices found.")
  
  tibble(
    person_id = as.numeric(V(g)$name[vids]),
    clust_co  = transitivity(g, type = "local", vids = vids),
    neigh_deg = knn(g)$knn[vids]
  )
}

topo_features <- get_node_features(coinv_net, broker_ids_char)

broker_features_df <- target_brokers %>%
  select(person_id, main_ctry, deg, constraint, kcore) %>%
  inner_join(topo_features, by = "person_id") %>%
  mutate(
    clust_co  = replace_na(clust_co, 0),
    neigh_deg = replace_na(neigh_deg, 0)
  )

# Cluster Brokers into Roles (K-Means)
features_mat <- broker_features_df %>%
  select(deg, constraint, kcore, clust_co, neigh_deg) %>%
  scale()

set.seed(2025)
km_res <- kmeans(features_mat, centers = 4, nstart = 25)
broker_features_df$role_id <- factor(km_res$cluster)

# Role centroids for interpretation
role_centers <- as.data.frame(km_res$centers) %>%
  rownames_to_column("role_id") %>%
  pivot_longer(-role_id, names_to = "feature", values_to = "z_score")

ggplot(role_centers, aes(x = feature, y = z_score, fill = role_id)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ role_id, nrow = 1) +
  scale_fill_brewer(palette = "Set2", name = "Role ID") +
  labs(
    title = "Broker roles based on local network features",
    y     = "Standardized feature value (z-score)"
  ) +
  theme_brokers() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Overall role distribution
role_summary_overall <- broker_features_df %>%
  group_by(role_id) %>%
  summarise(
    N               = n(),
    Share           = n() / nrow(broker_features_df),
    Mean_Degree     = mean(deg),
    Mean_Constraint = mean(constraint),
    Mean_KCore      = mean(kcore),
    .groups = "drop"
  )

print(role_summary_overall)




gc()

# ===================================================================== --
# 30. Visual Validation: Ego Networks (1-step & 2-step) ----
# ===================================================================== --

# Determine descriptive labels based on feature means
profile_summary <- broker_features_df %>%
  group_by(role_id) %>%
  summarise(
    m_deg    = mean(deg, na.rm = TRUE),
    m_clust  = mean(clust_co, na.rm = TRUE),
    m_constr = mean(constraint, na.rm = TRUE),
    .groups  = "drop"
  )

# Assign types logic
id_1 <- profile_summary %>% arrange(desc(m_clust)) %>% slice(1) %>% pull(role_id)
id_2 <- profile_summary %>% filter(role_id != id_1) %>% arrange(desc(m_deg)) %>% slice(1) %>% pull(role_id)
id_3 <- profile_summary %>% filter(!role_id %in% c(id_1, id_2)) %>% arrange(m_constr) %>% slice(1) %>% pull(role_id)
id_4 <- setdiff(profile_summary$role_id, c(id_1, id_2, id_3))

new_labels <- setNames(
  c("Type 1: High clustering", "Type 2: High degree", "Type 3: Low constraint", "Type 4: Other"),
  c(id_1, id_2, id_3, id_4)
)

broker_features_df <- broker_features_df %>%
  mutate(role_label = new_labels[as.character(role_id)])

# Ego plotting function (kept as your "advanced" visual check)
plot_ego_1_2 <- function(role_name, g_full = coinv_net, seed = 2025) {
  set.seed(seed)
  
  # 1. Pick an exemplar inventor
  ex_row <- broker_features_df %>%
    filter(!is.na(role_label)) %>%
    filter(role_label == role_name | grepl(role_name, role_label, fixed = TRUE)) %>%
    arrange(desc(deg)) %>%
    slice(1)
  
  if (nrow(ex_row) == 0L) stop("No brokers found for role: ", role_name)
  
  ex_id       <- ex_row$person_id[1]
  ex_deg      <- ex_row$deg[1]
  ex_kcore    <- ex_row$kcore[1]
  ex_constr   <- ex_row$constraint[1]
  ex_clust    <- ex_row$clust_co[1]
  ex_neighdeg <- ex_row$neigh_deg[1]
  
  # 2. Match exemplar to vertex
  v_idx <- which(V(g_full)$name == as.character(ex_id))
  if (length(v_idx) == 0L) stop("Exemplar ID not found in graph.")
  
  # 3. Create ego graphs
  ego1 <- make_ego_graph(g_full, order = 1, nodes = v_idx)[[1]]
  ego2 <- make_ego_graph(g_full, order = 2, nodes = v_idx)[[1]]
  
  # 4. Layer definition
  dist2 <- distances(
    ego2,
    v = which(V(ego2)$name == V(g_full)$name[v_idx])
  )[1, ]
  
  V(ego2)$layer <- case_when(
    dist2 == 1 ~ "1-step neighbours",
    dist2 == 2 ~ "2-step neighbours",
    TRUE ~ "ego"
  )
  
  # 5. Plotting
  op <- par(mfrow = c(1, 2), mar = c(1, 1, 4, 1), bg = "white")
  
  # Left: 1-step
  lay1   <- layout_with_fr(ego1)
  v_col1 <- if_else(V(ego1)$name == V(g_full)$name[v_idx], "#d73027", "#4575b4")
  v_size1 <- if_else(V(ego1)$name == V(g_full)$name[v_idx], 9, 4.5)
  
  plot(
    ego1, layout = lay1, vertex.size = v_size1, vertex.label = NA,
    vertex.color = v_col1, vertex.frame.color = NA, edge.color = "grey80",
    edge.width = 0.6, edge.curved = 0.1,
    main = paste0(role_name, " – 1-step ego\n"),
    sub  = paste0("Degree: ", ex_deg, "   K-core: ", ex_kcore,
                  "   Constraint: ", round(ex_constr, 3))
  )
  
  # Right: 2-step
  lay2 <- layout_with_fr(ego2)
  v_col2 <- case_when(
    V(ego2)$name == V(g_full)$name[v_idx] ~ "#d73027",
    V(ego2)$layer == "1-step neighbours" ~ "#4575b4",
    TRUE ~ "#91bfdb"
  )
  v_size2 <- case_when(
    V(ego2)$name == V(g_full)$name[v_idx] ~ 9,
    V(ego2)$layer == "1-step neighbours" ~ 4.5,
    TRUE ~ 3
  )
  
  plot(
    ego2, layout = lay2, vertex.size = v_size2, vertex.label = NA,
    vertex.color = v_col2, vertex.frame.color = NA, edge.color = "grey85",
    edge.width = 0.5, edge.curved = 0.1,
    main = paste0(role_name, " – 2-step ego\n"),
    sub  = paste0("Local clustering: ", round(ex_clust, 3),
                  "   Avg neighbor degree: ", round(ex_neighdeg, 1))
  )
  
  par(op)
}

# Generate plots (advanced visual examples only)
plot_ego_1_2("Type 1: High clustering")
plot_ego_1_2("Type 2: High degree")
plot_ego_1_2("Type 3: Low constraint")
plot_ego_1_2("Type 4: Other")

gc()

# ===================================================================== --
# 31. Descriptive Statistics for Broker Roles ----
# ===================================================================== --

library(kableExtra)

role_descriptives <- broker_features_df %>%
  group_by(role_label) %>%
  summarise(
    N_brokers         = n(),
    Share             = n() / nrow(broker_features_df),
    Mean_Degree       = mean(deg, na.rm = TRUE),
    Median_Degree     = median(deg, na.rm = TRUE),
    Mean_KCore        = mean(kcore, na.rm = TRUE),
    Median_KCore      = median(kcore, na.rm = TRUE),
    Mean_Constraint   = mean(constraint, na.rm = TRUE),
    Median_Constraint = median(constraint, na.rm = TRUE),
    Mean_Clustering   = mean(clust_co, na.rm = TRUE),
    Median_Clustering = median(clust_co, na.rm = TRUE),
    Mean_NeighDeg     = mean(neigh_deg, na.rm = TRUE),
    Median_NeighDeg   = median(neigh_deg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(role_label)

print(role_descriptives)

role_table <- role_descriptives %>%
  mutate(
    Share             = scales::percent(Share, accuracy = 0.1),
    Mean_Degree       = round(Mean_Degree, 1),
    Median_Degree     = round(Median_Degree, 1),
    Mean_KCore        = round(Mean_KCore, 2),
    Median_KCore      = round(Median_KCore, 2),
    Mean_Constraint   = round(Mean_Constraint, 3),
    Median_Constraint = round(Median_Constraint, 3),
    Mean_Clustering   = round(Mean_Clustering, 3),
    Median_Clustering = round(Median_Clustering, 3),
    Mean_NeighDeg     = round(Mean_NeighDeg, 1),
    Median_NeighDeg   = round(Median_NeighDeg, 1)
  ) %>%
  select(
    `Broker type`        = role_label,
    `N`                  = N_brokers,
    `Share`              = Share,
    `Mean degree`        = Mean_Degree,
    `Median degree`      = Median_Degree,
    `Mean k-core`        = Mean_KCore,
    `Median k-core`      = Median_KCore,
    `Mean constraint`    = Mean_Constraint,
    `Median constraint`  = Median_Constraint,
    `Mean clustering`    = Mean_Clustering,
    `Median clustering`  = Median_Clustering,
    `Mean neigh. deg.`   = Mean_NeighDeg,
    `Median neigh. deg.` = Median_NeighDeg
  )

role_table %>%
  kable("html", caption = "Descriptive Statistics of Broker Types",
        align = "lrrrrrrrrrrrr") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, bold = TRUE) %>%
  scroll_box(width = "100%", height = "400px")


# ===================================================================== --
# 32. Geographic Map: Dominant Broker Roles (NUTS2) ----
# ===================================================================== --

# 1. Link brokers to their NUTS2 region
map_nuts3_nuts2 <- regions_eu %>%
  distinct(reg_code, nuts2 = up_reg_code)

broker_locs <- broker_features_df %>%
  select(person_id, role_label) %>%
  inner_join(
    inv_epo_link %>% select(person_id, reg_code),
    by = "person_id",
    relationship = "many-to-many"
  ) %>%
  inner_join(map_nuts3_nuts2, by = "reg_code") %>%
  distinct(person_id, nuts2, role_label)

# 2. Regional role profiles + Location Quotient
reg_role_counts <- broker_locs %>%
  count(nuts2, role_label, name = "n_role") %>%
  group_by(nuts2) %>%
  mutate(n_total = sum(n_role)) %>%
  ungroup() %>%
  filter(n_total >= 50)

eu_baseline <- broker_features_df %>%
  count(role_label, name = "n_eu") %>%
  mutate(share_eu = n_eu / sum(n_eu)) %>%
  select(role_label, share_eu)

reg_specialization <- reg_role_counts %>%
  left_join(eu_baseline, by = "role_label") %>%
  mutate(
    share_reg = n_role / n_total,
    LQ        = share_reg / share_eu
  )

dominant_role_map <- reg_specialization %>%
  group_by(nuts2) %>%
  arrange(desc(LQ)) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(nuts2, dominant_role = role_label, max_lq = LQ, n_total)

# Prepare Map Data
nuts2_map_ready <- nuts2_sf %>%
  st_transform(3035) %>%
  left_join(dominant_role_map, by = c("nuts_id" = "nuts2"))

map_plot_data <- nuts2_map_ready %>%
  filter(!is.na(dominant_role)) %>%
  st_crop(xmin = 2500000, xmax = 7000000,
          ymin = 1500000, ymax = 5500000)

ggplot() +
  geom_sf(data = countries_laea, fill = "grey95", color = "white", linewidth = 0.2) +
  geom_sf(data = map_plot_data, aes(fill = dominant_role), color = "grey80", linewidth = 0.1) +
  scale_fill_brewer(palette = "Set2", name = "Dominant broker role") +
  labs(
    title    = "Dominant broker roles by NUTS2 region",
    subtitle = "Role with the highest location quotient (regions with ≥ 50 strong brokers)",
    caption  = "Location quotient relative to EU-wide role shares."
  ) +
  theme_void(base_family = "serif") +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 16),
    plot.subtitle   = element_text(color = "grey40")
  )

gc()

# ===================================================================== --
# 33. Clustering Countries by Broker Mix ----
# ===================================================================== --

# Select countries with sufficient data
valid_ctry <- broker_features_df %>%
  count(main_ctry) %>%
  filter(n > 100) %>%
  pull(main_ctry)

# Create matrix of shares
ctry_mix <- broker_features_df %>%
  filter(main_ctry %in% valid_ctry) %>%
  count(main_ctry, role_label) %>%
  group_by(main_ctry) %>%
  mutate(share = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = role_label, values_from = share, values_fill = 0) %>%
  column_to_rownames("main_ctry")

# Clustering (Ward, k=3)
hc      <- hclust(dist(ctry_mix), method = "ward.D2")
ctry_grps <- cutree(hc, k = 3)

cluster_df <- tibble(
  main_ctry  = names(ctry_grps),
  cluster_id = factor(ctry_grps)
) %>%
  mutate(cluster_name = paste0("Group ", as.numeric(cluster_id)))

# Stacked bar: broker mix by country group (advanced visualization)
p_bars_final <- broker_features_df %>%
  filter(main_ctry %in% cluster_df$main_ctry) %>%
  left_join(cluster_df, by = "main_ctry") %>%
  ggplot(aes(x = reorder(main_ctry, as.numeric(cluster_id)), fill = role_label)) +
  geom_bar(position = "fill", width = 0.8) +
  facet_grid(~ cluster_name, scales = "free_x", space = "free") +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set2", name = "Broker type") +
  labs(
    title    = "Broker role composition by country group",
    subtitle = "Countries grouped by Ward clustering (k = 3) of broker-type shares",
    x        = "",
    y        = "Share of strong brokers"
  ) +
  theme_brokers() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

print(p_bars_final)

# Map: NUTS2 regions coloured by country cluster
map_data <- nuts2_sf %>%
  st_transform(3035) %>%
  left_join(cluster_df, by = c("CNTR_CODE" = "main_ctry")) %>%
  filter(!is.na(cluster_name)) %>%
  st_crop(xmin = 2600000, xmax = 6200000, ymin = 1450000, ymax = 5400000)

bg_map <- gisco_get_countries(region = "Europe") %>% st_transform(3035)

p_map_final <- ggplot() +
  geom_sf(data = bg_map, fill = "#f0f0f0", color = "white", linewidth = 0.2) +
  geom_sf(data = map_data, aes(fill = cluster_name), color = "white", linewidth = 0.05) +
  scale_fill_manual(
    values = c("#1B9E77", "#D95F02", "#7570B3"),
    name   = "Country group"
  ) +
  coord_sf(xlim = c(2600000, 6200000), ylim = c(1450000, 5400000)) +
  labs(
    title    = "Geography of brokerage profiles (country clusters)",
    subtitle = "Group 1–3 based on broker-type mixture (Ward clustering, k = 3)",
    caption  = "Regions inherit their country's cluster."
  ) +
  theme_void(base_family = "serif") +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 16, hjust = 0.1, margin = margin(t = 10)),
    legend.title    = element_text(face = "bold")
  )

print(p_map_final)


# ===================================================================== --
# 34. Country-level Summaries for Export ----
# ===================================================================== --

# 1) Country-level descriptives
country_descriptives <- broker_features_df %>%
  group_by(main_ctry) %>%
  summarise(
    n_brokers         = n(),
    mean_deg          = mean(deg, na.rm = TRUE),
    median_deg        = median(deg, na.rm = TRUE),
    mean_kcore        = mean(kcore, na.rm = TRUE),
    median_kcore      = median(kcore, na.rm = TRUE),
    mean_constraint   = mean(constraint, na.rm = TRUE),
    median_constraint = median(constraint, na.rm = TRUE),
    mean_clust_co     = mean(clust_co, na.rm = TRUE),
    median_clust_co   = median(clust_co, na.rm = TRUE),
    mean_neigh_deg    = mean(neigh_deg, na.rm = TRUE),
    median_neigh_deg  = median(neigh_deg, na.rm = TRUE),
    .groups = "drop"
  )

# 2) Country-level broker-type shares (WIDE)
country_role_shares_wide <- broker_features_df %>%
  count(main_ctry, role_label, name = "n_role") %>%
  group_by(main_ctry) %>%
  mutate(share = n_role / sum(n_role)) %>%
  ungroup() %>%
  mutate(role_col = case_when(
    role_label == "Type 1: High clustering" ~ "share_type1_high_clust",
    role_label == "Type 2: High degree"     ~ "share_type2_high_deg",
    role_label == "Type 3: Low constraint"  ~ "share_type3_low_constr",
    role_label == "Type 4: Other"           ~ "share_type4_other",
    TRUE ~ paste0("share_", gsub("\\s+", "_", role_label))
  )) %>%
  select(main_ctry, role_col, share) %>%
  pivot_wider(
    names_from  = role_col,
    values_from = share,
    values_fill = 0
  )

# 3) Cluster assignment
country_clusters <- cluster_df %>%
  select(main_ctry, cluster_id, cluster_name) %>%
  distinct()

# 4) Final Summary
country_broker_summary_for_chatgpt <- country_descriptives %>%
  left_join(country_role_shares_wide, by = "main_ctry") %>%
  left_join(country_clusters, by = "main_ctry") %>%
  relocate(cluster_id, cluster_name, .before = main_ctry) %>%
  arrange(cluster_id, main_ctry)

print(country_broker_summary_for_chatgpt, n = nrow(country_broker_summary_for_chatgpt))

gc()

# ===================================================================== --
# 35. Maximum Spanning Tree (MST) Backbone Visualization (Straight Edges)
#      – keep only this "nice" backbone map ----
# ===================================================================== --

message("Generating Maximum Spanning Tree (MST) Backbone Map...")

# 1. Prepare the graph on NUTS2 edges (EPO-only)
if (!exists("edges_tl2")) stop("edges_tl2 object not found. Run Section 25 first.")

g_full <- graph_from_data_frame(
  d        = edges_tl2,
  directed = FALSE
)

g_full <- simplify(g_full, edge.attr.comb = "sum")

# 2. Maximum spanning tree via negative weights
E(g_full)$neg_weight <- -E(g_full)$weight
g_mst <- mst(g_full, weights = E(g_full)$neg_weight)

message(paste("MST calculated. Reduced edges from", ecount(g_full), "to", ecount(g_mst)))

# 3. Convert MST edges to straight spatial lines in LAEA
mst_edges_df <- igraph::as_data_frame(g_mst, what = "edges") %>%
  as_tibble() %>%
  select(from, to, weight)

if (!exists("nuts2_sf")) stop("nuts2_sf object not found. Run Section 26 first.")
if (!exists("nodes_tl2")) stop("nodes_tl2 object not found. Run Section 26 first.")

crs_laea <- 3035

# Region classification: high-activity vs others
nodes_tl2 <- nodes_tl2 %>%
  mutate(
    region_class = if_else(
      strength >= quantile(strength, 0.75, na.rm = TRUE),
      "High-activity",
      "Other region"
    )
  )

# NUTS2 polygons in LAEA + region_class
nuts2_laea <- nuts2_sf %>%
  st_transform(crs_laea) %>%
  left_join(nodes_tl2 %>% select(nuts_id, region_class), by = "nuts_id")

nuts2_col <- nuts2_laea %>% filter(!is.na(region_class))
nuts2_na  <- nuts2_laea %>% filter(is.na(region_class))

# Centroids in LAEA
centroids_laea <- st_point_on_surface(nuts2_laea) %>%
  select(nuts_id, geometry)

# Rebuild MST edges in LAEA (straight lines)
mst_edges <- mst_edges_df %>%
  left_join(centroids_laea, by = c("from" = "nuts_id")) %>%
  rename(geom_from = geometry) %>%
  left_join(centroids_laea, by = c("to" = "nuts_id")) %>%
  rename(geom_to = geometry) %>%
  filter(!st_is_empty(geom_from), !st_is_empty(geom_to))

mst_lines <- st_sfc(
  map(1:nrow(mst_edges), function(i) {
    coords <- rbind(
      st_coordinates(mst_edges$geom_from[i]),
      st_coordinates(mst_edges$geom_to[i])
    )
    st_linestring(coords)
  }),
  crs = crs_laea
)

mst_laea <- st_sf(
  from     = mst_edges$from,
  to       = mst_edges$to,
  weight   = mst_edges$weight,
  geometry = mst_lines
)

# Weight classes for slide-style graphic
mst_laea <- mst_laea %>%
  mutate(
    w_class = cut(
      weight,
      breaks = c(-Inf, 10, 100, 1000, Inf),
      labels = c("1", "10", "100", "1000"),
      include.lowest = TRUE
    )
  )

bg_eu <- gisco_get_countries(region = "Europe") %>% st_transform(crs_laea)

p_backbone <- ggplot() +
  geom_sf(data = bg_eu, fill = "grey95", color = "white", linewidth = 0.2) +
  geom_sf(data = nuts2_na, fill = "grey90", color = "white", linewidth = 0.1) +
  geom_sf(data = nuts2_col, aes(fill = region_class), color = "white", linewidth = 0.1) +
  geom_sf(data = mst_laea, aes(size = w_class), color = "black", alpha = 0.55) +
  scale_fill_manual(
    values = c("High-activity" = "#e85c47", "Other region" = "grey80"),
    name   = "Region type"
  ) +
  scale_size_manual(
    values = c("1" = 0.1, "10" = 0.25, "100" = 0.55, "1000" = 1.2),
    name   = "Collaboration ties"
  ) +
  coord_sf(
    xlim   = c(2500000, 7000000),
    ylim   = c(1500000, 5500000),
    expand = FALSE
  ) +
  labs(
    title    = "Backbone of the European co-inventor network",
    subtitle = "Maximum spanning tree (NUTS2-level collaboration ties)",
    caption  = "Edges classified by collaboration intensity."
  ) +
  theme_void(base_size = 16, base_family = "serif") +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 13),
    legend.text     = element_text(size = 11),
    plot.title      = element_text(face = "bold", size = 16, hjust = 0.5,
                                   margin = margin(t = 10, b = 5)),
    plot.subtitle   = element_text(color = "grey40", hjust = 0.5,
                                   margin = margin(b = 15)),
    plot.caption    = element_text(color = "grey60", size = 8, hjust = 0.9)
  )

print(p_backbone)
