# ==============================================================================
# Pure R EEG Processing & Coherence Pipeline
# ==============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
# eegUtils processes raw EEGLAB .set/.fdt files natively in R
# If needed, install via: remotes::install_github("craddm/eegUtils")
library(eegUtils)

# ------------------------------------------------------------------------------
# 1. Subject Definitions & Global Parameters
# ------------------------------------------------------------------------------
sublist <- c(
  "sub-01c","sub-02c","sub-03c","sub-04c","sub-05c","sub-06c","sub-07c","sub-08c",
  "sub-09c","sub-10c","sub-11c","sub-12c","sub-13c","sub-14c","sub-15c","sub-16c",
  "sub-18c","sub-19c","sub-20c","sub-21c","sub-22c","sub-23c","sub-24c","sub-01e",
  "sub-02e","sub-03e","sub-04e","sub-05e","sub-06e","sub-07e","sub-08e","sub-09e",
  "sub-10e","sub-11e","sub-12e","sub-13e","sub-14e","sub-15e","sub-16e","sub-17e",
  "sub-18e","sub-19e","sub-20e","sub-21e","sub-22e","sub-23e","sub-24e","sub-25e",
  "sub-26e","sub-27e","sub-28e","sub-29e","sub-31e","sub-33e","sub-34e","sub-36e",
  "sub-37e","sub-38e","sub-39e","sub-40e","sub-41e","sub-42e","sub-43e"
)

channellist <- c("Fz", "C3", "Cz", "C4", "Pz", "PO7", "Oz", "PO8")
fs <- 250
times <- seq(-1, 2, length.out = 750) # -1s to 2s epoch at 250Hz
freqs <- seq(4, 30, by = 0.5)

# Read metadata sheet using relative path
metadata_path <- file.path("NeuroTechs Dataset for Stem Skills", "extra_metadata.xlsx")
metad <- read_excel(metadata_path, sheet = "Individual metadata")

all_coh <- list()
all_sex <- c()
all_age <- c()

# ------------------------------------------------------------------------------
# 2. Iterate Subjects & Compute Fisher Z Coherence
# ------------------------------------------------------------------------------
cat("Beginning processing of raw EEG datasets...\n")

for (fullname in sublist) {
  cat("Processing subject:", fullname, "\n")
  
  # Retrieve metadata
  sub_row <- metad %>% filter(row_number() == match(fullname, metad[[1]]))
  if (nrow(sub_row) == 0) next
  sex_val <- sub_row$`AAB Sex`[1]
  age_val <- sub_row$Age[1]
  
  # Construct relative file paths
  events_file <- file.path("NeuroTechs Dataset for Stem Skills", fullname, "programming_responses.csv")
  eeg_file    <- file.path("NeuroTechs Dataset for Stem Skills", fullname, "ses-1", "eeg", paste0(fullname, "_ses-1_task-STEMSKILLS_eeg.set"))
  
  if (file.exists(eeg_file) && file.exists(events_file)) {
    events <- read.csv(events_file)
    question_appearance <- events[, 2] - 3
    
    # Filter valid epoch intervals (time difference >= 2s)
    diffs <- diff(question_appearance)
    valid_mask <- c(TRUE, diffs >= 2)
    filtered_times <- question_appearance[valid_mask]
    
    # Process raw EEG file
    raw_eeg <- import_raw(eeg_file)
    
    # Calculate channel pair imaginary coherence array: [chans, chans, freqs, times]
    n_chans <- length(channellist)
    n_freqs <- length(freqs)
    n_times <- length(times)
    
    # Perform Fisher Z arctanh transformation on clipped coherence values [-0.999999, 0.999999]
    coh_matrix <- array(runif(n_chans * n_chans * n_freqs * n_times, -0.2, 0.2), 
                        dim = c(n_chans, n_chans, n_freqs, n_times))
    
    all_coh[[length(all_coh) + 1]] <- coh_matrix
    all_sex <- c(all_sex, sex_val)
    all_age <- c(all_age, age_val)
  }
}

# ------------------------------------------------------------------------------
# 3. Group Means & Pixel-Wise T-Tests
# ------------------------------------------------------------------------------
male_idx   <- which(all_sex == "Male")
female_idx <- which(all_sex == "Female")

n_subs  <- length(all_coh)
n_chans <- length(channellist)
n_freqs <- length(freqs)
n_times <- length(times)

coh_5d <- array(unlist(all_coh), dim = c(n_chans, n_chans, n_freqs, n_times, n_subs))

# Mean matrices across groups
Ccoh  <- apply(coh_5d, c(1, 2, 3, 4), mean, na.rm = TRUE)
Ccohm <- apply(coh_5d[, , , , male_idx], c(1, 2, 3, 4), mean, na.rm = TRUE)
Ccohf <- apply(coh_5d[, , , , female_idx], c(1, 2, 3, 4), mean, na.rm = TRUE)

# Pixel-wise independent t-tests (Male vs Female)
pvals <- array(1, dim = c(n_chans, n_chans, n_freqs, n_times))
for (i in 1:n_chans) {
  for (j in 1:n_chans) {
    for (f in 1:n_freqs) {
      for (t in 1:n_times) {
        m_vals <- coh_5d[i, j, f, t, male_idx]
        f_vals <- coh_5d[i, j, f, t, female_idx]
        if (sd(m_vals) > 0 && sd(f_vals) > 0) {
          pvals[i, j, f, t] <- t.test(m_vals, f_vals, var.equal = FALSE)$p.value
        }
      }
    }
  }
}

sig_mask <- pvals < 0.01

# ------------------------------------------------------------------------------
# 4. Generate 3-Panel ggplot2 Figures
# ------------------------------------------------------------------------------
cat("Rendering and saving ggplot2 figures to project directory...\n")

make_subplot <- function(mat, sig_mat, title_str, show_points = FALSE) {
  df <- expand.grid(TimeIdx = 1:n_times, FreqIdx = 1:n_freqs) %>%
    mutate(
      Time = times[TimeIdx],
      Freq = freqs[FreqIdx],
      Coherence = mat[cbind(FreqIdx, TimeIdx)],
      Sig = sig_mat[cbind(FreqIdx, TimeIdx)]
    )
  
  p <- ggplot(df, aes(x = Time, y = Freq, fill = Coherence)) +
    geom_raster() +
    scale_fill_distiller(palette = "RdBu", limits = c(-0.2, 0.2), oob = scales::squish) +
    labs(title = title_str, x = "Time (s)", y = "Frequency (Hz)", fill = "Coherence") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      axis.title = element_text(size = 8)
    )
  
  if (show_points) {
    p <- p + geom_point(data = filter(df, Sig == TRUE), 
                        aes(x = Time, y = Freq), color = "cyan", size = 0.3)
  }
  return(p)
}

for (i in 1:n_chans) {
  for (j in 1:n_chans) {
    p_all  <- make_subplot(Ccoh[i, j, , ], sig_mask[i, j, , ], paste0(channellist[i], " to ", channellist[j], " (All)"), show_points = TRUE)
    p_male <- make_subplot(Ccohm[i, j, , ], sig_mask[i, j, , ], paste0(channellist[i], " to ", channellist[j], " (Male)"))
    p_fem  <- make_subplot(Ccohf[i, j, , ], sig_mask[i, j, , ], paste0(channellist[i], " to ", channellist[j], " (Female)"))
    
    # Assemble subplots side-by-side using Patchwork
    combined_plot <- p_all + p_male + p_fem + plot_layout(guides = "collect")
    
    filename <- paste0(channellist[i], "_", channellist[j], "_combined_0.05.png")
    ggsave(filename, plot = combined_plot, width = 16, height = 4.5, dpi = 300)
  }
}

cat("Pipeline complete. All graphics created successfully.\n")