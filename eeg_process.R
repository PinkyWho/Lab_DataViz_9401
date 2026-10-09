# ==============================================================================
# R EEG Processing & Coherence Pipeline (Using MNE via Reticulate)
# ==============================================================================

library(reticulate)

# Specify python version to use. Originally ran with Python 3.9
# 1. Define local venv path in your working directory
venv_dir <- file.path(getwd(), "r_eeg_env")

# 2. Create the environment if it doesn't exist yet
if (!virtualenv_exists(venv_dir)) {
  cat("Creating new R-dedicated virtual environment...\n")
  virtualenv_create(envname = venv_dir)
  
  # 3. Install required Python dependencies into this environment
  cat("Installing required Python packages (mne, mne-connectivity, numpy, scipy)...\n")
  virtualenv_install(
    envname = venv_dir,
    packages = c("mne", "mne-connectivity", "numpy", "scipy", "h5py", "pandas"),
    ignore_installed = FALSE
  )
}

# 4. Lock reticulate to this specific virtual environment
use_virtualenv(venv_dir, required = TRUE)

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# 1. Setup Python Environment & Parameters (No native R Coherence method)
# ------------------------------------------------------------------------------

# Import Python modules into R
mne <- import("mne")
mne_conn <- import("mne_connectivity")
np <- import("numpy")

# Editable to remove more outliers if need be
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
freqs <- seq(4, 30, by = 0.1) # 4 to 30 Hz
n_cycles <- freqs / 2         # Time-frequency tradeoff

# Convert R vectors to explicit NumPy arrays for MNE compatibility
freqs_np <- np$array(freqs)
n_cycles_np <- np$array(n_cycles)

# Read metadata
metad <- read_excel("NeuroTechs Dataset for Stem Skills/extra_metadata.xlsx", sheet = "Individual metadata")

all_coh <- list()
all_sex <- c()
all_age <- c()

# ------------------------------------------------------------------------------
# 2. Iterate Subjects, Epoch, and Compute Coherence
# ------------------------------------------------------------------------------
cat("Beginning processing...\n")

for (fullname in sublist) {
  cat("Processing:", fullname, "\n")
  
  # Metadata
  sub_row <- metad %>% filter(row_number() == match(fullname, metad[[1]]))
  if (nrow(sub_row) == 0) next
  age <- sub_row$Age[1]   
  sex <- sub_row$`AAB Sex`[1]
  
  # Paths
  events_file <- file.path("NeuroTechs Dataset for Stem Skills", fullname, "programming_responses.csv")
  eeg_file    <- file.path("NeuroTechs Dataset for Stem Skills", fullname, "ses-1", "eeg", paste0(fullname, "_ses-1_task-STEMSKILLS_eeg.set"))
  chan_file   <- file.path("NeuroTechs Dataset for Stem Skills", fullname, "ses-1", "eeg", paste0(fullname, "_ses-1_electrodes.tsv"))
  
  if (file.exists(eeg_file) && file.exists(events_file)) {
    # Read events and filter
    events_df <- read.csv(events_file)
    question_appearance <- events_df[, 3] - 3 # - 3 because event timeline is off by 3 seconds from the preprocessing pipeline
    
    # Questions that were answered in less than 2 seconds are not considered
    diffs <- diff(question_appearance)
    valid_mask <- c(TRUE, diffs >= 2)
    filtered_times <- question_appearance[valid_mask]
    event_samples <- as.integer(filtered_times * fs)
    
    # Create MNE-compatible event matrix in R
    mne_events <- cbind(event_samples, integer(length(event_samples)), rep(1L, length(event_samples)))
    storage.mode(mne_events) <- "integer" # MNE requires strictly integers
    
    # Load and preprocess raw EEG via MNE
    raw <- mne$io$read_raw_eeglab(eeg_file, preload = TRUE, verbose = FALSE)
    montage <- mne$channels$read_custom_montage(chan_file)
    raw$set_eeg_reference('average', verbose = FALSE)
    raw$filter(4, 30, verbose = FALSE)
    
    # Epoch data (-1 to 2 seconds)
    epochs <- mne$Epochs(raw, mne_events, event_id = 1L, tmin = -1, tmax = 2, verbose = FALSE)
    
    # Calculate Coherence (Python mne)
    con <- mne_conn$spectral_connectivity_epochs(
      epochs,
      method = 'imcoh',
      mode = 'cwt_morlet',
      sfreq = fs,
      cwt_freqs = freqs_np,
      cwt_n_cycles = n_cycles_np,
      fmin = 4,
      fmax = 30,
      tmin = -1.0,
      verbose = FALSE
    )
    
    # Extract dense array back into R environment
    con_dense <- con$get_data(output = 'dense')
    
    # Clip and Fisher Z Transform (1's and 0's are disregarded)
    con_dense <- pmax(pmin(con_dense, 0.999999), -0.999999)
    con_dense_z <- atanh(con_dense)
    
    all_coh[[length(all_coh) + 1]] <- con_dense_z
    all_sex <- c(all_sex, sex)
    all_age <- c(all_age, age)
  }
}

# ------------------------------------------------------------------------------
# 3. Group Means & Statistics
# ------------------------------------------------------------------------------
cat("Testing significance...\n")

male_idx   <- which(all_sex == "Male")
female_idx <- which(all_sex == "Female")

# Stack into 5D array: [Channel, Channel, Freq, Time, Subject]
# Note: Python shapes come into R slightly differently. MNE's dense output is [Ch, Ch, Freq, Time]
coh_5d <- array(unlist(all_coh), dim = c(dim(all_coh[[1]]), length(all_coh)))

# Averages across subjects
Ccoh  <- apply(coh_5d, c(1, 2, 3, 4), mean, na.rm = TRUE)
Ccohm <- apply(coh_5d[, , , , male_idx], c(1, 2, 3, 4), mean, na.rm = TRUE)
Ccohf <- apply(coh_5d[, , , , female_idx], c(1, 2, 3, 4), mean, na.rm = TRUE)

# Stats: t-tests (p < 0.05 threshold)
n_times <- dim(Ccoh)[4]
n_freqs <- dim(Ccoh)[3]
pvals <- array(1, dim = dim(Ccoh)[1:4])

for (i in 1:length(channellist)) {
  for (j in 1:length(channellist)) {
    if (i == j) next # Skip self-coherence
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
sig_mask <- pvals < 0.05

# Save the heavy calculated objects to a file
save(Ccoh, Ccohm, Ccohf, sig_mask, time_axis, freq_axis, n_times, n_freqs, channellist, file = "coherence_results.RData")

# ------------------------------------------------------------------------------
# 4. Generate 3-Panel ggplot2 Figures
# ------------------------------------------------------------------------------
cat("Rendering plots...\n")

# Recreate time and frequency axes based on final matrix dimensions
time_axis <- seq(-1, 2, length.out = n_times)
freq_axis <- seq(4, 30, length.out = n_freqs)

make_subplot <- function(mat, sig_mat, title_str, show_points = FALSE) {
  # FIXED: FreqIdx must come first so the grid matches matrix flattening
  df <- expand.grid(FreqIdx = 1:n_freqs, TimeIdx = 1:n_times) %>%
    mutate(
      Time = time_axis[TimeIdx],
      Freq = freq_axis[FreqIdx],
      Coherence = as.vector(mat),
      Sig = as.vector(sig_mat)
    )
  
  p <- ggplot(df, aes(x = Time, y = Freq, fill = Coherence)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_distiller(palette = "RdBu", limits = c(-0.2, 0.2), oob = scales::squish) +
    labs(title = title_str, x = "Time (s)", y = "Frequency (Hz)", fill = "Coherence") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 9)
    ) +
    coord_cartesian(expand = FALSE)
  
  if (show_points) {
    p <- p + geom_point(data = filter(df, Sig == TRUE), 
                        aes(x = Time, y = Freq), color = "cyan", size = 0.5, alpha = 0.7)
  }
  return(p)
}

# Plot all electrode pairs
for (i in 2:(length(channellist))) {
  for (j in 1:(i - 1)) { 
    
    p_all  <- make_subplot(Ccoh[i, j, , ], sig_mask[i, j, , ], paste0(channellist[i], " to ", channellist[j], " (All)"), show_points = TRUE)
    p_male <- make_subplot(Ccohm[i, j, , ], sig_mask[i, j, , ], paste0(channellist[i], " to ", channellist[j], " (Male)"))
    p_fem  <- make_subplot(Ccohf[i, j, , ], sig_mask[i, j, , ], paste0(channellist[i], " to ", channellist[j], " (Female)"))
    
    # Assemble side-by-side using Patchwork, sharing one colorbar!
    combined_plot <- p_all + p_male + p_fem + plot_layout(guides = "collect")
    
    filename <- paste0(channellist[i], "_", channellist[j], "_combined_0.05.png")
    ggsave(filename, plot = combined_plot, path = "figs",width = 14, height = 2, dpi = 300)
  }
}
cat("Pipeline complete.\n")