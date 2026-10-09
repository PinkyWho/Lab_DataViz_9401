library(ggplot2)
library(dplyr)
library(patchwork)

# Load the pre-calculated data back into memory instantly
load("coherence_results.RData")

# [Paste the corrected make_subplot function here]
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
# [Paste the Step 4 for-loop here to generate the plots]
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