# SeqRate, a novel pipeline for RNA polymerase II elongation rate estimation

library(devtools)
library(gsignal)
library(GenomicFeatures)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rtracklayer)
library(glue)
library(DescTools)
library(data.table)
library(zoo)
library(MASS)
library(stringr)
library(scales)
#-----------------------------------------------------------------------------

#Set directories


# Download GTEx median TPMs:
gtex_url <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
gtex_data <- fread(gtex_url)

#Get genes expressed in many tissues

#Set a threshold for tpm count
tpm_threshold <- 0.5

gtex_data_filtered <- gtex_data %>%
  rowwise() %>%
  dplyr::mutate(
    tissues_expressed = sum(c_across(3:70) > tpm_threshold),
    median_tpm = median(c_across(3:70)),
    gene_id = sub("\\.\\d+$", "", Name)
  ) %>%
  ungroup() %>%
  dplyr::mutate(
    percentile_score = percent_rank(tissues_expressed) * 0.5 + 
      percent_rank(median_tpm) * 0.5,
    rank = ntile(percentile_score, 100)
  ) %>%
  dplyr::filter( rank > 60)
    
#Plot percentile scores
ggplot(gtex_data_filtered, aes(x = percentile_score)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Percentile Scores", x = "Percentile Score", y = "Frequency") +
  theme_minimal()

#Load annotation from ENCODE
full_annotation <- import("~/seqrate/code/annotation.gtf.gz")
full_annotation <- as.data.frame(full_annotation)


#Keep only annotation on protein coding genes, with no duplicate genes
full_annotation <- full_annotation %>%
  dplyr::mutate(seqnames = paste0("chr", seqnames)) %>%
  dplyr::filter(!is.na(gene_name),
         gene_biotype == "protein_coding"
         ) %>%
semi_join(
  gtex_data_filtered, by = "gene_id" 
)

#Create a dataframe which contains annotation on whole genes
genes_annotation <- full_annotation %>%
  dplyr::filter(type == "gene")  
  
#Create a df which only contains annotation on exons
introns_annotation <- full_annotation %>%
  dplyr::filter(type == "exon") %>%

  #Keep only exons on the same strand with the strand on the gene annotation
  group_by(gene_id) %>%
  dplyr::filter(strand == genes_annotation$strand[match(gene_id, genes_annotation$gene_id)]) %>%
  arrange(if_else(strand == "-",dplyr::desc(row_number()),row_number())) %>%
  ungroup() %>%
  
  #Generate intron coordinates from exons of each transcript
  group_by(gene_id, transcript_id) %>%
  dplyr::mutate(
    intron_annotation_start = end + 1,
    annotation_intron_end = lead(start)-1,
    annotation_intron_length = annotation_intron_end - intron_annotation_start,
  ) %>%
 
 #Remove introns with NA length values 
  dplyr::filter(
    !is.na(annotation_intron_length)
  )
  
  #Calculate mean intron length and stn deviation for each transcript   
  transcripts_annotation <- introns_annotation %>%
  summarise(
    intron_annotation_count = n(),
    annotation_mean_length = mean(annotation_intron_length),
    annotation_sd_length = sd(annotation_intron_length)
  ) %>%
  ungroup() %>% 
  
  #Keep only transcripts with a good number of introns, where many are long  
  dplyr::filter(!is.na(annotation_sd_length)) %>%
  dplyr::filter(
    intron_annotation_count >= 3,
    intron_annotation_count <= 10,
    annotation_mean_length >= quantile(annotation_mean_length, 0.5),
    annotation_sd_length <= quantile(annotation_sd_length, 0.8)
  )%>%
  group_by(gene_id) %>%
  
  summarise(long_intron_transcripts = n())

#Density plot
library(ggplot2)

transcripts_annotation %>%
  ggplot(aes(x = long_intron_transcripts)) +
  geom_density(fill = "#56B4E9", alpha = 0.6, color = "black", linewidth = 0.6) +
  scale_x_continuous(limits = c(0, 5), breaks = 0:5) +
  labs(
    title = "Distribution of Long Intron Transcripts per Gene",
    x = "Number of Long Intron Transcripts",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray85")
  )

#Keep genes with high possibility of exhibiting long intron transcripts with alt splicing
#transcripts_annotation <- transcripts_annotation %>%
#  dplyr::filter(
#    long_intron_transcripts >= 2
#  )

genes_annotation_filtered <- genes_annotation %>%
  semi_join(transcripts_annotation, by = "gene_id") %>%
  dplyr::select(-source,
                -score,
                -phase,
                -width,
                -type,
                -gene_biotype,
                -transcript_id,
                -transcript_biotype,
                -transcript_name,
                -exon_number,
                -exon_id,
                -exon_version,
                -transcript_support_level,
                -protein_id,
                -ccds_id,
                -gene_version,
                -transcript_version,
                -tag,
                -transcript_source,
                -exon_version,
                -protein_version,
                -gene_source
                ) %>%
  dplyr::rename(chromosome = seqnames) %>%

#Remove any overlaping genes
  dplyr::group_by(chromosome) %>%
  dplyr::filter(start > cummax(lag(end, default = first(start)))) %>%
  
  ungroup()

#Keep only annotation of selected genes
introns_annotation <- introns_annotation %>%
  semi_join(genes_annotation_filtered, by = "gene_id") 


#-----------------------------------------------------------------------------------------
#Set bin length and sorting parameters
bin_length <- 100
chromosome_order <- paste0("chr", c(1:22,"X", "Y"))

#Sort by chromosome and by gene starting coordinates
df_for_counts <- genes_annotation_filtered %>%
  mutate(chromosome = factor(as.character(chromosome), levels = chromosome_order)) %>%
  arrange(chromosome, start) %>%
  


#If a gene is smaller that bin size keep it as is, else create bins of desired length 
 rowwise() %>% 
    dplyr::mutate(
    bin_starts = ifelse(end-start < bin_length, list(start),list(seq(from = as.numeric(start), 
                          to = as.numeric(end), 
                          by = bin_length + 1))),
    bin_ends = ifelse(end - start < bin_length, 
                      list(end), 
                      list(c(bin_starts[-1] - 1, end)))
  ) %>%
  dplyr::select(-start, -end) %>%
  unnest(cols = c(bin_starts, bin_ends)) %>%
  dplyr::mutate(
    X2 = bin_starts,
    X3 = bin_ends
  ) %>%
  dplyr::select(-bin_starts, -bin_ends) %>%
  relocate(X2, .after = chromosome) %>%
  relocate(X3, .after = X2) %>%
  mutate(X2 = as.integer(X2), X3 = as.integer(X3)) %>%
  dplyr::rename(
    bin_start = X2,
    bin_end = X3
  )


#Convert bins and introns to GRanges
bins_gr <- GRanges(
  seqnames = df_for_counts$chromosome,
  ranges = IRanges(start = df_for_counts$bin_start, 
                   end = df_for_counts$bin_end)
)

introns_gr <- GRanges(
  seqnames = introns_annotation$seqnames,
  ranges = IRanges(start = introns_annotation$intron_annotation_start,
                   end = introns_annotation$annotation_intron_end)
)

#Find overlaps between bins and introns
overlaps <- findOverlaps(bins_gr, introns_gr, type = "any")

#Keep only unique bin rows wich are intronic
overlapping_bins <- unique(queryHits(overlaps))
intronic_bins <- df_for_counts[overlapping_bins, ]

#For bins overlapping more than one intron assign them to the first intron.
first_overlaps <- overlaps[!duplicated(queryHits(overlaps))]
intronic_bins$intron_index <- subjectHits(first_overlaps)

#Save to BED file
write.table(
  intronic_bins,
  "~/seqrate/code/bins.bed",
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

#---------------------------------------------------------------------------------

#Open counts files of the first sample from bedtools
all_genes_table<-read_delim(file="~/seqrate/counts/huvec_early3_counts.txt",
                            delim = "\t",
                            col_types = "cnncccnn",
                            comment="",
                            col_names = FALSE)

#add the name of the first sample (start with the word counts)
f_sample_name = "counts_huvec_early"

#Rename columns 
df_bins <- as.data.frame(all_genes_table) %>%
  dplyr::rename(
    chromosome = X1,
    start = X2,
    end = X3,
    strand = X4,
    gene_id = X5,
    gene_name = X6,
    intron_id = X7,
    !!f_sample_name := X8
  )


#Add counts from other samples
add_counts_columns <- function(file_paths,column_names,counts_data) {
  for (i in 1:length(file_paths)) {
    
    counts_data_temp <- read_delim(
      file = file_paths[i],
      delim = "\t",
      col_types = "nnncccnn",
      comment = "",
      col_names = FALSE
    )
    
    counts_data[[column_names[i]]] <- counts_data_temp$X8
    
  }
  
  return(counts_data)
}  

df_bins <- add_counts_columns(
  file_paths = "~/seqrate/counts/huvec_late3_counts.txt",
  column_names = "counts_huvec_late",
  counts_data = df_bins
)

#Redefine bin size
bin_size <- 100

#Removing "flat" introns    
  df_bins <- df_bins %>%
  group_by(gene_name,intron_id) %>%
    dplyr::mutate(across(starts_with("counts"), sd, .names = "stn_dev_{.col}"))%>%
    dplyr::filter(if_all(starts_with("stn_dev"), ~!is.na(.x))) %>%
    dplyr::filter(if_all(starts_with("stn_dev"), ~ .x >= 1.5)) %>%
    ungroup() %>%
    

#Set location in gene for each intron    
  group_by(gene_name) %>%
  
#Ignore genes with no counts  
  dplyr::filter(if_all(starts_with("counts"), ~ all(.x != 0))) %>%

#Ignore genes with less than 200 bins left
  dplyr::filter(length(start) > 200) %>%
    
#Ignore genes with less low read coverage counts
    dplyr::filter(if_all(starts_with("counts"), ~ mean(.x >40) >= 0.8 )) %>%

#Set a specific id for each gene based on gene position
  dplyr::mutate(
    intron_rank = dense_rank(intron_id)) %>%
  ungroup()
  
 #Create a df with info on each intron coords   
df_introns_coords <- df_bins %>%
  group_by(chromosome, gene_name , intron_rank,intron_id) %>%
  summarize(
    intron_begin = min(start),     
    intron_finish = max(end),
    intron_length = intron_finish - intron_begin,
    .groups ='drop'
)

#Trimming
df_bins <- df_bins %>%
  left_join(
    df_introns_coords %>% dplyr::select(gene_name, intron_id, intron_length),
    by = c("gene_name", "intron_id")
  ) %>%
  dplyr::mutate(
    bins_in_intron = ceiling(intron_length / bin_size),
    intron_index = 1:n(),
    across(starts_with("counts"), ~ predict(MASS::rlm(.x ~ intron_index,
                                     psi = psi.bisquare,
                                     maxit=5000)),
           .names = "linear_trend_{.col}")
  )%>%
  dplyr::mutate(
    across(starts_with("counts"),
           ~ .x - get(paste0("linear_trend_", cur_column())),
           .names = "detrended_{.col}")
  ) %>%
  ungroup()

  trim_peaks <- function(data,window_multiplier){
    data %>%
      group_by(gene_name, intron_id) %>%
      dplyr::mutate(
        window_size = max(7, ceiling(window_multiplier * bins_in_intron)),
        across(starts_with("detrended_counts"),
            ~{
        window_median = zoo::rollapply(
            .x,
            width = window_size,
            FUN = median,
            fill = NA,
            align = "center",
            partial = TRUE
          )
        
         window_mad = zoo::rollapply(
            .x,
            width = window_size,
            FUN = \(x) mad(x, constant = 1),
            fill = NA,
            align = "center",
            partial = TRUE
        )
         
        if_else(
          .x > (window_median + 3 * window_mad),
        window_median,
        .x
          )
        },
    .names= "{.col}"
      )    
    ) %>%
    dplyr::select(-window_size)
  }
  
  window_multiplier <- c(0.10, 0.05, 0.02, 0.01)
for (multiplier in window_multiplier) {
  df_bins <- trim_peaks(df_bins,multiplier)
}

df_bins <- df_bins %>%
  dplyr::mutate(
    across(starts_with("detrended_counts"),
      ~ .x + get(paste0("linear_trend_", str_remove(cur_column(), "detrended_"))),
    .names = "trimmed_{str_remove(.col, 'detrended_')}" 
    )
  )%>%
  dplyr::select(-intron_index, -starts_with("linear_trend"), -starts_with("detrended_counts"))

#Create a complete df for introns
df_introns <- df_bins %>%
  group_by(chromosome,gene_name, intron_rank,intron_id) %>%
  summarize(
    intron_begin = min(start),     
    intron_finish = max(end),
    across(starts_with("trimmed_counts"),
           ~ max(.x),
           .names = "max_{.col}"),
    across(starts_with("trimmed_counts"),
           ~min(.x),
           .names = "min_{.col}"
           ),
    across(starts_with("stn_dev"),
           ~min(.x),
           .names = "{.col}_intron"
          ),
    strand = first(strand),
    .groups ='drop'
  ) %>%
  dplyr::mutate(intron_lengths = intron_finish - intron_begin) %>%
  dplyr::filter(if_all(
    starts_with("stn_dev") & ends_with("intron"),
    ~.x >= 1 | !is.na(.x))) %>%
  group_by(chromosome) %>%
  arrange(intron_begin, .by_group = TRUE) %>%
  ungroup()

#--------------------------------------------------------------------------------------
#Generate the sawtooth model for the expression of each gene and create plots
gene_names_from_counts <- c(unique(df_introns$gene_name))  

df_genes <- data.frame(gene_names_from_counts)
score_sample <- vector("numeric",length=length(gene_names_from_counts))
gene_number <- 1

#Vector to store the MRE
count_cols <- grep("^trimmed_counts", names(df_bins), value = TRUE)
short_names <- sub("^trimmed_counts", "", count_cols)

mre_scores <- matrix(  
  nrow = length(gene_names_from_counts),
  ncol = length(count_cols),
  dimnames = list(
    NULL,
    short_names
  )
)

#Setup a dir for the plots
output_dir <- "~/seqrate/plots/huvec3"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
pdf(file.path(output_dir, "all_plots.pdf"), width = 12, height = 9)

#Set the grid for the plots
par(mar=c(4, 4, 2, 2))
par(mfrow = c(3,4))
plot_counter <- 0

ggplot(df_introns, aes(x = intron_lengths)) +
  geom_density(fill = "skyblue", alpha = 0.7) +
  scale_x_continuous(
    limits = c(0, 1e5),  # 0-100K bp range
    labels = label_number(scale = 1e-3, suffix = "K"),  # Show as "0K" to "100K"
    breaks = seq(0, 1e5, by = 2e4)  # Ticks every 20K
  ) +
  labs(x = "Intron Length (0-100K bp)", y = "Density") +
  theme_minimal()





#Loop through all of the genes

gene_counter <- 1
for(g in gene_names_from_counts) {
  g_unique <- which(as.character(df_bins$gene_name) == g) 
  i_unique <- which(as.character(df_introns$gene_name) == g) 
  
  intron_length <- as.numeric(df_introns$intron_lengths[i_unique])
  introns_per_gene <- as.numeric(max(df_introns$intron_rank[i_unique]))  
  one_g_unique <- min(as.numeric(g_unique))
  strands <- df_bins[one_g_unique,4]
  
  intron_starts_gene <- df_introns$intron_begin[i_unique]
  intron_starts_gene <- as.numeric(intron_starts_gene)
  intron_ends_gene <-df_introns$intron_finish[i_unique]
  intron_ends_gene <- as.numeric(intron_ends_gene)
  

  #Process each counts column
  for (count_col in count_cols ) {
    max_col <- paste0("max_",count_col)
    min_col <- paste0("min_",count_col)
    max_counts_gene <- as.numeric(df_introns[[max_col]][i_unique])
    min_counts_gene <- as.numeric(df_introns[[min_col]][i_unique])
  
  #Create a sawtooth model
  waveform <- unlist(lapply(1:introns_per_gene, function(x){ 
    intron_length_set <- intron_length[x]
    max_count_intron <- max_counts_gene[x]
    min_count_intron <- min_counts_gene[x]
    intron_start <- intron_starts_gene[x]
    intron_end <- intron_ends_gene[x]
    
    t <- seq(0, 1,length.out = intron_length_set)
    strand <- strands
    if(strand == "+"){
      sawtooth_wave <- max_count_intron - ((max_count_intron - min_count_intron) * (sawtooth(2 * pi * t) + 1) / 2)
    }else{
      sawtooth_wave <-min_count_intron + ((max_count_intron-min_count_intron) * (sawtooth(2 * pi  * t)+1)/2)
    }
    return(sawtooth_wave)
  }))
  

  
  #PLot RNA-seq data 
  z<-seq(1:length(g_unique))
  y<-as.numeric(unlist(df_bins[g_unique, count_col]))
  
  plot(z,y, 
       type="h"
       ,col = "#1f77b4",
       main=g,
       xlab = "Position in Gene (Bins)",
       ylab = "RNA-seq read counts")
  
  start_length <- length(waveform)
  end_length <- length(y)
  

      #Fit the sawtooth wave model on the data
  downscaled_wave <- approx(
    x = seq(0, 1, length.out = start_length), 
    y = waveform,                                    
    xout = seq(0, 1, length.out = end_length))$y    
  lines(downscaled_wave,col="black")

  plot_counter <- plot_counter + 1  
    # Add a new page if the grid is filled
  if (plot_counter %% 12 == 0) {  
    par(mfrow = c(3, 4))
    plot_counter <- 0 
  }
  #Score each gene based on how well it fits on the sawtooth model
  
  #Normalization to [0,1]
  
  y_normalized <- (y - min(y) + 1) / (max(y) - min(y) + 1)  # Add pseudocount
  sawtooth_normalized <- (downscaled_wave - min(downscaled_wave) + 1) / 
    (max(downscaled_wave) - min(downscaled_wave) + 1)
  mre = mean(abs(y_normalized-sawtooth_normalized)/ (abs(sawtooth_normalized) + 1e-10))
  mre_scores[gene_counter, sub("^trimmed_counts", "", count_col)] <- mre
  
  gene_idx <- which(gene_names_from_counts == g)
  col_idx <- which(short_names == sub("^trimmed_counts_", "", count_col))
  mre_scores[gene_idx, col_idx] <- mre
  
  } 
  gene_counter <- gene_counter + 1  
}



df_genes <- cbind(df_genes, as.data.frame(mre_scores))
mre_cutoff <- 0.5

#Filter genes based on MRE values
df_genes <- df_genes %>%
  dplyr::rename_with(
    ~paste0("mre", .x), !starts_with("gene")
  ) %>%
  dplyr::rename(gene_name = gene_names_from_counts)

  df_genes_selection <- df_genes %>%
    dplyr::filter(across(-1, ~ . < mre_cutoff)) 

df_genes_long <- df_genes %>%
  pivot_longer(
    cols = starts_with("mre"),
    names_to = "sample",
    values_to = "mre_score"
  )


# Calculate max density from the data
density_data <- density(df_genes_long$mre_score, na.rm = TRUE)
max_density <- max(density_data$y)

ggplot(df_genes_long, aes(x = mre_score, fill = sample)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = mre_cutoff, linetype = "dashed") +
  annotate(
    "text",
    x = mre_cutoff,
    y = max_density * 0.9,  
    label = paste("MRE Threshold =", round(mre_cutoff, 2)),
    color = "#FF6600",
    hjust = -0.1,
    size = 4.5
  ) +
  # Aesthetics
  scale_fill_viridis_d() + 
  scale_color_viridis_d() + 
  theme_minimal(base_size = 14) +
  labs(
    title = "MRE Score Distributions Across Samples",
    x = "MRE Score",
    y = "Density",
    fill = "Sample",  # Legend titles
    color = "Sample"
  )


#Keep only the introns of selected genes
df_introns <- df_introns %>%
  dplyr::filter(gene_name %in% df_genes_selection$gene_name)

df_bins <- df_bins %>%
  dplyr::filter(gene_name %in% df_genes_selection$gene_name)

#Linear regression for each intron and calculation of slope and R^2

results <- purrr::map_dfr(df_introns$intron_id, function(v) {
  v_intron <- which(df_introns$intron_id == v)
  v_gene <- which(df_bins$intron_id == v)
  
  length_of_intron <- as.numeric(df_introns[v_intron,"intron_lengths"])
  x<-seq_along(v_gene)
  
  purrr::map_dfr(count_cols, function(col) {
  
  q<-unlist(df_bins[v_gene, col])
  rlm_intron <- MASS::rlm(q ~ x)
  slope <- coef(rlm_intron)[[2]]
  residual <- rlm_intron$resid
  fitted_value <- rlm_intron$fitted.values
  total_ss <- sum((q-mean(q))^2)
  residual_ss <- sum(residual^2)
  R_squared <- 1 - (residual_ss / total_ss)
  
  data.frame(
    intron_id = v,
    sample = sub("^trimmed_counts_", "", col),
    slope = slope,
    R_squared = R_squared,
    length = length_of_intron
    )
  })
})

df_introns_with_results <- left_join(
  df_introns,
  results,
  by = "intron_id"
)

df_introns_with_results <- df_introns_with_results %>%
  dplyr::mutate(
    slope = ifelse(strand == "-", slope * -1 ,slope)) %>%
  dplyr::select(-matches("counts"),
                -intron_lengths)

#Set cutoff values for length and R^2
R_squared_cutoff <- 0.5
length_cutoff <- 1000

df_introns_selection <- df_introns_with_results %>%
  pivot_wider(
    id_cols = c(chromosome,gene_name,intron_rank,intron_id,
                intron_begin,intron_finish,strand,length
    ),
    names_from = sample,
    values_from = c(slope, R_squared)
)
df_introns_selection <- df_introns_selection %>%
  dplyr::filter(
    length > length_cutoff,
    if_all(starts_with("slope"), ~ .x < 0),
    if_all(starts_with("R_squared"), ~ .x >= R_squared_cutoff)
  ) %>%
  dplyr::mutate(
    across(starts_with("slope"), sd, .names = "stn_dev_{.col}"),
    across(starts_with("slope"), mean, .names = "mean_{.col}"),
    )

df_introns_selection <- df_introns_selection %>%
  rowwise() %>%
  filter(
    any(
      c_across(starts_with("slope")) > 
        c_across(starts_with("mean_slope")) - 3 * c_across(starts_with("stn_dev_slope")) &
        c_across(starts_with("slope")) < 
        c_across(starts_with("mean_slope")) + 3 * c_across(starts_with("stn_dev_slope"))
    )
  ) %>%
  ungroup()

df_introns_selection_long <- df_introns_selection %>%
  pivot_longer(
    cols = starts_with("slope"),
    names_to = "sample",
    names_prefix = "slope_",
    values_to ="slope"
  )

# Compute mode per sample using density
get_mode <- function(v) {
  d <- density(v, adjust = 1.5, na.rm = TRUE)
  d$x[which.max(d$y)]
}


#Calculate modes of slopes for plotting
modes <- df_introns_selection_long %>%
  group_by(sample) %>%
  summarise(
    mode_val = get_mode(slope),
    max_dens = max(density(slope, adjust = 1.5, na.rm = TRUE)$y),
    .groups = "drop"
  )

# Plot
ggplot(df_introns_selection_long, aes(x = slope, color = sample, fill = sample)) +
  geom_density(alpha = 0.3, adjust = 1.5, linewidth = 0.8) +
  
  # Vertical lines for modes
  geom_vline(
    data = modes,
    aes(xintercept = mode_val, color = sample),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  # Mode labels 
  geom_label(
    data = modes,
    aes(
      x = mode_val,
      y = max_dens * 0.9,
      label = sprintf("Mode_slope: %.2f", mode_val),  # Changed here
      color = sample
    ),
    fill = "white",
    alpha = 0.8,
    size = 4,
    show.legend = FALSE,
    hjust = 0.5,
    vjust = 0,
    label.size = 0.5,
    label.padding = unit(0.2, "lines")
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold")
  ) +
  
  scale_color_viridis_d(aesthetics = c("color", "fill"), end = 0.9) +
  labs(
    title = "Distribution of Intron Slopes by Sample",
    x = "Slope Value",
    y = "Density",
    color = "Sample",
    fill = "Sample"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_rug(alpha = 0.1, sides = "b", length = unit(0.02, "npc"))


dev.off()