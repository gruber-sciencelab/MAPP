# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: create_design_table
# ARGUMENTS: sample_design_file_path
# DESCRIPTION: Reads in the design table from a given path          
# -----------------------------------------------------------------------------
create_design_table <- function(sample_design_file_path, 
                                group_col = 'group',
                                sample_id_col = 'sample')
{
  # read the table and create row names
  design_table = 
    read.table(sample_design_file_path, 
               h=TRUE, 
               as.is=T, 
               sep="\t",
               comment.char="",
               check.names=FALSE)
  rownames(design_table) = design_table[,sample_id_col]
  
  # return the design table
  return(design_table)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: create_contrast_pairs
# ARGUMENTS: design_table
# DESCRIPTION: From a given design table, the function creates a list of
#              vectors, each of which contains a control and a treatment 
#              sample name. The list names are created from the two sample
#              names.
# -----------------------------------------------------------------------------
create_contrast_pairs <- function(design_table,
                                  contrast_col_name = 'contrast',
                                  contrast_col_control_term = 'CNTRL',
                                  #pair_col = 'name',
                                  sample_id_col = 'sample')
{
  # Add the contrast information
  if (contrast_col_name %in% colnames(design_table))
  {
    # create the contrast pairs
    contrast_pairs = list()
    #TREAT_name = design_table[(design_table$contrast != contrast_col_control_term),sample_id_col][1]
    for (TREAT_name in design_table[(design_table$contrast != contrast_col_control_term),sample_id_col])
    {
      if (design_table[TREAT_name, contrast_col_name] %in% design_table[,sample_id_col])
      {
        CNTRL_name = design_table[(design_table[,sample_id_col] == design_table[TREAT_name, contrast_col_name]), sample_id_col]
        contrast_pairs[[paste(TREAT_name, 'vs', CNTRL_name, sep="_")]] = c(TREAT_name, CNTRL_name)
      } else {
        stop(paste("[ERROR] Control sample", 
                   design_table[TREAT_name,'contrast'], "specified as contrast for sample", 
                   design_table[TREAT_name,'sample_id'], "could not be found in the design table. Sample skipped from analysis!", sep=" "))
      }
    }
    return(contrast_pairs)
  } else {
    return(NULL)
  }
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: create_windows
# ARGUMENTS: region_upstream_nt
#            region_downstream_nt
#            window_size_nt
# DESCRIPTION: create_windows
# -----------------------------------------------------------------------------
create_windows <- function(region_upstream_nt = 100,
                           region_downstream_nt = 100,
                           window_size_nt = 50,
                           step_size_nt = 25)
{
  # Create the start positions
  start_positions = 
    seq(from=-1*region_upstream_nt, to=region_downstream_nt-window_size_nt, by=step_size_nt)
  
  # Create a matrix for the regions
  regions = 
    matrix(0, nrow = length(start_positions), ncol = 5)
  
  start_pos_idx = 2
  for (start_pos_idx in 1:length(start_positions))
  {
    # Calculate the end position
    pos_start = start_positions[start_pos_idx]
    pos_end = pos_start + window_size_nt
    
    # Add the start pos to the matrix
    if ((pos_start < 0)&(pos_end <= 0)) {
      regions[start_pos_idx,1] = pos_start*-1
      regions[start_pos_idx,2] = pos_end*-1
    } else if ((pos_start < 0)&(pos_end > 0)) {
      regions[start_pos_idx,1] = pos_start*-1
      regions[start_pos_idx,4] = pos_end
    } else if ((pos_start >= 0)&(pos_end > 0)) {
      regions[start_pos_idx,3] = pos_start
      regions[start_pos_idx,4] = pos_end
    } else {
      stop("ERROR: The provided numbers do not make sense to define regions.")
    }
  }
  
  # Create the column names
  colnames(regions) = c("u.s", "u.e", "d.s", "d.e", "region_id")
  
  # u0to0.d0to50   u0to0.d50to100  u125to75.d0to0   u175to125.d0to0  u25to0.d0to25  u75to25.d0to0
  # u0to0.d25to75  u100to50.d0to0  u150to100.d0to0  u200to150.d0to0  u50to0.d0to0
  
  # Convert to data frame
  regions = 
    as.data.frame(regions, stringsAsFactors = FALSE)
  
  # Convert to numberic
  regions = transform(regions, u.s = as.numeric(u.s))
  regions = transform(regions, u.e = as.numeric(u.e))
  regions = transform(regions, d.s = as.numeric(d.s))
  regions = transform(regions, d.e = as.numeric(d.e))
  
  # Create the region ids
  regions[,c("region_id")] = paste("u", regions[,c("u.s")],
                                   "to", regions[,c("u.e")], 
                                   ".",
                                   "d", regions[,c("d.s")],
                                   "to", regions[,c("d.e")],
                                   sep="")
  
  # Set the rownames
  rownames(regions) = regions[,c("region_id")]
  
  # Define the positions
  regions[,c("pos")] = start_positions
  
  # Return the matrix
  return(regions)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: correct_by_sitecounts_expected_per_chance
# ARGUMENTS: sitecount_matrix.raw_counts
#            considered_region_length
#            nt_frequency_vector
# DESCRIPTION: Removes for each motiv the number of counts that are expected
#              to be found within a region of length 
#              'considered_region_length' given the nucleotide frequencies 
#              'nt_frequency_vector'.
#              Remove motif counts expected per chance and ensure that depleted 
#              k-mers are set to 0 counts. (Hint: it would not make sense to set them to
#              a minus value, since after centering the sitecounts per exon, they would
#              again have more or less the same distribution. That is, the k-mer counts
#              would still count, even though each poly(A) site on an exon was depleted in 
#              the motif. For instance, if we find TGT at one poly(A) site but in no other
#              one, then TGT is depleted for ALL the poly(A) sites, since we expect to 
#              find a TGT all 4^3=64 nucleotides and we expect to find 1-2 TGTs per chance
#              for each poly(A) site).
# -----------------------------------------------------------------------------
correct_by_sitecounts_expected_per_chance <- function(sitecount_matrix.raw_counts, 
                                                      considered_region_length,
                                                      nt_frequency_vector)
{
  # create a new matrix that we can correct
  sitecount_matrix.background_corrected = 
    as.data.frame(sitecount_matrix.raw_counts,
                  stringsAsFactors=FALSE)
  ##kmer = "TCTC"
  for (kmer in colnames(sitecount_matrix.raw_counts))
  {
    # calculate the probability to see the current k-mer
    kmer_nt_probs = numeric()
    for (nt in strsplit(kmer, "")[[1]]) {
      kmer_nt_probs = c(kmer_nt_probs,nt_frequency_vector[nt])
    }
    prob_to_see_kmer = prod(kmer_nt_probs)
    
    # TGTA: (0.3006*0.1996*0.3006*0.3004)*97 = 0.5255453
    nr_kmers_expected_in_region = prob_to_see_kmer*(considered_region_length-nchar(kmer)+1)
    
    # remove the expected counts from each region in the sitecount matrix
    sitecount_matrix.background_corrected[,kmer] = 
      (sitecount_matrix.background_corrected[,kmer] - nr_kmers_expected_in_region)
  }
  
  # set negative counts to zero and round the sitecounts
  sitecount_matrix.background_corrected[(sitecount_matrix.background_corrected < 0)] = 0
  
  # return the background corrected counts
  return(sitecount_matrix.background_corrected)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: center.rows
# ARGUMENTS: X
# DESCRIPTION: Center the rows of a matrix X.
# -----------------------------------------------------------------------------
center.rows <- function(X) 
{
  return( t( scale( t(X), scale=FALSE) ) )
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: center.cols
# ARGUMENTS: X
# DESCRIPTION: Center the columns of a matrix X.
# -----------------------------------------------------------------------------
center.cols = function(X) 
{
  # Subtract the column-means of a matrix
  return( scale( X, scale=FALSE))
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: calculate_activity_difference
# ARGUMENTS: activity_ctrl
#            activity_treat
# DESCRIPTION: Calculate the activity difference of two given activities.
# -----------------------------------------------------------------------------
calculate_activity_difference <- function(activity_ctrl,
                                          activity_treat)
{
  activity_difference = activity_treat - activity_ctrl
  return(activity_difference)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: calculate_delta_difference
# ARGUMENTS: delta_ctrl
#            delta_treat
# DESCRIPTION: Calculate the delta difference of two given deltas.
# -----------------------------------------------------------------------------
calculate_delta_difference <- function(delta_ctrl,
                                       delta_treat)
{ 
  delta_difference = sqrt(delta_treat^2.0 + delta_ctrl^2.0)
  return(delta_difference)
} 

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: calcultate_zvalue
# ARGUMENTS: coefficient
#            error
# DESCRIPTION: Calculate a z-value for a given coefficient and error.
# -----------------------------------------------------------------------------
calcultate_zvalue <- function(coefficient,
                              error)
{ 
  zvalue = coefficient / error
  return(zvalue)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: calculate_activity_differences
# ARGUMENTS: activities
#            deltas
#            contrast_pairs
#            treat_idx
#            ctrl_idx
# DESCRIPTION: Calculate activity differences, errors and z-scores.
# -----------------------------------------------------------------------------
calculate_activity_differences <- function(activities,
                                           deltas,
                                           contrast_pairs,
                                           treat_idx,
                                           ctrl_idx)
{  
  # specify the order of kmers we want to work with
  kmers = colnames(activities)
  
  # create matrices for the activities and the deltas
  activities.mat = matrix(nrow=length(names(contrast_pairs)),ncol=length(kmers),
                          dimnames=list(names(contrast_pairs), kmers))
  deltas.mat = matrix(nrow=length(names(contrast_pairs)),ncol=length(kmers),
                      dimnames=list(names(contrast_pairs), kmers))
  zscores.mat = matrix(nrow=length(names(contrast_pairs)),ncol=length(kmers),
                       dimnames=list(names(contrast_pairs), kmers))
  
  # calculate activity differences for the current contrast
  for (contrast in names(contrast_pairs)) 
  {
    activities.mat[contrast, kmers] = 
      calculate_activity_difference(activity_ctrl=activities[contrast_pairs[[contrast]][ctrl_idx],kmers,drop=F],
                                    activity_treat=activities[contrast_pairs[[contrast]][treat_idx],kmers,drop=F])
    deltas.mat[contrast, kmers] = 
      calculate_delta_difference(delta_ctrl=deltas[contrast_pairs[[contrast]][ctrl_idx],kmers,drop=F],
                                 delta_treat=deltas[contrast_pairs[[contrast]][treat_idx],kmers,drop=F])
    zscores.mat[contrast, kmers] = 
      activities.mat[contrast, kmers] / deltas.mat[contrast, kmers]
  }
  
  # calculate the mean z-score
  zscores.mat = 
    rbind(apply(zscores.mat, 2, mean), zscores.mat)
  rownames(zscores.mat)[1] = c("mean_diff_zscore")
  
  # finally create a zorted data frame from it
  mean_diff_zscore_df = as.data.frame(t(zscores.mat))
  mean_diff_zscore_df.sorted = mean_diff_zscore_df[order(abs(mean_diff_zscore_df$mean_diff_zscore), decreasing=T), ]
  
  # write all the results into one list  
  results_list = list()
  results_list[["activity_differences"]] = t(activities.mat)
  results_list[["activity_difference_deltas"]] = t(deltas.mat)
  results_list[["activity_difference_zscores"]] = t(zscores.mat)
  results_list[["activity_difference_zscores.sorted"]] = mean_diff_zscore_df.sorted
  
  # return the results
  return(results=results_list)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: activity_profile_plots
# ARGUMENTS: motifs
#            activities.mat
#            deltas.mat
#            results_dir=NULL
#            mar=c(10,10,10,10)
#            mgp=c(6,2,0)
#            colors=NULL
#            x_labels.vec=NULL
#            y_label=NULL
#            cex_value=1.0
#            cex_axis=1.0
#            cex_lab=1.0
#            cex_main=1.0
#            lwd=1.0
#            pdf_height=10
#            pdf_width=10
#            main.text=''
#            dot_type=19
#            additional_code_to_be_executed=NULL
# DESCRIPTION: Create a plot for multiple motifs.
# -----------------------------------------------------------------------------
activity_profile_plots <- function(motifs,
                                   activities.mat,
                                   deltas.mat,
                                   results_dir=NULL,
                                   mar=c(10,10,10,10),
                                   mgp=c(6,2,0),
                                   colors=NULL,
                                   x_labels.vec=NULL,
                                   y_label=NULL,
                                   cex_value=1.0,
                                   cex_axis=1.0,
                                   cex_lab=1.0,
                                   cex_main=1.0,
                                   lwd=1.0,
                                   pdf_height=10,
                                   pdf_width=10,
                                   main.text='',
                                   dot_type=19,
                                   additional_code_to_be_executed=NULL)
{
  # ---------------------------------------------------------------------------
  # create one plot per motif
  for (motif in motifs)
  {
    # create the plots
    activity_profile_plot(motifs_to_plot.vec=motif,
                          activities.mat=activities.mat,
                          deltas.mat=deltas.mat,
                          results_dir=results_dir,
                          file_name=paste("activity_difference_profile", sep="_"),
                          mar=mar,
                          mgp=mgp,
                          colors=c("blue"),
                          x_labels.vec=NULL,
                          y_label=y_label,
                          cex_value=cex_value,
                          cex_axis=cex_axis,
                          cex_lab=cex_lab,
                          cex_main=cex_main,
                          lwd=lwd,
                          pdf_height=pdf_height,
                          pdf_width=pdf_width,
                          main.text=main.text,
                          dot_type=dot_type,
                          additional_code_to_be_executed=NULL)        
  }
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: activity_profile_plot
# ARGUMENTS: motifs_to_plot.vec
#            activities.mat
#            deltas.mat
#            results_dir=NULL
#            file_name=NULL
#            mar=c(10,10,10,10)
#            mgp=c(6,2,0)
#            colors=NULL
#            x_labels.vec=NULL
#            y_label=NULL
#            cex_value=1.0
#            cex_axis=1.0
#            cex_lab=1.0
#            cex_main=1.0
#            lwd=1.0
#            pdf_height=10
#            pdf_width=10
#            main.text=''
#            dot_type=1
#            additional_code_to_be_executed=NULL
# DESCRIPTION: Create a motif activity profile plot.
# -----------------------------------------------------------------------------
activity_profile_plot <- function(motifs_to_plot.vec,
                                  activities.mat,
                                  deltas.mat,
                                  results_dir=NULL,
                                  file_name=NULL,
                                  mar=c(10,10,10,10),
                                  mgp=c(6,2,0),
                                  colors=NULL,
                                  x_labels.vec=NULL,
                                  y_label=NULL,
                                  cex_value=1.0,
                                  cex_axis=1.0,
                                  cex_lab=1.0,
                                  cex_main=1.0,
                                  lwd=1.0,
                                  pdf_height=10,
                                  pdf_width=10,
                                  main.text='',
                                  dot_type=1,
                                  additional_code_to_be_executed=NULL)
{
  # ---------------------------------------------------------------------------
  # check if the requirements that we need to do the plots are fulfilled
  # (we cannot put more than two activity profiles into one plot)
  if (length(motifs_to_plot.vec) > 2) {
    stop(paste("This function cannot put more than two motifs into one plot! ",
               "Please choose < 3 motifs!", sep=''))
  }
  
  # ---------------------------------------------------------------------------
  # create some colors, if we did not specify them
  if (is.null(colors))
    colors = rainbow(length(motifs_to_plot.vec))
  
  # ---------------------------------------------------------------------------
  # create the plotting dirctory if it does not exist already
  # if there is given a results dir, we write the plot into
  # a file
  if (!is.null(file_name))
  {
    if (!is.null(results_dir))
    {
      # create the directory if it does not exist
      dir.create(results_dir, 
                 showWarnings = FALSE)
    }
    # split the filename
    file_name.split = unlist(strsplit(x=file_name, split='.', fixed=TRUE))
    # create the file
    if ((length(file_name.split) > 1) & (tail(file_name.split, n=1) == "png")) {
      png(paste(results_dir, '/', file_name, sep=''), 
          height=pdf_height, width=pdf_width, res=300, units="in") # we use inches, because pdf does this also per default
    } else if ((length(file_name.split) > 1) & (tail(file_name.split, n=1) == "pdf")) {
      pdf(paste(results_dir, '/', file_name, sep=''), 
          height=pdf_height, width=pdf_width, useDingbats=F)
    } else {
      pdf(paste(results_dir, '/', file_name, '.pdf', sep=''), 
          height=pdf_height, width=pdf_width, useDingbats=F)
    }
  }
  
  # ---------------------------------------------------------------------------
  # check if there are lables given
  if (is.null(x_labels.vec)) {
    labels = rownames(activities.mat)
  } else {
    if (length(x_labels.vec) == nrow(activities.mat))
    {
      labels = x_labels.vec
    } else {
      stop('[ERROR] given lables vector length differs from activities matrix length!\n')
    }
  }
  
  # create the x_axis datapoints
  x_axis_datapoints = seq(from=1,to=length(labels),by=1)
  
  # ---------------------------------------------------------------------------
  # I think this was done by Piotr
  profileHeight = 960;
  if( max(x_axis_datapoints) > 50 )
  {
    profileHeight = 120 + max(x_axis_datapoints) * 15;
  }
  width = 1000;
  width = width/100;
  profileHeight = profileHeight/100;
  
  # ---------------------------------------------------------------------------
  # define the x-axis boarders
  xmin <- min(x_axis_datapoints)
  xmax <- max(x_axis_datapoints) #*1.05
  
  # ---------------------------------------------------------------------------
  # create the plot for the first motif
  for (motif_idx in 1:min(c(length(motifs_to_plot.vec),2)))
  {
    motif = motifs_to_plot.vec[motif_idx]
    y_axis_activities = activities.mat[,motif]
    y_axis_deltas = deltas.mat[,motif]
    
    # -------------------------------------------------------------------------
    # define the y-axis boarders
    scaling_factor = 1000.0
    ymax = ceiling(max(activities.mat[,motif] + deltas.mat[,motif])*scaling_factor) / scaling_factor
    ymin = 0.0 - ymax
    if (ymin > min(activities.mat[,motif] - deltas.mat[,motif]))
    {
      ymin = floor(min(activities.mat[,motif] - deltas.mat[,motif])*scaling_factor) / scaling_factor
      ymax = 0.0 - ymin
    }  
    
    # -------------------------------------------------------------------------
    # create the plot
    if (motif_idx == 1) {
      par(mar=mar, mgp=mgp)
      side=2
      padj=1
    } else {
      par(new=T)
      side=4
      padj=-1
    }
    plot(x_axis_datapoints, 
         y_axis_activities,
         xlim = c(xmin, xmax),
         ylim = c(ymin, ymax),
         type = "o",
         lwd=lwd,
         #cex.axis=cex_axis,
         cex.axis=cex_axis,
         cex.lab=cex_lab,
         cex.main=cex_main,
         main = main.text, 
         col = colors[motif_idx],
         pch = dot_type,
         yaxt="n", ylab="", xaxt="n", xlab="", bty="n",
         axes = F);
    
    # create the error bars
    arrows(x_axis_datapoints, y_axis_activities - y_axis_deltas,
           x_axis_datapoints, y_axis_activities + y_axis_deltas,
           code = 3, angle = 90, length= 0.05,
           col = colors[motif_idx],
           lty = 1, lwd=lwd)
    
    # add the y axis (try also the adj=0.1 param)
    axis(side=side, line = 0, padj=padj, cex.axis=cex_axis)
    
    # add the y-axis text
    if (is.null(y_label)) {
      ylab_text = paste(unlist(strsplit(x=motif, split='.', fixed=TRUE))[1], "activities", sep=' ')
    } else {
      ylab_text = y_label
    }
    mtext(side=side, line = 3, col=colors[motif_idx], cex=cex_lab,
          text = ylab_text)
  }
  
  # ---------------------------------------------------------------------------
  # add the rest to the plot
  
  # add a line at 0
  abline(0, 0, lty=2, col = "black")
  
  # add the x-axis
  axis(1, 
       at=c(xmax, xmin),
       labels=FALSE,
       lwd=1,
       lwd.ticks=0)
  
  # add the tick marks for the x-axis
  axis(1, 
       at=unique(x_axis_datapoints), 
       labels=FALSE,
       tcl=-0.5,        # how long are the ticks lines and in which direction do they go
       lwd.ticks=1)
  
  # add the labels for the x-axis
  axis(1, 
       at=unique(x_axis_datapoints), #c(1.0,0.5,0.0,-0.5,-1.0,-1.5,-2.0),
       lwd=0,
       cex.axis=cex_axis,
       lwd.ticks=0,       
       labels=labels,
       las=2)
  
  # add a box to the plot
  box()
  
  # if we should add some additional staff to the plot, we do that now
  if (!is.null(additional_code_to_be_executed))
  {
    eval(parse(text=additional_code_to_be_executed))
  }
  
  # if we have created a new file, we have also to close it now
  if (!is.null(file_name))
  {
    dev.off();
  }  
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: write_incl_rownames
# ARGUMENTS: data
#            col_name
#            filename
# DESCRIPTION: Write a dataframe or matrix including its rownames.
# -----------------------------------------------------------------------------
write_incl_rownames <- function(data, col_name, filename, verbose=FALSE)
{
  # merge the rownames to the matrix
  data.incl_rownames =
    cbind(rownames(data),
          data)
  
  # create a column name
  colnames(data.incl_rownames)[1] = col_name
  
  # create the full filename
  full_filename = paste(filename, '.tsv', sep='')
  
  # write the data to the file
  write.table(data.incl_rownames, file=full_filename,
              quote=F, row.names=F, col.names=T, sep="\t")
  
  # give some feedback to the user
  if (verbose) {
    message(paste('[INFO] Wrote successfully: ', full_filename, '\n', sep=''))
  }
}
# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: fit_model
#           Fits a model using a specified model.
# ARGUMENTS: contrast_pairs
#            treat_idx
#            ctrl_idx
#            sitecount_matrix.group_centered
#            expr_table
#            PAS_overlap_col
#            polyA2exon_mapping
#            exon_col
#            center_data = TRUE
#            compare_kmer_results_to_full_model = FALSE
#            randomize_data = FALSE
#            verbose = FALSE
# DESCRIPTION: fits a linear model to the given data
# -----------------------------------------------------------------------------
fit_model <- function(sitecount_matrix,
                      expr_table,
                      samples_to_consider,
                      results_dir,
                      nr_tests_for_multiple_testing_correction = NULL,
                      contrast_pairs = NULL,
                      contrast_pairs.treat_idx = NULL,
                      contrast_pairs.ctrl_idx = NULL,
                      prepare_data_function_name = "prepare_data_for_model",
                      prepare_data_function_attr_list = NULL,
                      model_function_name = "fit_simple_linear_model_for_each_motif",
                      create_files_for_each_motif = FALSE,
                      create_plots_for_each_motif = FALSE,
                      center_data = TRUE,
                      run_permutation_test = FALSE,
                      number_of_randomized_runs = 1000,
                      store_all_details_of_random_runs = TRUE,
                      docs_80_underlines = paste(rep("_", 80), collapse = ""),
                      docs_80_dashes = paste(rep("-", 80), collapse = ""),
                      verbose = FALSE)
{
  # Source the model function file and check if the function exists
  model_function = eval(parse(text = model_function_name))
  if (!is.function(model_function)) 
    stop("Argument 'model_function' is not a function!")
  
  # Provide some feedback to the user 
  if (verbose) {
    message(paste(docs_80_underlines, sep=""))
    message(paste(docs_80_dashes, sep=""))
    message(paste("[INFO] FITTING MOTIF ACTIVITIES.", sep=""))
    message(paste(docs_80_dashes, sep=""))
  }
  
  # Create the results directory (in case it does not yet exist)
  dir.create(results_dir, showWarnings = FALSE)
  
  # -----------------------------------------------------------------------------
  # we require to perform at least 30 random runs in order to get p-values 
  # calculated based on a z-statistic (see: https://en.wikipedia.org/wiki/Z-test)
  min_zstatistic_sample_size = 30
  recommended_zstatistic_sample_size = 1000
  
  # -----------------------------------------------------------------------------
  # The threshold that will be used in order to consider a distribution to be
  # normal.
  not_norm_pval_threshold = 0.05
  
  # _____________________________________________________________________________
  # -----------------------------------------------------------------------------
  # Now fit the model to the actual data (not randomized)
  # -----------------------------------------------------------------------------
  result = 
    model_function(sitecount_matrix = sitecount_matrix, 
                   expr_table = expr_table,
                   samples_to_consider = samples_to_consider,
                   results_dir = results_dir,
                   contrast_pairs = contrast_pairs,
                   contrast_pairs.treat_idx = treat_idx,
                   contrast_pairs.ctrl_idx = ctrl_idx,
                   prepare_data_function_name = prepare_data_function_name,
                   prepare_data_function_attr_list = prepare_data_function_attr_list,
                   create_files_for_each_motif = create_files_for_each_motif,
                   create_plots_for_each_motif = create_plots_for_each_motif,
                   center_data = center_data,
                   randomize_data = FALSE,
                   verbose = verbose)
  
  # get the result matrix
  fit.result_matrix = 
    result[["results_matrix"]]
  
  # get the rownames, so that we can bind afterwards the random result to it
  fit.result_matrix = as.data.frame(fit.result_matrix)
  fit.result_matrix.rownames = rownames(fit.result_matrix)
  
  # get the criterion for ranking the final results
  # HINT: The used model knows what can be used as ranking criterion
  #       (also based on the information about having contrast pairs)
  motif_ranking_criterion = result[["motif_ranking_criterion"]]
  
  # _____________________________________________________________________________
  # -----------------------------------------------------------------------------
  # Fit the model to randomized data (if wanted) 
  # -----------------------------------------------------------------------------
  # Create a list in which the results of all the random runs will be collected
  result.RandRuns_list = list()
  
  # Create a vector for the column names of the random runs
  fit.result_matrix.rand_run_cols = character(0)
  
  # Flag for running with background correction
  fit_model_with_background_rand_correction = TRUE
  
  if (number_of_randomized_runs < min_zstatistic_sample_size) {
    
    # in case the number of randomized runs is smaller than the recommended
    # number, drop a warning.
    warning(paste("[WARNING] The chosen number of randomized runs ",
                  "'--number_of_randomized_runs' (=", 
                  number_of_randomized_runs, ") ",
                  "is smaller than the minimum number required (=",
                  min_zstatistic_sample_size, ") for calculating ",
                  "z-scores on activity differences. ",
                  "Therefore z-scores and p-values can not be provided ",
                  "in the final results. ",
                  sep=""))
    
    fit_model_with_background_rand_correction = FALSE
    
  } else {
    
    if (verbose) {
      message(paste(docs_80_underlines, sep=""))
      message(paste(docs_80_dashes, sep=""))
      message("[INFO] FITTING MOTIF ACTIVITIES FOR RANDOMIZED EXPRESSION DATA.")
      message(paste(docs_80_dashes, sep=""))
    }
    
    # in case the number of randomized runs is smaller than the recommended
    # number, drop a warning.
    if (number_of_randomized_runs < recommended_zstatistic_sample_size) {
      warning(paste("[WARNING] The chosen number of randomized runs ",
                    "'--number_of_randomized_runs' (=", 
                    number_of_randomized_runs, ") ",
                    "is relatively small and might lead to variable results. ",
                    "The recommended minimum number of randomized runs is ",
                    recommended_zstatistic_sample_size, ".",
                    sep=""))
    }
    
    # Specify the prefix for the random result runs
    rand_runs_prefix = "Random_run"
    #rand_run_nr = 1
    for (rand_run_nr in (1:number_of_randomized_runs))
    {
      message(paste("[INFO] Fitting activities for randomized expression data: ", rand_run_nr, sep=""))

      # run fit on randomized data
      fit.rand_result = 
        model_function(sitecount_matrix = sitecount_matrix, 
                       expr_table = expr_table,
                       samples_to_consider = samples_to_consider,
                       results_dir = NULL,
                       contrast_pairs = contrast_pairs,
                       contrast_pairs.treat_idx = treat_idx,
                       contrast_pairs.ctrl_idx = ctrl_idx,
                       prepare_data_function_name = prepare_data_function_name,
                       prepare_data_function_attr_list = prepare_data_function_attr_list,
                       create_files_for_each_motif = FALSE,
                       create_plots_for_each_motif = FALSE,
                       center_data = center_data,
                       randomize_data = TRUE,
                       verbose = FALSE)
      
      # Store all the results if requested
      if (store_all_details_of_random_runs) {
        result.RandRuns_list[[paste(rand_runs_prefix, rand_run_nr, sep="_")]] = 
          fit.rand_result
      }
      
      # get the results matrix  
      fit.rand_result_matrix = 
        fit.rand_result[["results_matrix"]]
      
      # bind the result to the original result
      fit.result_matrix[fit.result_matrix.rownames,
                        paste(rand_runs_prefix, rand_run_nr, motif_ranking_criterion, sep="_")] = 
        fit.rand_result_matrix[fit.result_matrix.rownames,motif_ranking_criterion, drop=FALSE]
    }
    
    # If there has not been provided a specific number of test for multiple testing
    # we have to at least correct for the number if fits done here.
    if (is.null(nr_tests_for_multiple_testing_correction)) {
      nr_tests_for_multiple_testing_correction = ncol(fit.rand_result[["Ns"]])
    }
    
    # Provide some feedback to the user 
    if (verbose) {
      message(paste(docs_80_underlines, sep=""))
      message(paste(docs_80_dashes, sep=""))
      message(paste("[INFO] Number of motifs ",
                    "used for multiple testing correction: ", 
                    nr_tests_for_multiple_testing_correction, sep=""))
    }
    
    # _____________________________________________________________________________
    # -----------------------------------------------------------------------------
    # Calculate mean, stdev, z-scores for the runs on random data
    # -----------------------------------------------------------------------------
    # get the colnames of the random run results
    fit.result_matrix.rand_run_cols = 
      grep(x=colnames(fit.result_matrix), pattern = rand_runs_prefix, value=TRUE)
    
    # calculate the mean of the random runs
    fit.result_matrix[,paste(rand_runs_prefix, "mean", sep="_")] = 
      apply(fit.result_matrix[,fit.result_matrix.rand_run_cols], 1, mean)
    
    # calculate the standard deviation of the random runs
    fit.result_matrix[,paste(rand_runs_prefix, "stdev", sep="_")] = 
      apply(fit.result_matrix[,fit.result_matrix.rand_run_cols], 1, sd)
    
    # calculate the z-score for the motif_ranking_criterion
    motif_ranking_criterion_zscore_colname = paste(motif_ranking_criterion, "ZSCORE", sep="_")
    fit.result_matrix[,motif_ranking_criterion_zscore_colname] = 
      apply(fit.result_matrix, 1, 
            function(x) {
              diff_from_mean = (x[motif_ranking_criterion] - x[paste(rand_runs_prefix, "mean", sep="_")])
              zscore = diff_from_mean / x[paste(rand_runs_prefix, "stdev", sep="_")]
              return( zscore )
            })
    
    # calculate the pval for the motif_ranking_criterion
    motif_ranking_criterion_pval_colname = paste(motif_ranking_criterion, "PVAL", sep="_")
    fit.result_matrix[,motif_ranking_criterion_pval_colname] = 
      apply(fit.result_matrix, 1, 
            function(x) {
              pval.two_sided = 2 * pnorm(-abs( x[motif_ranking_criterion_zscore_colname] ))
              return( pval.two_sided )
            })
    
    # calculate the adjusted pval 
    motif_ranking_criterion_padj_colname = paste(motif_ranking_criterion, "PADJ", sep="_")
    fit.result_matrix[,motif_ranking_criterion_padj_colname] = 
      p.adjust(p = fit.result_matrix[,motif_ranking_criterion_pval_colname], 
               method = "bonferroni", 
               n = nr_tests_for_multiple_testing_correction)
    
    # test whether the values are normally distributed
    not_norm_pval_colname = paste("rand_not_normal_PVAL", sep="_")
    fit.result_matrix[,not_norm_pval_colname] = 
      apply(fit.result_matrix[,fit.result_matrix.rand_run_cols], 1, 
            function(x) {
              return(shapiro.test(x)$p.value)
            })
    
    # calculate the adjusted pval for not normal distribution.
    not_norm_padj_colname = paste("rand_not_normal_PADJ", sep="_")
    fit.result_matrix[,not_norm_padj_colname] = 
      p.adjust(p = fit.result_matrix[,not_norm_pval_colname], 
               method = "bonferroni", 
               n = nr_tests_for_multiple_testing_correction)
    
    # _____________________________________________________________________________
    # -----------------------------------------------------------------------------
    # Check how many of the z-score distributions are not normally distributed
    # (in this case a z-statistics might not be good to use...)
    # -----------------------------------------------------------------------------
    fit.result_matrix.not_norm_distributed = 
      fit.result_matrix[(fit.result_matrix[,not_norm_padj_colname] < not_norm_pval_threshold),]
    
    # calculate how much percent of the k-mer z-scores are not normal distributed
    percent_of_not_norm_distributed_kmers = 
      (nrow(fit.result_matrix.not_norm_distributed) / nrow(fit.result_matrix))*100
    
    # give some user feedback
    if (percent_of_not_norm_distributed_kmers > 0.0)
    {
      warning(paste("[WARNING] ", percent_of_not_norm_distributed_kmers, 
                    " percent of the z-scores are not normal distributed",
                    " (using a p-value threshold of ", not_norm_pval_threshold, ").", 
                    " Reported z-scores and p-values for activity differences",
                    " are only valid for motifs (k-mers) that are",
                    " normally distributed (as indicated by the '",
                    not_norm_pval_colname, "' column in the results).",
                    sep=""))
    }
  }
  
  # _____________________________________________________________________________
  # -----------------------------------------------------------------------------
  # Write out the k-mers sorted by the motif_ranking_criterion
  # -----------------------------------------------------------------------------
  if (verbose) {
    message(paste(docs_80_underlines, sep=""))
    message(paste(docs_80_dashes, sep=""))
    message(paste("[INFO] SORTING RESULTS AND WRITING THEM TO THE RESULTS FILE.", sep=""))
  }
  
  # _________________________________________________________________________
  # -------------------------------------------------------------------------
  # Assemble the final results matrix
  # -------------------------------------------------------------------------
  # Get first the ranking criterion and the p-values (in case we have them)
  if (fit_model_with_background_rand_correction == TRUE)
  {
    final_results_matrix = 
      fit.result_matrix[,c(motif_ranking_criterion,
                           motif_ranking_criterion_padj_colname,
                           not_norm_padj_colname)]
  } else {
    final_results_matrix = 
      fit.result_matrix[,c(motif_ranking_criterion), drop=F]
  }
  # colnames(final_results_matrix)[1] =
  #   paste(colnames(final_results_matrix)[1], "RANKING", sep="_")
  
  # Define the metrics that should be add to the final results
  metrics_to_add = c("activities","deltas","zscores")
  
  # If the fit was done on contrast pairs, we have to add the results as well
  if (!is.null(contrast_pairs))
  {
    # extend the later merge by the metrices for contrasts
    metrics_to_add = 
      c("diff_activities",
        "diff_deltas", "diff_zscores", metrics_to_add)
    
    # Add the mean metrics for contrast runs
    final_results_matrix = 
      merge(x=final_results_matrix,
            y=fit.result_matrix[,c("mean_diff_activity","mean_diff_delta")],
            by.x="row.names",
            by.y="row.names",
            all=FALSE)
    rownames(final_results_matrix) = final_results_matrix$Row.names
    final_results_matrix$Row.names <- NULL
  }
  
  # Finally, add the individual activities, errors and zscores
  ##metric = "activities"
  for (metric in metrics_to_add) 
  {
    # Get the activities
    result.metric = t(result[[metric]])
    colnames(result.metric) = 
      paste(metric, colnames(result.metric), sep = "_")
    
    # Merge it to the results
    final_results_matrix = 
      merge(x=final_results_matrix,
            y=result.metric,
            by.x="row.names",
            by.y="row.names",
            all=FALSE)
    rownames(final_results_matrix) = final_results_matrix$Row.names
    final_results_matrix$Row.names <- NULL
  }
  
  # sort it also by zscore
  final_results_matrix.sorted = 
    final_results_matrix[(order(abs(final_results_matrix[,motif_ranking_criterion]),
                                decreasing=T)),]
  
  if (verbose) {
    message(paste(docs_80_underlines, sep=""))
    message(paste(docs_80_dashes, sep=""))
    message(paste("[INFO] fit run finished.", sep=""))
  }
  
  return(list(motif_ranking_criterion=motif_ranking_criterion,
              fit_results.summary_statistics=final_results_matrix.sorted,
              fit_results.actual_data=result,
              fit_results.random_runs=result.RandRuns_list))
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: fit_simple_linear_model_for_each_motif
#           Fits a simple linear model.
# ARGUMENTS: contrast_pairs
#            treat_idx
#            ctrl_idx
#            sitecount_matrix.group_centered
#            expr_table
#            PAS_overlap_col
#            polyA2exon_mapping
#            exon_col
#            center_data = TRUE
#            compare_kmer_results_to_full_model = FALSE
#            randomize_data = FALSE
#            verbose = FALSE
# DESCRIPTION: fits a linear model to the given data
# -----------------------------------------------------------------------------
fit_simple_linear_model_for_each_motif <- function(sitecount_matrix,
                                                   results_dir,
                                                   expr_table,
                                                   samples_to_consider,
                                                   contrast_pairs = NULL,
                                                   contrast_pairs.treat_idx = NULL,
                                                   contrast_pairs.ctrl_idx = NULL,
                                                   prepare_data_function_name = "prepare_data_for_model",
                                                   prepare_data_function_attr_list = NULL,
                                                   create_files_for_each_motif=FALSE,
                                                   create_plots_for_each_motif=FALSE,
                                                   center_data = TRUE,
                                                   randomize_data = FALSE,
                                                   verbose = FALSE)
{
  # ___________________________________________________________________________
  # ---------------------------------------------------------------------------
  # Source the dependencies
  # ---------------------------------------------------------------------------
  # Source the model dependencies
  
  # Source the prepare data function file and check if the function exists
  prepare_data_function = eval(parse(text = prepare_data_function_name))
  if (!is.function(prepare_data_function)) 
    stop("Argument 'prepare_data_function' is not a function!")
  
  # ___________________________________________________________________________
  # ---------------------------------------------------------------------------
  # Create the finally returned resutls list
  # ---------------------------------------------------------------------------
  # Create the results list that we can then returns
  results_list = list()
  
  # ___________________________________________________________________________
  # ---------------------------------------------------------------------------
  # Prepare the matrix for the end results
  # ---------------------------------------------------------------------------
  if (!is.null(contrast_pairs))
  {
    # Define the samples we need to consider
    samples_to_consider = unname(unlist(contrast_pairs))
    
    # Check if we have expression values for all samples we want to consider
    samples_to_consider.available = 
      (samples_to_consider %in% colnames(expr_table))
    if (!all(samples_to_consider.available, na.rm = FALSE))
    {
      # Get the missing samples
      stop(paste("ERROR: The following samples are requested for ",
                 "sample comparison, but could not be found in the ",
                 "expression table provided: ", 
                 paste(samples_to_consider[!samples_to_consider.available], 
                       collapse = ", "), sep=""))
    }
    
    # Add the motif ranking criterion
    results_list[["motif_ranking_criterion"]] = "mean_diff_zscores"
    
  } else {
    # Add the motif ranking criterion
    results_list[["motif_ranking_criterion"]] = "combined.Zscore"
  }
  
  # Prepare the data for fitting a linear model
  linear_model_fit_data = 
    prepare_data_function(samples_to_consider = samples_to_consider,
                          sitecount_matrix = sitecount_matrix,
                          expr_table = expr_table,
                          additional_attr_list = prepare_data_function_attr_list,
                          center_data = center_data,
                          randomize_data = randomize_data,
                          verbose = verbose)
  
  # Get the matrixes
  Ns = linear_model_fit_data$Ns
  Es = linear_model_fit_data$Es
  
  # ___________________________________________________________________________
  # ---------------------------------------------------------------------------
  # Prepare the matrixes for the end results
  # ---------------------------------------------------------------------------
  # Add to the results list that we can then returns
  results_list[["expression_table_not_centered"]] = 
    linear_model_fit_data$expression_table_not_centered
  results_list[["Es"]] = Es
  results_list[["Ns"]] = Ns
  
  # ___________________________________________________________________________
  # ---------------------------------------------------------------------------
  # Create matrixes for collecting all the fits
  # ---------------------------------------------------------------------------
  # Get the motifs we are interested in
  motifs = colnames(Ns)
  
  # create matrices for the activities and the deltas
  activities.mat = matrix(nrow=length(samples_to_consider),ncol=length(motifs),
                          dimnames=list(samples_to_consider, motifs))
  deltas.mat = matrix(nrow=length(samples_to_consider),ncol=length(motifs),
                      dimnames=list(samples_to_consider, motifs))
  zscores.mat = matrix(nrow=length(samples_to_consider),ncol=length(motifs),
                       dimnames=list(samples_to_consider, motifs))
  combined_zscores.mat = matrix(nrow=1,ncol=length(motifs),
                                dimnames=list("combined.Zscore", motifs))
  
  # _____________________________________________________________________________
  # -----------------------------------------------------------------------------
  # If contrast pairs are provided, we will calculate the contrast statistics
  # -----------------------------------------------------------------------------
  if (!is.null(contrast_pairs))
  {
    # ---------------------------------------------------------------------------
    # Get the names of the contrasts
    # ---------------------------------------------------------------------------
    contrasts = names(contrast_pairs)
    
    # ---------------------------------------------------------------------------
    # Collect all the fits
    # ---------------------------------------------------------------------------
    # create matrices for the diff activities and the deltas
    diff_activities.mat = matrix(nrow=length(contrasts),ncol=length(motifs),
                                 dimnames=list(contrasts, motifs))
    diff_deltas.mat = matrix(nrow=length(contrasts),ncol=length(motifs),
                             dimnames=list(contrasts, motifs))
    diff_zscores.mat = matrix(nrow=length(contrasts),ncol=length(motifs),
                              dimnames=list(contrasts, motifs))
    
    # ---------------------------------------------------------------------------
    # create a matrix for the final results
    # ---------------------------------------------------------------------------
    result_cols = c('mean_diff_zscores', 
                    'mean_diff_activity', 
                    'mean_diff_delta')
    mean_diff_results.mat =
      matrix(ncol=length(result_cols),
             nrow=length(motifs))
    rownames(mean_diff_results.mat) = motifs
    colnames(mean_diff_results.mat) = result_cols
  }
  
  # ___________________________________________________________________________
  # ---------------------------------------------------------------------------
  # Get the fit for each k-mer
  # ------------------w---------------------------------------------------------
  # motif = motifs[10]
  # motif = "u25to0.d0to25.TCTC"
  # TODO: Parallelization of this loop is possible
  ##motif = motifs[1]
  for (motif in motifs)
  {
    if (verbose) {
      message(paste("[INFO] FITTING LINEAR MODEL FOR MOTIF: ", motif, sep=""))
    }
    
    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Fit the linear model to the current k-mer
    # -------------------------------------------------------------------------
    results = fit_linear_model(N=Ns[,motif,drop=F],E=Es,lambda=0.0)
    
    # Store the results
    activities.mat[samples_to_consider, motif] = 
      results[["Ahat"]][motif, samples_to_consider]
    deltas.mat[samples_to_consider, motif] = 
      results[["AhatSE"]][motif, samples_to_consider]
    zscores.mat[samples_to_consider, motif] = 
      results[["Zscore"]][motif, samples_to_consider]
    combined_zscores.mat["combined.Zscore", motif] = 
      results[["combined.Zscore"]]
    
    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # Calculate the contrast
    # -------------------------------------------------------------------------
    if (!is.null(contrast_pairs))
    {
      # -----------------------------------------------------------------------
      # calculate activity differences as well as the corresponding deltas and 
      # z-scores
      diff_results = calculate_activity_differences(activities=t(results$Ahat),
                                                    deltas=t(results$AhatSE),
                                                    contrast_pairs=contrast_pairs,
                                                    treat_idx=contrast_pairs.treat_idx,
                                                    ctrl_idx=contrast_pairs.ctrl_idx)
      
      # _________________________________________________________________________
      # -------------------------------------------------------------------------   
      # Store the results
      # -------------------------------------------------------------------------   
      # Activity differences, etc.
      diff_activities.mat[contrasts,motif] = 
        diff_results[["activity_differences"]][motif, contrasts]
      diff_deltas.mat[contrasts,motif] = 
        diff_results[["activity_difference_deltas"]][motif, contrasts]
      diff_zscores.mat[contrasts,motif] = 
        diff_results[["activity_difference_zscores"]][motif, contrasts]
      
      # MEAN activity differences, etc.
      mean_diff_results.mat[motif,'mean_diff_zscores'] = 
        diff_results$activity_difference_zscores[motif,'mean_diff_zscore']
      mean_diff_results.mat[motif,'mean_diff_activity'] = 
        mean(diff_results$activity_differences[motif,])
      mean_diff_results.mat[motif,'mean_diff_delta'] = 
        mean(diff_results$activity_difference_deltas[motif,])
      
      # _________________________________________________________________________
      # -------------------------------------------------------------------------   
      # Create PLOTS, if wanted
      # -------------------------------------------------------------------------   
      if (create_plots_for_each_motif)
      {
        # create the results directory for the current k-mer
        path_motif_results_dir = paste(results_dir, motif, sep="/")
        dir.create(path_motif_results_dir, showWarnings=FALSE, recursive=TRUE)
        
        # create activity difference plots
        activity_profile_plots(motifs=motif,
                               activities.mat=t(diff_results$activity_differences),
                               deltas.mat=t(diff_results$activity_difference_deltas),
                               results_dir=path_motif_results_dir,
                               mar=c(15,4,4,2),
                               mgp=c(3,1,0),
                               colors=NULL,
                               x_labels.vec=NULL,
                               y_label=paste("activity difference"),
                               cex_value=1.0,
                               cex_axis=1.0,
                               cex_lab=1.0,
                               cex_main=1.0,
                               lwd=1.0,
                               pdf_height=6,
                               pdf_width=1+length(contrast_pairs),
                               main.text=paste("activity difference profile\nof motif ", motif, sep=""),
                               dot_type=1,
                               additional_code_to_be_executed=NULL)
      }
    }
    
    # ___________________________________________________________________________
    # ---------------------------------------------------------------------------  
    # Create files and/or plots if wanted
    # ---------------------------------------------------------------------------  
    if (create_files_for_each_motif | create_plots_for_each_motif)
    {
      # create the results directory for the current k-mer
      path_motif_results_dir = paste(results_dir, motif, sep="/")
      dir.create(path_motif_results_dir, showWarnings=FALSE, recursive=TRUE)
      
      # Create FILES, if wanted
      if (create_files_for_each_motif)
      {
        # write out the activities
        write_incl_rownames(data = t(results[["Ahat"]]),
                            col_name = 'sample_id', 
                            filename = paste(path_motif_results_dir, 
                                             'activities', 
                                             sep='/'))
        
        # write out the activities
        write_incl_rownames(data = t(results[["AhatSE"]]),
                            col_name = 'sample_id', 
                            filename = paste(path_motif_results_dir, 
                                             'deltas', 
                                             sep='/'))
      }
    }
  }
  
  # Add the results to the final results list
  results_list[["activities"]] = activities.mat
  results_list[["deltas"]] = deltas.mat
  results_list[["zscores"]] = zscores.mat
  results_list[["combined_zscores"]] = combined_zscores.mat
  
  # _________________________________________________________________________
  # -------------------------------------------------------------------------
  # Calculate the contrast
  # -------------------------------------------------------------------------
  if (!is.null(contrast_pairs))
  {
    # Add the results to the final results list
    results_list[["diff_activities"]] = diff_activities.mat
    results_list[["diff_deltas"]] = diff_deltas.mat
    results_list[["diff_zscores"]] = diff_zscores.mat
    # Return the MEAN diff results as result
    results_list[["results_matrix"]] = mean_diff_results.mat
  } else {
    # Return the combined z-scores matrix as result
    results_list[["results_matrix"]] = t(combined_zscores.mat)
  }
  
  # return the results
  return(results_list)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: prepare_relative_usage_data
#           Randomizes the expression table (if randomize_data = TRUE).
# DESCRIPTION: prepare relative usage data for model fitting.
# -----------------------------------------------------------------------------
prepare_relative_usage_data_for_model <- function(samples_to_consider,
                                                  sitecount_matrix,
                                                  expr_table,
                                                  additional_attr_list, # additional_attr_list = list(pas_overlap_col="PAS_overlap")
                                                  center_data = TRUE,
                                                  randomize_data = FALSE,
                                                  verbose = FALSE)
{
  # ---------------------------------------------------------------------------
  # randomize the expression table in case wanted
  if (randomize_data == TRUE)
  {
    # randomize the rownames
    rows.original = rownames(expr_table)
    rows.randomized = sample(rows.original)
    
    # get the rows in a very different order
    expr_table.for_model_fitting = expr_table[rows.randomized,]
    
    # however, give back the original names
    rownames(expr_table.for_model_fitting) = rows.original
    
    # give also back the original overlapping info
    # (so that we conserve how many exons will be filterd out
    #  and the only thing that changed are the expression values)
    expr_table.for_model_fitting[rows.original,
                                 additional_attr_list[['pas_overlap_col']] ] = 
      expr_table[rows.original, additional_attr_list[['pas_overlap_col']] ]
    
  } else {
    # use the expression table as it is (no randomization)
    expr_table.for_model_fitting = expr_table
  }
  
  # ---------------------------------------------------------------------------
  # get the expression table
  expr_table.centered.filtered = 
    create_relative_usage_table(expr_table = expr_table.for_model_fitting,
                                PAS_overlap_col = additional_attr_list[['pas_overlap_col']],
                                polyA2exon_mapping = additional_attr_list[['polyA2exon_mapping']],
                                exon_col = additional_attr_list[['exon_col']])
  
  # ---------------------------------------------------------------------------
  # get only those entries for which we also have expression values 
  # and ensure that the matrices match
  ##expr_table.contrast_pairs = 
  expr_table.samples_to_consider = 
    expr_table.centered.filtered[,(colnames(expr_table.centered.filtered) %in% samples_to_consider)]
  
  # ---------------------------------------------------------------------------
  # get the sitecounts for the poly(A) sites we consider
  sitecount_matrix.incl_zero_cols = 
    sitecount_matrix[rownames(expr_table.samples_to_consider),]
  
  # ---------------------------------------------------------------------------
  # filter out motifs having 0 counts (since we cannot fit a model for such a
  # motif)
  sitecount_matrix.curated = 
    sitecount_matrix.incl_zero_cols[,apply(sitecount_matrix.incl_zero_cols, 2, 
                                           function(x) {
                                             !(all(x == 0.0))
                                           })]
  
  # ___________________________________________________________________________
  # ---------------------------------------------------------------------------
  # Prepare the matrixes
  # ---------------------------------------------------------------------------
  # center the columns of the sitecount matrix
  Ns = sitecount_matrix.curated
  
  # column and row center the expression matrix
  Es = expr_table.samples_to_consider
  
  # center the rows, if wanted
  if (center_data) {
    Es = center.rows(Es)
  } 
  
  # _____________________________________________________________________________
  # -----------------------------------------------------------------------------
  # Now fit the model to the actual data (not randomized)
  # -----------------------------------------------------------------------------
  return(list(expression_table_not_centered=expr_table.samples_to_consider,
              Es=Es,
              Ns=Ns))
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: create_realative_usage_table
# ARGUMENTS: expr_table
#            PAS_overlap_col
#            polyA2exon_mapping
#            exon_col
# DESCRIPTION: creates the relative usage table
# -----------------------------------------------------------------------------
create_relative_usage_table <- function(expr_table,
                                        PAS_overlap_col,
                                        polyA2exon_mapping,
                                        exon_col)
{
  # remove all exons that are overlapping with other PAS
  expr_table.non_overlapping = 
    keep_only_non_overlapping_multipas_exons(expr_table=expr_table,
                                             polyA2exon_mapping=polyA2exon_mapping,
                                             exon_col=exon_col,
                                             PAS_overlap_col=PAS_overlap_col)
  
  # add a pseudocounts
  expr_table.incl_pseudocount = add_pseudocount(expr_table=expr_table.non_overlapping,
                                                cols_to_remove=PAS_overlap_col,
                                                pseudocount=1)
  # calculate the relative usage
  rel_usage_table = calculate_relative_usage(expr_table=expr_table.incl_pseudocount,
                                             cols_to_remove=NULL,
                                             polyA2exon_mapping=polyA2exon_mapping,
                                             exon_col=exon_col,
                                             in_percent=TRUE)
  
  # add a pseudocount and then go into log space
  expr_table.log2 = get_log_expression(expr_table=rel_usage_table,
                                       cols_to_remove=NULL,
                                       pseudocount=0.0)
  
  # center the matrix columns by groups (=exons in our case)
  expr_table.centered = 
    center_matrix_rows_by_groups(matrix_to_center=expr_table.log2, 
                                 rowname_to_group_mappings=polyA2exon_mapping,
                                 group_col=exon_col)
  
  # return the expression table
  return(expr_table.centered)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: get_log_expression
# ARGUMENTS: expr_table
#            cols_to_remove=NULL
#            pseudocount=1
# DESCRIPTION: Get the log2 expression of theexpression table except the
#              'cols_to_remove'. Add a pseudocount before going to log-space.
# -----------------------------------------------------------------------------
keep_only_non_overlapping_multipas_exons  <- function(expr_table,
                                                      polyA2exon_mapping,
                                                      exon_col,
                                                      PAS_overlap_col)
{
  expr_table.incl_groups = 
    merge(x=expr_table, 
          y=polyA2exon_mapping[,exon_col,drop=F],
          by.x="row.names",
          by.y="row.names",
          all=FALSE)
  rownames(expr_table.incl_groups) = expr_table.incl_groups$Row.names
  factor_cols = sapply(expr_table.incl_groups, is.factor)
  expr_table.incl_groups[factor_cols] = 
    lapply(expr_table.incl_groups[factor_cols], as.character)
  expr_table.incl_groups.filtered = 
    do.call(rbind,
            lapply(
              split(expr_table.incl_groups, expr_table.incl_groups[,exon_col], drop=F),
              function(group_df) { 
                if (dim(group_df)[1] > 1) {
                  if (length(unique(group_df[,PAS_overlap_col])) == 1) {
                    if (unique(group_df[,PAS_overlap_col]) == "OK") {
                      return(group_df)
                    }
                  }
                }
              }))
  rownames(expr_table.incl_groups.filtered) = expr_table.incl_groups.filtered$Row.names
  expr_table.incl_groups.filtered$Row.names <- NULL
  return(expr_table.incl_groups.filtered[,!(colnames(expr_table.incl_groups.filtered) %in% exon_col)])
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: add_pseudocount
# ARGUMENTS: expr_table
#            cols_to_remove=NULL
#            pseudocount=1
# DESCRIPTION: adds a pseudocount to a dataframe or matrix within the 
#              columns that are not in 'cols_to_remove'.
# -----------------------------------------------------------------------------
add_pseudocount <- function(expr_table,
                            cols_to_remove=NULL,
                            pseudocount=1)
{
  # add the pseudocount to the columns we are interested in
  expr_table.incl_pseudo_count = 
    expr_table[,!(colnames(expr_table) %in% cols_to_remove)] + pseudocount
  
  # return the expression table
  return(expr_table.incl_pseudo_count)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: calculate_relative_usage
# ARGUMENTS: expr_table
#            cols_to_remove
#            polyA2exon_mapping
#            exon_col
#            in_percent=TRUE
# DESCRIPTION: Given a poly(A) site expression table, the function calculates
#              the relative usage for all columns except 'cols_to_remove' 
#              in percent (if in_percent=TRUE) or as fractions otherwise.
#              The usage is relative to other poly(A) sites located on the 
#              same exon as specified by the matrix 'polyA2exon_mapping', 
#              which contains poly(A) site ids as rownames and associated
#              exon ids in the 'exon_col' column.
# -----------------------------------------------------------------------------
calculate_relative_usage  <- function(expr_table,
                                      cols_to_remove,
                                      polyA2exon_mapping,
                                      exon_col,
                                      in_percent=TRUE)
{
  # map the rownames to the groups
  expr_table.incl_groups = 
    merge(x=expr_table[,!(colnames(expr_table) %in% cols_to_remove)], 
          y=polyA2exon_mapping[,exon_col,drop=F],
          by.x="row.names",
          by.y="row.names",
          all=FALSE)
  rownames(expr_table.incl_groups) = expr_table.incl_groups$Row.names
  expr_table.incl_groups$Row.names <- NULL
  
  # convert the factor back to character
  factor_cols = sapply(expr_table.incl_groups, is.factor)
  expr_table.incl_groups[factor_cols] = 
    lapply(expr_table.incl_groups[factor_cols], as.character)
  
  # center all read data by group
  relative_usage = 
    do.call(rbind,
            lapply(
              split(expr_table.incl_groups, expr_table.incl_groups[,exon_col], drop=F),
              function(group_df) { 
                return(as.matrix(sweep(group_df[,!(colnames(group_df) %in% exon_col)], 
                                       2, 
                                       colSums(group_df[,!(colnames(group_df) %in% exon_col)]), FUN='/'))) } ))
  
  # calculate in percent if wanted
  if (in_percent)
  {
    relative_usage = relative_usage * 100
  }
  
  # return the expression table
  return(relative_usage)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: get_log_expression
# ARGUMENTS: expr_table
#            cols_to_remove=NULL
#            pseudocount=1
# DESCRIPTION: Get the log2 expression of theexpression table except the
#              'cols_to_remove'. Add a pseudocount before going to log-space.
# -----------------------------------------------------------------------------
get_log_expression <- function(expr_table,
                               cols_to_remove=NULL,
                               pseudocount=1)
{
  # add the pseudocount to the columns we are interested in
  expr_table.incl_pseudo_count = 
    expr_table[,!(colnames(expr_table) %in% cols_to_remove)] + pseudocount
  
  # get the log2 values
  expr_table.log2 = log2(expr_table.incl_pseudo_count)
  
  # return the expression table
  return(expr_table.log2)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: center_matrix_rows_by_groups
# ARGUMENTS: matrix_to_center
#            rowname_to_group_mappings
#            group_col
#            scale_by_sd
# DESCRIPTION: centers the rows of a matrix ('matrix_to_center') according
#              to the group to which they belong to as specified by 
#              another matrix ('rowname_to_group_mappings') that has the 
#              same rownames as 'matrix_to_center' and the group to which 
#              the rownames belong to specified in a column with the 
#              colum nname specified by 'group_col'. The results can be
#              scaled by the standard deviation if wanted ('scale_by_sd').
# -----------------------------------------------------------------------------
center_matrix_rows_by_groups <- function(matrix_to_center, 
                                         rowname_to_group_mappings,
                                         group_col,
                                         scale_by_sd=FALSE)
{  
  # map the rownames to the groups
  matrix_to_center.incl_groups = 
    merge(x=matrix_to_center, 
          y=rowname_to_group_mappings[,group_col,drop=F],
          by.x="row.names",
          by.y="row.names",
          all=FALSE)
  rownames(matrix_to_center.incl_groups) = matrix_to_center.incl_groups$Row.names
  matrix_to_center.incl_groups$Row.names <- NULL
  
  # convert the factor back to character
  factor_cols = sapply(matrix_to_center.incl_groups, is.factor)
  matrix_to_center.incl_groups[factor_cols] = 
    lapply(matrix_to_center.incl_groups[factor_cols], as.character)
  
  # center all read data by group
  kmers.aggregated = 
    do.call(rbind,
            lapply(
              split(matrix_to_center.incl_groups, matrix_to_center.incl_groups[,group_col], drop=F),
              function(group_df) { return(scale(group_df[,!(colnames(group_df) %in% group_col), drop=F], center = TRUE, scale=scale_by_sd)) } ))
  
  # return the group centered matrix
  return(kmers.aggregated)
}
# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: fit_linear_model
# ARGUMENTS: N
#            E
#            lambda
# DESCRIPTION: Fit a linear model and return a list
#              containing a matrix with coefficients ('Ahat') 
#              and errors ('AhatSE').
# -----------------------------------------------------------------------------
fit_linear_model <- function(N,E,lambda)
{
  
  # check whether the matrixes fit
  if ( (!identical(rownames(N),rownames(E))) 
       | (is.null(rownames(N))) 
       | (is.null(rownames(E)))) {
    stop('Rownames of N and E matrixes are not equal. Please ensure that the matrices match!')
  }
  
  # do a SVD on matrix N
  Ns = perform_svd_fast(N)
  
  # calculate the right hand side
  rhs = crossprod(Ns$u,E)
  
  # calculate diagonal matrix
  dia = Ns$d/(Ns$d^2 + nrow(N)*lambda)
  
  # calculate the activities
  Ahat = sweep(Ns$v,2,dia,FUN='*') %*% rhs
  dimnames(Ahat) = list(colnames(N), colnames(E))
  
  # calculate the errors
  C = tcrossprod(sweep(Ns$v, 2, 1 / (Ns$d^2 + nrow(N)*lambda), FUN='*'), Ns$v)
  AhatSE = sqrt( diag(C) %x% t( colSums((E-N %*% Ahat)^2) / nrow(E)) )
  colnames(AhatSE) = colnames(Ahat)
  rownames(AhatSE) = rownames(Ahat)
  
  # calculate the z-scores
  Zscore = Ahat / AhatSE
  
  # calculate the combined z-score
  combined.Zscore = sqrt(rowMeans(Zscore^2))
  
  # return the results in form of a list
  fit = list(Ahat=Ahat,
             Zscore=Zscore,
             AhatSE=AhatSE,
             combined.Zscore=combined.Zscore)
  return(fit)
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: optimize.lambda
# ARGUMENTS: N
#            E
#            lambda
# DESCRIPTION: Optimize Lambda for a linear model.
# -----------------------------------------------------------------------------
optimize.lambda <- function(N,E) 
{
  # optimize lambda by generalized cross-validation
  # The input arguments N and E must be already centered
  
  Ns = perform_svd_fast(M=N) # a list with entries: u, d, v  
  rhs = crossprod(Ns$u,E)
  
  lambda.bnd = 10^c(-12,-6) * nrow(N) * ncol(N)
  
  gcv.error <- function(lambda,E,Ns,rhs) {
    D = Ns$d^2/(Ns$d^2 + nrow(N)*lambda) # Hat matrix: H = VDV^t
    resid = E - sweep(Ns$u,2,D,FUN='*') %*% rhs
    GCV = sum((resid/(nrow(E)-sum(D)))^2)
    return(GCV)
  }
  
  opt = optimize(gcv.error,lambda.bnd,E,Ns,rhs)
  lambda.opt = opt$minimum
  gcv.opt = opt$objective
  
  return(list(lambda.opt = lambda.opt,
              gcv.opt = gcv.opt))
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FUNCTION: perform_svd_fast
# ARGUMENTS: M
# DESCRIPTION: Perform an SVD (fast) on matrix M and return the result.
# SOURCE: Adapted from the R 'corpcor' package.
#         The package is developed by the Strimmer Lab 
#         (http://strimmerlab.org/software/corpcor/) and 
#         distributed under the GNU General Public License. It is 
#         available at the CRAN archive (https://cran.r-project.org).
# -----------------------------------------------------------------------------
perform_svd_fast <- function(M) 
{
  if (nrow(M) > 2*ncol(M)) {
    s = svd(crossprod(M)) 
    s$d = sqrt(s$d)
    s$u = M %*% sweep(s$v,2,s$d,FUN='/')
  } else if (ncol(M) > 2*nrow(M)) {
    s = svd(tcrossprod(M)) 
    s$d = sqrt(s$d)
    s$v = sweep(crossprod(M,s$u),2,s$d,FUN='/')
  } else {
    s = svd(M)
  }
  
  return(s)
}
# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# KAPAC (K-mer activity on poly(A) site choice)
# -----------------------------------------------------------------------------
# Developed in R version: R version 3.3.1
# -----------------------------------------------------------------------------
# Please run with KAPACv2.R with Rscript --vanilla. 
# If running manually clean up and then start with a fresh run
##rm(list=ls())

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Read options
# -----------------------------------------------------------------------------
# load libraries
library(optparse)

option_list <- list(
  make_option(c("--sample_design"), action="store", type="character", help="The sample design file. If it contains a 'contrast' column, this column will be used to create contrasts based on samples with the entry 'CNTRL'. Otherwise it will run without contrasts."),
  make_option(c("--expression_matrix"), action="store", type="character", help="The poly(A) site expression matrix."),
  make_option(c("--sitecount_matrix"), action="store", type="character", help="The site count matrix."),
  make_option(c("--pas_exon_associations"), action="store", type="character", help="The file that contains the associations of poly(A) sites to the exon on which they are located."),
  make_option(c("--selected_motifs"), default="all", action="store", type="character", help="Per default, all k-mers present in the sitecount matrix (see option --sitecount_matrix) found in a high enough fraction of poly(A) sites (see option --min_kmer_abundance_fraction) will be considered in the KAPAC run. Alternatively, a file can be provided which contains the k-mers that should be considered in a column named 'motif'."),
  make_option(c("--results_dir"), action="store", type="character", help="The directory to which the results will be written."),
  make_option(c("--create_plots_for_each_motif"), default=FALSE, action="store", type="logical", help="If set TRUE, for each kmer plots will be created (NOTE: dependent on the sitecount matrix (see option --sitecount_matrix) and the selected k-mers (see option --selected_motifs), thousands of directories and files might be created)."),
  make_option(c("--create_files_for_each_motif"), default=FALSE, action="store", type="logical", help="If set TRUE, for each kmer detailed files (activities, errors, z-scores)  will be created (NOTE: dependent on the sitecount matrix (see option --sitecount_matrix) and the selected k-mers (see option --selected_motifs), thousands of directories and files might be created)."),
  make_option(c("--row_center_expression"), default=TRUE, action="store", type="logical", help="If set TRUE, the expression matrix (see option --expression_matrix) will be row-centered. That is, changes in relative usage across samples will be explained."),
  make_option(c("--treat_minus_ctrl"), default=FALSE, action="store", type="logical", help="If true, KAPAC will consider treatment minus control for the fit and control minus treatment otherwise."),
  make_option(c("--expression_pseudocount"), default=1.0, action="store", type="double", help="The pseudocount that should be add to the expression values prior going to log-space."),
  make_option(c("--consider_excess_counts_only"), default=FALSE, action="store", type="logical", help="Do only use when running in k-mer mode (not in PWM mode). If set TRUE, background correction will be performed. That is, site counts will be corrected by the number of counts that are expected to be found per chance given the considered region length (see option --considered_region_length)."),
  make_option(c("--considered_region_length"), action="store", type="double", help="The length of the region in which the k-mers have been counted (e.g. needed for background correction, see option --consider_excess_counts_only)."),
  make_option(c("--min_kmer_abundance_fraction"), default=0.01, action="store", type="double", help="The fraction of poly(A) sites that needs to contain counts for a specific k-mer in order to be considered. That is, k-mers that have counts for a smaller fraction of poly(A) sites in the sitecount matrix (see option --sitecount_matrix) will not be considered in the KAPAC run."),
  make_option(c("--report_NAs_for_kmers_below_min_kmer_abundance_fraction"), default=TRUE, action="store", type="logical", help="If set TRUE, in the final output there will be an entry for each kmer present in the sitecount matrix (see option --sitecount_matrix), even if there was not performed any model fit for it because it does e.g. not fulfill the minimum abundance fraction (see option --min_kmer_abundance_fraction), in which case NAs are reported."),
  make_option(c("--number_of_randomized_runs"), default=1000, action="store", type="double", help="The number of runs done with randomized expression to k-mer count associations in order to calculate adjusted p-values for the reported activity difference z-scores."),
  make_option(c("--verbose"), default=TRUE, action="store", type="logical", help="If set TRUE, the script will be verbose (reporting detailed infos on what is done)."))
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS]",
                           option_list = option_list, add_help_option=TRUE)
opt <- parse_args(opt_parser)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# for debugging
# -----------------------------------------------------------------------------
# do you want to run in debug mode?
debugging_mode = FALSE

# in case we want to run in debug mode, we set all variables
if (debugging_mode == TRUE)
{
  opt = list()
  opt$sample_design = "../DATA/kapac_design.tsv"
  opt$expression_matrix = "../DATA/pas_expression.tsv"
  opt$sitecount_matrix = "../DATA/kmer_counts.tsv"
  opt$pas_exon_associations = "../DATA/pas2exon.tsv"
  opt$selected_motifs = "../DATA/selected_kmers.tsv"
  #opt$selected_motifs = "all"
  opt$results_dir = "RESULTS"
  opt$create_files_for_each_motif = TRUE
  opt$create_plots_for_each_motif = TRUE
  opt$row_center_expression = TRUE
  opt$treat_minus_ctrl = FALSE
  opt$expression_pseudocount = 1
  opt$consider_excess_counts_only = TRUE
  opt$considered_region_length = 50
  opt$min_kmer_abundance_fraction = 0.1
  opt$report_NAs_for_kmers_below_min_kmer_abundance_fraction = TRUE
  opt$number_of_randomized_runs = 30
  opt$verbose = TRUE
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# variables for writing messages, warnings and errors
# -----------------------------------------------------------------------------
docs_80_underlines = paste(rep("_", 80), collapse = "")
docs_80_dashes = paste(rep("-", 80), collapse = "")

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Create the nucleotides frequency vector
# -----------------------------------------------------------------------------
# In case we want to background correct, we will use the following
# frequencies and report it to the user
nt_frequency_vector = c(0.2973, 0.1935, 0.2007, 0.3084)
names(nt_frequency_vector) = c("A","C","G","T")

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# The name of the pas overlapping column in the expression file.
pas_overlap_col = "PAS_overlap"

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Source additional functions
# -----------------------------------------------------------------------------

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# create the results dir
# -----------------------------------------------------------------------------
dir.create(opt$results_dir, showWarnings = FALSE, recursive = TRUE)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# read in the design table and create contrasts
# -----------------------------------------------------------------------------
# create the design table from the design file
design_table = create_design_table(sample_design_file_path = opt$sample_design)

# create contrast pairs
contrast_pairs = create_contrast_pairs(design_table = design_table)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Read in the PAS to exon mapping so that we are able to center by exons
# later on
# -----------------------------------------------------------------------------
polyA2exon_mapping = 
  as.matrix(read.table(opt$pas_exon_associations, 
                       h=TRUE, 
                       as.is=T, 
                       sep="\t",
                       row.names=1,
                       comment.char="",
                       stringsAsFactors=FALSE,
                       check.names=FALSE))

# Specify the column in which the exons are specified
exon_col="exon"

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Read in the expression table
# -----------------------------------------------------------------------------
expr_table = 
  read.table(opt$expression_matrix, 
             h=TRUE, 
             as.is=T, 
             sep="\t",
             row.names=1,
             comment.char="",
             check.names=FALSE)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Read in the full sitecount matrix
# -----------------------------------------------------------------------------
sitecount_matrix.all = 
  as.matrix(read.table(opt$sitecount_matrix,
                       h=TRUE, 
                       as.is=T, 
                       sep="\t",
                       row.names=1,
                       comment.char="",
                       check.names=FALSE))

# -----------------------------------------------------------------------------
# Determine how many motifs are considered in total (can be used for 
# Bonferroni correction). We are conservative here and take the total number
# prior any filtering and selection (which makes also comparissons easier).
# -----------------------------------------------------------------------------
total_nr_of_kmers = ncol(sitecount_matrix.all)
total_set_of_kmers = colnames(sitecount_matrix.all)

# -----------------------------------------------------------------------------
# Select specific k-mers, if wanted
if (opt$selected_motifs != "all")
{
  # read in the selected kmers
  selected_kmers = 
    read.table(opt$selected_motifs, 
               h=TRUE, 
               as.is=T, 
               sep="\t",
               comment.char="",
               check.names=FALSE)[,1]
  
  # drop all other k-mers
  sitecount_matrix.all = sitecount_matrix.all[,selected_kmers,drop=F]
}

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Filter out k-mers that have counts in less than x [%] of all poly(A) sites
# in the genome.
# In case we want to run on selected k-mers, we want to get estimated
# -----------------------------------------------------------------------------
cols_to_keep = 
  apply(sitecount_matrix.all, 2, 
        function(x) { (sum(x > 0) / length(x)) >= opt$min_kmer_abundance_fraction } )

# -----------------------------------------------------------------------------
# get the columns we are interested in
sitecount_matrix.raw_counts = 
  sitecount_matrix.all[, cols_to_keep, drop = FALSE]

# -----------------------------------------------------------------------------
# get the k-mers that are filtered out
filtered_out_kmers = 
  colnames(sitecount_matrix.all[ ,!(colnames(sitecount_matrix.all) %in% colnames(sitecount_matrix.raw_counts))])

# -----------------------------------------------------------------------------
# calculate the fraction of filtered out k-mers
fraction_of_filtered_out_kmers = 
  length(filtered_out_kmers) / ncol(sitecount_matrix.all)

if (opt$verbose) {
  message(paste(docs_80_underlines, sep=""))
  message(paste(docs_80_dashes, sep=""))
  message(paste("[INFO] FILTERING OUT MOTIFS / K-MERS HAVING NO/LOW ABUNDANCE.", sep=""))
  
  # report how many k-mers will be dropped
  message(paste("[INFO] Minimum motif / k-mer abundance fraction to be considered (--min_kmer_abundance_fraction): ", 
                opt$min_kmer_abundance_fraction, "\n", sep=""))
  
  # report how many k-mers will be dropped
  message(paste("[INFO] Percentage of filtered out motifs / k-mers: ", 
                (fraction_of_filtered_out_kmers*100), "[%] (=", length(filtered_out_kmers), " motifs / k-mers). ", sep=""))
}

# background correct, if wanted
if (opt$consider_excess_counts_only == TRUE) {
  
  if (opt$verbose) {
    message(paste(docs_80_underlines, sep=""))
    message(paste(docs_80_dashes, sep=""))
    message(paste("[INFO] BACKGROUND CORRECTION: 'Excess' motif counts will be removed.", sep=""))
    message(paste("[INFO] The following nucleotide frequencies will be used:", sep=""))
    message(paste("[INFO] ", names(nt_frequency_vector), " ", nt_frequency_vector, "\n", sep=" "))
  }
  
  # perform the background correction
  sitecounts_fit = 
    correct_by_sitecounts_expected_per_chance(sitecount_matrix.raw_counts=sitecount_matrix.raw_counts, 
                                              considered_region_length=opt$considered_region_length,
                                              nt_frequency_vector=nt_frequency_vector)
} else {
  
  # give some feedback to the user
  if (opt$verbose) {
    message(paste(docs_80_underlines, sep=""))
    message(paste(docs_80_dashes, sep=""))
    message(paste("[INFO] NO BACKGROUND CORRECTION will be applied.", sep=""))
    message(paste("[INFO] Raw motif counts will be used.", sep=""))
  }
  
  # in case we do not want to perform background correction, we use the raw counts
  sitecounts_fit = sitecount_matrix.raw_counts
}

# -----------------------------------------------------------------------------
# center the motif counts per exon
# -----------------------------------------------------------------------------
if (opt$verbose) {
  message(paste(docs_80_underlines, sep=""))
  message(paste(docs_80_dashes, sep=""))
  message(paste("[INFO] CENTERING MOTIF COUNTS PER EXON.", sep=""))
}

# center the matrix
sitecount_matrix.group_centered = 
  center_matrix_rows_by_groups(matrix_to_center=sitecounts_fit, 
                               rowname_to_group_mappings=polyA2exon_mapping,
                               group_col=exon_col)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# FIT MODEL
# -----------------------------------------------------------------------------
# determine if we should plot treatment minus control or
# control minus treatment
# -----------------------------------------------------------------------------
if (opt$treat_minus_ctrl == TRUE) {
  treat_idx = 1
  ctrl_idx = 2
} else {
  treat_idx = 2
  ctrl_idx = 1
}

# -----------------------------------------------------------------------------
# fit the model
# -----------------------------------------------------------------------------

fit_results = 
  fit_model(sitecount_matrix = sitecount_matrix.group_centered,
            expr_table = expr_table,
            samples_to_consider = rownames(design_table),
            results_dir = opt$results_dir,
            nr_tests_for_multiple_testing_correction = total_nr_of_kmers,
            contrast_pairs = contrast_pairs,
            contrast_pairs.treat_idx = treat_idx,
            contrast_pairs.ctrl_idx = ctrl_idx,
            prepare_data_function_name = "prepare_relative_usage_data_for_model",
            prepare_data_function_attr_list = list(pas_overlap_col="PAS_overlap", 
                                                   polyA2exon_mapping=polyA2exon_mapping,
                                                   exon_col=exon_col),
            model_function_name = "fit_simple_linear_model_for_each_motif",
            create_files_for_each_motif = opt$create_files_for_each_motif,
            create_plots_for_each_motif = opt$create_plots_for_each_motif,
            center_data = TRUE,
            run_permutation_test = TRUE,
            number_of_randomized_runs = opt$number_of_randomized_runs,
            store_all_details_of_random_runs = FALSE,
            docs_80_underlines = paste(rep("_", 80), collapse = ""),
            docs_80_dashes = paste(rep("-", 80), collapse = ""),
            verbose = opt$verbose)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Write out the results
# -----------------------------------------------------------------------------
# Create a dataframe for the final results
final_results_df = fit_results[['fit_results.summary_statistics']]

# -----------------------------------------------------------------------------
# In case requested by the user, determine the motifs that have been excluded 
# from the fit (e.g. because the min. abundant fraction is too small)
# and add NA rows to the output for these motifs.
if (opt$report_NAs_for_kmers_below_min_kmer_abundance_fraction)
{
  # Determine the kmers that have been excluded from the final results
  kmers_excluded_from_final_results = 
    setdiff(total_set_of_kmers, rownames(final_results_df))
  
  # Create the final rownames
  final_results_df.rownames = c(rownames(final_results_df),kmers_excluded_from_final_results)
  
  # Add NA rows to the end of the final resutls data frame
  final_results_df = final_results_df[final_results_df.rownames, ,drop=F]
  rownames(final_results_df) = final_results_df.rownames
}

# -----------------------------------------------------------------------------
# Add the fraction of each k-mer to the results
regions_present_fraction = 
  apply(sitecount_matrix.all, 2, 
        function(x) { return(sum(x > 0) / length(x)) } )

# Convert to matrix
regions_present_fraction.mat = as.matrix(regions_present_fraction)
colnames(regions_present_fraction.mat)[1] = "regions_present_fraction"

# Add the fraction of motif / k-mer abundance
final_results_df = 
  cbind(final_results_df, 
        regions_present_fraction.mat[rownames(final_results_df),,drop=F])

# -----------------------------------------------------------------------------
# write out the zscores
write_incl_rownames(data=final_results_df,
                    col_name='ID', 
                    filename=paste(opt$results_dir, '/KAPAC_results',
                                   sep=''))
