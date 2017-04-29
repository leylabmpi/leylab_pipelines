# Shiny server
library(shiny)
library(dplyr)
#library(tidyr)
#library(ggplot2)
library(propagate)
library(rhandsontable)
library(readxl)


#-- functions --#
prop_total_GB = function(N_lanes, GB_per_run, #GB_per_run_sd, 
                         Lanes_per_run, alpha=0.05){
  # n_runs = numeric
  # GB_per_run = vector(numeric...)
  df = cbind(N_lanes, Lanes_per_run, GB_per_run) 
  ex = expression(GB_per_run / Lanes_per_run * N_lanes)
  ret = propagate(ex, df, type='stat', alpha=alpha)
  ret = ret$prop[c(1, 3)] %>% as.data.frame
  colnames(ret) = c('Total_GB')
  return(ret)
}

prop_cost_per_GB = function(Total_cost, Total_GB, alpha=0.05){
  # CALC: cost_per_GB = total_cost / total_GB 
  df = cbind(Total_cost, Total_GB)
  ex = expression(Total_cost / Total_GB)
  ret = propagate(ex, df, type='stat', alpha=alpha)
  ret = ret$prop[c(1, 3)] %>% as.data.frame
  colnames(ret) = ('Cost_per_GB')
  return(ret) 
}

prop_GB_per_sample = function(Total_GB, N_samples, alpha=0.05){
  #.$Total_GB_mean, .$Total_GB_sd, .$N_samples
  # n_samples = numeric
  # total_GB = prop_object
  # CALC: GB_per_sample = total_GB / n_samples
  df = cbind(Total_GB, N_samples)
  ex = expression(Total_GB / N_samples)
  ret = propagate(ex, df, type='stat', alpha=alpha)
  ret = ret$prop[c(1, 3)] %>% as.data.frame
  colnames(ret) = c('GB_per_sample')
  return(ret)  
}

prop_target_coverage = function(GB_per_sample, Target_rel_abund,
                                Target_genome_size, alpha=0.05){
  # GB_per_sample = prop_object
  # target_rel_abund = numeric (mean,sd)
  # target_genome_size = numeric (mean,sd)
  # CALC: target_coverage = GB_per_sample * (target_rel_abund / 100) / target_genome_size
  df = cbind(GB_per_sample, Target_rel_abund, Target_genome_size)
  ex = expression(GB_per_sample * (Target_rel_abund / 100) / (Target_genome_size / 1000))
  ret = propagate(ex, df, alpha=alpha)
  ret = ret$prop[c(1, 3)] %>% as.data.frame
  colnames(ret) = c('Target_coverage')
  return(ret) 
}


make_sum_table = function(input, df_seq, df_lib){
  
  # filtering df
  if(input$sequencer == 'HiSeq_3000'){
    sequencer_reagents = input$HiSeq_sequencer_reagents
  } else
  if(input$sequencer == 'MiSeq'){
    sequencer_reagents = input$MiSeq_sequencer_reagents
  } else{
    stop('Sequencer not recognized')
  }
  df_seq = df_seq %>%
      filter(Sequencer == input$sequencer,
             Seq_reagents == sequencer_reagents) 
    
  df_lib = df_lib %>%
      filter(Lib_prep_kit == input$library_prep_kit)
  df_seq = cbind(df_seq, df_lib)
  
  # calculating 
  df_seq = df_seq %>%
      mutate(# number of sequencing lanes
             N_samples = input$n_samples,
             N_multiplex = input$n_multiplex,
             N_lanes = ceiling(N_samples / N_multiplex),
             N_samples_per_lane = N_samples / N_lanes,
             N_lanes = N_lanes * input$n_runs_per_sample,
             N_seq_reagent_kits = N_lanes,
             N_lib_prep_kits = ceiling(N_samples) / Lib_prep_kit_multiplex,
             # costs
             Total_cost = N_lib_prep_kits * Lib_prep_kit_cost + 
                          N_seq_reagent_kits * Lane_cost,
             Cost_per_lane = Total_cost / N_lanes,
             Cost_per_sample = Total_cost / N_samples
             )
  df_seq = rbind(df_seq, rep(0, length(df_seq)))
  df_seq[2, 'GB_per_run'] = df_seq[1,'GB_per_run_sd']
  df_seq = df_seq %>%
      dplyr::select(-GB_per_run_sd)
  
  # error propagation
  ## Total GB (all lanes)
  df_tmp = df_seq %>%
    do(prop_total_GB(.$N_lanes, .$GB_per_run, .$Lanes_per_run)) 
  df_seq = cbind(df_seq, df_tmp)
  ## Cost per GB
  df_tmp = df_seq %>%
    do(prop_cost_per_GB(.$Total_cost, .$Total_GB))
  df_seq = cbind(df_seq, df_tmp)  
  ## GB per sample
  df_tmp = df_seq %>%
    do(prop_GB_per_sample(.$Total_GB, .$N_samples))
  df_seq = cbind(df_seq, df_tmp)
  ## Coverage of target genome
  if(input$target_genome_bool == TRUE){
    df_tmp = data.frame(Target_rel_abund = c(input$target_rel_abund, 
                                             input$target_rel_abund_sd),
                        Target_genome_size = c(input$target_genome_size,
                                               input$target_genome_size_sd))
    df_seq = cbind(df_seq, df_tmp)
    df_tmp = df_seq %>%
      do(prop_target_coverage(.$GB_per_sample, .$Target_rel_abund, .$Target_genome_size))
  } else {
    df_tmp = data.frame(Target_coverage = c(NA, NA))
  }
  df_seq = cbind(df_seq, df_tmp)

  # formatting output  
  df_seq = df_seq %>%
    dplyr::select(GB_per_run, Lanes_per_run, 
                    Lane_cost, Lib_prep_kit_cost,
                    N_lanes, N_seq_reagent_kits, N_lib_prep_kits, 
                    N_samples_per_lane,
                    Total_cost, Cost_per_lane, 
                    Cost_per_sample, Total_GB, Cost_per_GB,
                    GB_per_sample, Target_coverage) 
  df = df_seq %>% t %>% as.data.frame
  colnames(df) = c('Mean', 'SD')
  df$Variable = gsub('_', ' ', colnames(df_seq))
  df$Variable = gsub('^N ', '# of ', df$Variable)
  df = df %>%
      dplyr::select(Variable, Mean, SD) %>%
      mutate('Mean + SD' = Mean + SD,
             'Mean - SD' = Mean - SD) %>%
      dplyr::select(-SD)
  
  
  return(df)
}


#-- server --#
shinyServer(function(input, output, session) {
  
  values = reactiveValues()
  
  df_seq = read_excel('data/seq_costs.xlsx', sheet='sequencer')
  df_lib = read_excel('data/seq_costs.xlsx', sheet='lib_prep')
  
  # summary table
  output$summaryTable = renderTable(make_sum_table(input, df_seq, df_lib),
                                    digits=1)
})
