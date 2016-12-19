# Shiny server
library(shiny)
#library(dplyr)
#library(tidyr)
#library(ggplot2)
library(propagate)


#-- data table --#
# costs
reagent_costs = list(
  'HiSeq_3000' = list('1x50_bp' = 776.73,
                      '1x75_bp' = 1155.54,
                      '2x150_bp' = 1990.76),
  'MiSeq' = list('2x250_bp'= 1090.82,
                 '2x300_bp' = 1469.78)
)
reagent_GB = list(
  'HiSeq_3000' = list('1x50_bp' = c(115, 20),    # 105-125
                      '1x75_bp' = c(350, 25),    # 325-375
                      '2x150_bp' = c(700, 50)),  # 650-750
  'MiSeq' = list('2x250_bp'= c(3.55, 0.5),         # 3.3-3.8
                 '2x300_bp' = c(14.1, 1.8))        # 13.2-15
)
library_prep_kit_costs = c('nextera_24_rxns' = 1921.19,
                           'nextera_96_rxns' = 6500,
                           'lite' = 0,    # not sure about the actual cost
                           'PCR' = 0)
library_prep_kit_multi = c('nextera_24_rxns' = 24,
                           'nextera_96_rxns' = 96,
                           'lite' = 96,
                           'PCR' = 96)


prop_total_GB = function(n_runs, GB_per_run, alpha=0.05){
  # n_runs = numeric
  # GB_per_run = vector(numeric...)
  n_runs = c(n_runs, 0)
  df = cbind(n_runs, GB_per_run)
  ex = expression(n_runs * GB_per_run)
  ret = propagate(ex, df, alpha=alpha)
  return(ret$prop)
}

prop_cost_per_GB = function(total_cost, total_GB, alpha=0.05){
  # total_cost = numeric
  # total_GB = prop_object
  # CALC: cost_per_GB = total_cost / total_GB 
  total_GB = c(total_GB['Mean.1'], total_GB['sd.1'])
  total_cost = c(total_cost, 0)
  df = cbind(total_GB, total_cost)
  ex = expression(total_cost / total_GB)
  ret = propagate(ex, df, alpha=alpha)
  return(ret$prop) 
}

prop_GB_per_sample = function(n_samples, total_GB, alpha=0.05){
  # n_samples = numeric
  # total_GB = prop_object
  # CALC: GB_per_sample = total_GB / n_samples
  total_GB = c(total_GB['Mean.1'], total_GB['sd.1'])
  n_samples = c(n_samples, 0)
  df = cbind(total_GB, n_samples)
  ex = expression(total_GB / n_samples)
  ret = propagate(ex, df, alpha=alpha)
  return(ret$prop)  
}

prop_target_coverage = function(GB_per_sample, target_rel_abund, target_genome_size, alpha=0.05){
  # GB_per_sample = prop_object
  # target_rel_abund = numeric (mean,sd)
  # target_genome_size = numeric (mean,sd)
  # CALC: target_coverage = GB_per_sample * (target_rel_abund / 100) / target_genome_size
  GB_per_sample = c(GB_per_sample['Mean.1'], GB_per_sample['sd.1'])
  df = cbind(GB_per_sample, target_rel_abund, target_genome_size)
  ex = expression(GB_per_sample * (target_rel_abund / 100) / (target_genome_size / 1000))
  ret = propagate(ex, df, alpha=alpha)
  return(ret$prop) 
}


make_sum_table = function(input){
  # getting sequencer specifics
  reagent_cost = NA
  GB_per_run = NA
  if(input$sequencer == 'HiSeq_3000'){
    reagent_cost = reagent_costs[[input$sequencer]][[input$reagents_HiSeq_3000]]
    GB_per_run = reagent_GB[[input$sequencer]][[input$reagents_HiSeq_3000]]
  } else 
  if(input$sequencer == 'MiSeq'){
    reagent_cost = reagent_costs[[input$sequencer]][[input$reagents_MiSeq]]
    GB_per_run = reagent_GB[[input$sequencer]][[input$reagents_MiSeq]]
  }
  
  # calculating 
  ## number of runs
  n_runs = ceiling(input$n_samples / input$n_multiplex)
  n_runs = n_runs * input$n_runs_per_sample
  ## number of kits 
  n_lib_prep_kits = ceiling(input$n_samples / 
                            library_prep_kit_multi[[input$library_prep_kit]])
  n_reagent_kits = n_runs
  ## total cost
  lib_cost = library_prep_kit_costs[[input$library_prep_kit]]
  total_cost = n_lib_prep_kits * lib_cost + n_reagent_kits * reagent_cost
  ## cost per run  
  cost_per_run = total_cost / n_runs
  ## cost per sample
  cost_per_sample = total_cost / input$n_samples
  ## total bp produced
  total_GB = prop_total_GB(n_runs, GB_per_run)
  ## cost per GB
  #cost_per_GB = total_cost / total_GB 
  cost_per_GB = prop_cost_per_GB(total_cost, total_GB)
  ## GB per sample
  #GB_per_sample = total_GB / input$n_samples
  GB_per_sample = prop_GB_per_sample(input$n_samples, total_GB)
  ## coverage of target genome
  target_coverage = c('Mean.1'=NA, 'sd.1'=NA)
  if(input$target_genome_bool == TRUE){
    target_rel_abund = c(input$target_rel_abund, input$target_rel_abund_sd)
    target_genome_size = c(input$target_genome_size, input$target_genome_size_sd)
    target_coverage = prop_target_coverage(GB_per_sample, 
                                           target_rel_abund, 
                                           target_genome_size)
  }
  
  # data.frame init
  cats = c('Total cost',
           'Number of runs',
           'Cost per sample',
           'Total GB',
           'Cost per GB',
           'GB per sample',
           'Coverage of target genome')
  n_cats = length(cats)
  df_sum = data.frame('Category' = cats,
                      'Mean' = c(total_cost, n_runs, cost_per_sample, 
                                 total_GB['Mean.1'], 
                                 cost_per_GB['Mean.1'], 
                                 GB_per_sample['Mean.1'],
                                 target_coverage['Mean.1']),
                      'Mean_minus_sd' = c(NA, NA, NA,
                                        total_GB['Mean.1'] - total_GB['sd.1'],
                                        cost_per_GB['Mean.1'] - cost_per_GB['sd.1'], 
                                        GB_per_sample['Mean.1'] - GB_per_sample['sd.1'],
                                        target_coverage['Mean.1'] - target_coverage['sd.1']),
                      'Mean_plus_sd' = c(NA, NA, NA,
                                        total_GB['Mean.1'] + total_GB['sd.1'],
                                        cost_per_GB['Mean.1'] + cost_per_GB['sd.1'], 
                                        GB_per_sample['Mean.1'] + GB_per_sample['sd.1'],
                                        target_coverage['Mean.1'] + target_coverage['sd.1'])
                      )  
  
  return(df_sum) 
}


#-- server --#
shinyServer(function(input, output, session) {
  
  # summary table
  output$summaryTable = renderTable(make_sum_table(input))
     
})
