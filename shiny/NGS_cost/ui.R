# Shiny UI
library(shiny)


#-- data --#
# sequencers
sequencers = list('HiSeq 3000' = 'HiSeq_3000',
                  'MiSeq' = 'MiSeq')
HiSeq3000_reagents = list('2 x 150 bp' = '2x150_bp',
                          '1 x 75 bp' = '1x75_bp',
                          '1 x 50 bp' = '1x50_bp')
MiSeq_reagents = list('2 x 250 bp' = '2x250_bp',
                      '2 x 300 bp' = '2x300_bp')
library_prep_kits = c('LITE (Nextera-hacked)' = 'lite',
                      'Nextera 96 rxns' = 'nextera_96_rxns',
                      'Nextera 24 rxns' = 'nextera_24_rxns',
                      'Barcoded PCR (amplicon)' = 'PCR')


#-- shiny --#
# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("NGS cost (EUR)"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      h3('Set variables'),
      br(),
      # number of samples
      numericInput("n_samples",
                   label = "Number of samples",
                   value=96,
                   min=1,
                   max=NA),
      # amount of multiplexing
      numericInput("n_multiplex",
                  label = "Number of samples per run (multiplexing)",
                  value=96,
                  min=1,
                  max=NA),
      # amount of runs 
      sliderInput("n_runs_per_sample",
                  label = "Number of runs per sample (increasing coverage)",
                  value=1,
                  min=1,
                  max=10),
      # sequencer
      selectInput("sequencer", 
                  label = "Select Sequencer", 
                  choices = sequencers, 
                  selected = 1),
      # reagent kit
      conditionalPanel(
        condition = "input.sequencer == 'HiSeq_3000'",
        selectInput("reagents_HiSeq_3000", 
                    label = "Sequencing Reagents",
                    HiSeq3000_reagents)),
      conditionalPanel(
        condition = "input.sequencer == 'MiSeq'",
        selectInput("reagents_MiSeq", 
                    label = "Sequencing Reagents",
                    MiSeq_reagents)),
      # library prep
      selectInput("library_prep_kit", 
                  label = "Select library prep kit", 
                  choices = library_prep_kits, 
                  selected = 1),
      # target genome
      checkboxInput('target_genome_bool',
                    label = strong('Target genome?')),
      conditionalPanel(
        condition = "input.target_genome_bool == true",
        numericInput('target_genome_size',
                     label = 'Genome size (Mbp)',
                     value = 5.0,
                     min = 0),
        numericInput('target_genome_size_sd',
                     label = 'Stdev of genome size',
                     value = 0.1,
                     min = 0),
        numericInput('target_rel_abund',
                     label = 'Taxon relative abundance (%)',
                     value = 1.0,
                     min = 0,
                     max = 100),
        numericInput('target_rel_abund_sd',
                     label = 'Stdev of taxon rel. abund.',
                     value = 0.1,
                     min = 0,
                     max = 100))
    ),
    
    # Show a plot of the generated distribution
    mainPanel(tableOutput("summaryTable"))
  )
))