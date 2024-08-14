setwd("~/Dropbox/projects_WORKING/temporal_signal_comment/temporal_signal_comment/data_simulations/")

xml_analysis_template <- readLines("feast_template.xml")
fasta_files <- dir(pattern = "fasta$")
fasta_files

for(f in fasta_files){
  input_aln_name <- gsub(".fasta", "", f)
  xml_temp <- gsub("INPUT_ALN_NAME", input_aln_name, xml_analysis_template)
  writeLines(xml_temp, con = paste0(input_aln_name, ".xml"))
}

