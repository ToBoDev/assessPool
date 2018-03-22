#Runs PoPoolation2

runPopoolation <- function(use_parallel, no_cores, working_dir, project_name){
  setwd(paste(working_dir, project_name, "popoolation",sep="/"))
  
  if (length(list.files(pattern=".sync")>0)){ 
    if (use_parallel){
      system(paste("parallel --jobs", no_cores, "< popool2_run.sh"))
    } else{
      system(command="sh popool2_run.sh", wait=TRUE)
    }
    message("\n\nPopoolation2 run complete.")
  } else {
    message("\n\nERROR: No .sync files found. Please run sync file setup step.")
  }
  
}