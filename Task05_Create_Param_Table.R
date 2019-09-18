source("ams_initialize_script.R")
dirs$Rscript.name = "Task05_Create_Param_Table.R"
dirs$output.prefix= str_extract(dirs$Rscript.name,"^Task\\d\\d\\w?_")

#Drug list to loop through for finding file names
drugs_list = list("Pembro","VEGFR1","VEGFR2","Atezoli","Ramuc","Siltux","Tociliz") #ADD THIS LINE
source("ivsc_2cmt_RR_V1.R")
model  = ivsc_2cmt_RR_v1()

# Create paths to data files for each drug.
param = list()
i = 0
for (drug in drugs_list) { #ADD THIS LOOP
  i = i + 1
  filename = list.files(path = "parameters/",pattern = drug) #change filename line to this
  
  if (length(filename)>1)
    stop("check and see if you have any temp files open or something.  maybe close excel")
  

  param[[i]] =  file.path("./parameters/",filename) %>%
    read.param.file() %>%
    t() %>%
    as.data.frame() %>%
    mutate(drug = drug)
}

param = bind_rows(param)

param_summ = param %>%
  bind_rows() %>%
  select(-F,-ka) %>%
  select(drug,everything()) %>%
  t() %>%
  as.data.frame() 
names(param_summ) = param$drug

write.csv(param_summ, "parameters/Task05_Param_Summary.csv",quote = FALSE, row.names = TRUE)

  
