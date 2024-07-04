library(dplyr)
library(readxl)

dt_raw <- read_xlsx("./data/OxCal_input_data.xlsx") %>% arrange(desc(z_middle_depth_cm))

dt <- dt_raw %>% filter(!outlier)

# function for print a line of code for each date
print_code <- function(df) {
  str1 = ifelse(df$sample_type=='14C',
         paste('R_Date("', df$sample_name_short, '",', df$R_age_yr,',',df$R_sigma_yr,')', sep=''),
         paste('Date("', df$sample_name_short, '",', 'N(calBP(', df$age_yr_BP, '),', df$uncertainty,'))', sep=''))
  str2 = paste('{Outlier(0.05);z=',df$z_middle_depth_cm,';};', sep='')
  
  return(paste(str1, str2))
}




sink("./model/oxcal_model.txt")
code <- file("./model/oxcal_model.txt")
# write the beginning
writeLines(c('Plot()',
             '{', 
             'Outlier_Model("General",T(5),U(0,4),"t");',
             'P_Sequence("",1,0.04,U(-2,2))',
             '{',
             'Boundary();'), code)
close(code)
# append codes for the dates
code <- file("./model/oxcal_model.txt", open='a')
for (i in 1:nrow(dt)) {
  cat(print_code(dt[i,]), file=code, append = T, sep = '\n')
}

close(code)

# write the ending
code <- file("./model/oxcal_model.txt", open='a')
cat('Boundary();','};','};', file=code, append = T, sep = '\n')
close(code)


