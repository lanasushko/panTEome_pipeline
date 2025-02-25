# Load the optparse library
library(optparse)

# Define the options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input file path", metavar = "character")
)

# Parse the options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input_file <- opt$input

data=read.table(input_file, header=T)

# subset unclassified
data_unclassified=subset(data,data$cl_uncl=='unclassified')

# subset classified
data_classified=subset(data,data$cl_uncl=='classified')


### GET LENGTHS FROM CLASSIFIED TRANSPOSONS ###
# Here we are getting the 10-90% percentile of the average classified sequence length

data=data_classified

##### DNA TIR
tir.left=quantile(subset(data, data$order=='DNA' & data$superfamily!='Helitron')$length, probs=c(0.10))[[1]]
tir.right=quantile(subset(data, data$order=='DNA' & data$superfamily!='Helitron')$length, probs=c(0.90))[[1]]
##### DNA Helitron
helitron.left=quantile(subset(data, data$order=='DNA' & data$superfamily=='Helitron')$length, probs=c(0.10))[[1]]
helitron.right=quantile(subset(data, data$order=='DNA' & data$superfamily=='Helitron')$length, probs=c(0.90))[[1]]
##### LINE
line.left=quantile(subset(data, data$order=='LINE')$length, probs=c(0.10))[[1]]
line.right=quantile(subset(data, data$order=='LINE')$length, probs=c(0.90))[[1]]
##### LTR
ltr.left=quantile(subset(data, data$order=='LTR')$length, probs=c(0.10))[[1]]
ltr.right=quantile(subset(data, data$order=='LTR')$length, probs=c(0.90))[[1]]



### FILTERING UNCLASSIFIED FAMILIES ###
# Here we are filtering the unclassified families according to the lengths of protein-containing families or literature-based thresholds (for MITEs and SINEs)

data=data_unclassified


# LTR
print('LTR')
quantile(subset(data, data$order=='LTR')$length, probs=c(0.20,0.80))
LTR_filtered_cst <- subset(data,data$order=='LTR' & data$length > ltr.left & data$length < ltr.right)

# DNA helitron
print('Helitron')
quantile(subset(data, data$order=='DNA' & data$superfamily=='Helitron')$length, probs=c(0.20,0.95))
Helitron_filtered_cst <- subset(data,data$order=='DNA' & data$superfamily=='Helitron' & data$length > helitron.left & data$length < helitron.right)


# TIR rest
print('TIR')
quantile(subset(data, data$order=='DNA' & data$superfamily!='Helitron')$length, probs=c(0.20,0.95))
TIR_filtered_cst <- subset(data,data$order=='DNA' & data$superfamily!='Helitron' & data$length > tir.left & data$length < tir.right)

# MITE
print('MITE')
quantile(subset(data, data$order=='MITE')$length, probs=c(0.20,0.80))
q.left=quantile(subset(data, data$order=='MITE')$length, probs=0.20)
q.right=quantile(subset(data, data$order=='MITE')$length, probs=0.80)

MITE_filtered <- subset(data,data$order=='MITE' & data$length > q.left & data$length < q.right)

# LINE
print('LINE')
quantile(subset(data, data$order=='LINE')$length, probs=c(0.20,0.80))
LINE_filtered_cst <- subset(data,data$order=='LINE' & data$length > line.left & data$length < line.right)

# SINE
print('SINE')
quantile(subset(data, data$order=='SINE' | data$order=='SINE?')$length, probs=c(0.20,0.80))
q.left=quantile(subset(data, data$order=='SINE' | data$order=='SINE?')$length, probs=0.20)
q.right=quantile(subset(data, data$order=='SINE' | data$order=='SINE?')$length, probs=0.80)

SINE=subset(data,data$order=='SINE' | data$order=='SINE?')
SINE_filtered <- subset(SINE, SINE$length > q.left & SINE$length < q.right)

# Gethering classified and filtered unclassified data
filtered_TEs=rbind(data_classified,LTR_filtered_cst,Helitron_filtered_cst,TIR_filtered_cst,MITE_filtered,LINE_filtered_cst,SINE_filtered)

# writing output
write.table(filtered_TEs, paste0(input_file,'.filtered_TEs.txt'))

