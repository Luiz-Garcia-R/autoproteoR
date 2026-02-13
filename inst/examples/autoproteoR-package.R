## ---- autoproteoR-package.R ----

# Small synthetic dataset
proteins <- paste0("P", sprintf("%03d", 1:50))

 raw_data <- data.frame(
   ProteinID = proteins,
   Control_1   = rnorm(50, mean = 10000, sd = 2000),
   Control_2   = rnorm(50, mean = 10000, sd = 2000),
   Control_3   = rnorm(50, mean = 10000, sd = 2000),
   Control_4   = rnorm(50, mean = 10000, sd = 2000),
   Control_5   = rnorm(50, mean = 10000, sd = 2000),
   Treatment_1 = rnorm(50, mean = 20000, sd = 2500),
   Treatment_2 = rnorm(50, mean = 20000, sd = 2500),
   Treatment_3 = rnorm(50, mean = 20000, sd = 2500),
   Treatment_4 = rnorm(50, mean = 30000, sd = 2500),
   Treatment_5 = rnorm(50, mean = 30000, sd = 2500)
 )

 metadata <- data.frame(
   Sample = colnames(raw_data)[-1],
   Group  = rep(c("Control", "Treatment"), each = 5)
 )

 # Import data
 proteo.import(raw_data, metadata)

 # Normalize data
 proteo.normalize(raw_data, metadata, method = "none")

 # Check normalized data
 proteo.qc(normalized_data)

 # Dimension reduction
 proteo.dimred(normalized_data)

 # Volcano plot
 proteo.volcano(normalized_data)
