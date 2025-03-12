# Generate random data
data <- rnorm(1000, mean = 0, sd = 1)

# Plot histogram
hist(data, 
     main = "Simple Histogram", 
     xlab = "Value", 
     col = "skyblue", 
     breaks = 20)
     