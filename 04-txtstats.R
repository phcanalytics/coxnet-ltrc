# This file simply stores estimates from the analysis so that they can be used
# in the text
txt <- list() # List for in-line text statistics
txt$propTruncatedSim <- formatC(simdata_summary_p11$prop_truncated, 
                                format = "f", digits = 2)
txt$nPatients <- formatC(nrow(data), format = "d", big.mark = ",")
txt$nDeaths <- formatC(n_deaths, format = "d", big.mark = ",")
txt$maxP <- formatC(as.integer(max_p), format = "d", big.mark = ",")
txt$corrDxYear <- formatC(cor(data$index_date_year,
                              data$entry_days_dx), 
                          format = "f", digits = 2)
txt$entryMoreOneYear <- paste0(formatC(100 * mean(data$entry_days_dx > 365), 
                                       digits = 1, format = "f"),
                               "\\%")
txt$nCovsLarge <- formatC(ncol(train_large$x), format = "d", big.mark = ",")
txt$nTrainSmall <- formatC(nrow(train_small$x), format = "d", big.mark = ",")

# Convert statistics to data frame
txtstats <- data.frame(do.call(rbind, txt))

# Output to text file to input into latex
txtstats$def <-  "\\def"
names(txtstats)[1] <- "value"
txtstats$value <- as.character(txtstats$value)
txtstats <- data.frame(def = txtstats$def, name = rownames(txtstats), value =  txtstats$value)
txtstats$output <- paste(txtstats[, 1], " ", "\\", txtstats[, 2],
                         "{", txtstats[, 3], "}", sep = "")
fileConn <- file("txtstats.txt")
writeLines(txtstats$output, fileConn)
close(fileConn)