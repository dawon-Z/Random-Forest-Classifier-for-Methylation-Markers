# main.R

# Load necessary libraries
library(randomForest)
library(caret)
library(ROCR)
library(pROC)

# Read data
df1 <- read.table("input.txt", header = TRUE)
df1$Probe_ID <- factor(df1$Probe_ID, levels = c(1, 2), labels = c(FALSE, TRUE))
data <- df1

sink("output.txt")

set.seed(12345)

marker_names <- colnames(data)[-1]
marker_combinations <- combn(marker_names, 1, simplify = FALSE)

results_df <- data.frame(
  Combination = character(), Iteration = integer(), Fold = integer(),
  AUC = numeric(), Accuracy = numeric(), Sensitivity = numeric(), Specificity = numeric(),
  stringsAsFactors = FALSE
)

for (combo_index in seq_along(marker_combinations)) {
  selected_markers <- marker_combinations[[combo_index]]
  data_subset <- data[c("Probe_ID", selected_markers)]

  cat("Processing combination", combo_index, "of", length(marker_combinations),
      ": Markers", paste(selected_markers, collapse = ", "), "\n")

  for (iteration in 1:1000) {
    cat("  Iteration", iteration, "of 1000\n")

    cv <- createMultiFolds(data_subset$Probe_ID, k = 5, times = 1)

    for (fold in 1:5) {
      train_indices <- cv[[paste0("Fold", fold, ".Rep1")]]
      test_indices <- setdiff(seq_len(nrow(data_subset)), train_indices)

      data_train <- data_subset[train_indices, ]
      data_test <- data_subset[test_indices, ]

      RF_model <- randomForest(Probe_ID ~ ., data = data_train, importance = TRUE, ntree = 1000)

      RF_predict <- predict(RF_model, data_test, type = "class")
      RF_predict_prob <- predict(RF_model, data_test, type = "prob")

      if (length(unique(RF_predict)) > 1) {
        if (any(RF_predict_prob[, "TRUE"] != 1 & RF_predict_prob[, "TRUE"] != 0)) {
          RF_ROC_predict <- prediction(predictions = RF_predict_prob[, "TRUE"], labels = data_test$Probe_ID)
          RF_AUC <- performance(RF_ROC_predict, measure = "auc")
          auc_value <- as.numeric(RF_AUC@y.values)
        } else {
          auc_value <- NA
        }

        confusion <- confusionMatrix(RF_predict, data_test$Probe_ID)
        acc_value <- confusion$overall["Accuracy"]
        sens_value <- confusion$byClass["Sensitivity"]
        spec_value <- confusion$byClass["Specificity"]
      } else {
        auc_value <- acc_value <- sens_value <- spec_value <- NA
      }

      results_df <- rbind(results_df, data.frame(
        Combination = paste(selected_markers, collapse = ", "),
        Iteration = iteration, Fold = fold,
        AUC = auc_value, Accuracy = acc_value,
        Sensitivity = sens_value, Specificity = spec_value
      ))
    }
  }

  combo_filename <- paste("metrics_for_", paste(markers, collapse = "_"), ".txt", sep = "")
  write.table(results_df, file = combo_filename, row.names = FALSE, sep = "\t", quote = FALSE)
}

sink()
