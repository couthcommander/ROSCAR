#' @keywords internal

train_rf <- function(y, x) {
    if(!requireNamespace('caret', quietly = TRUE)) {
        stop('package "caret" is required to utilize random forest approach')
    }
    caret::train(y = y,
        x = x,
        method = "rf",
        trControl = caret::trainControl(method = "cv", number = 5),
        tuneGrid = data.frame(mtry = floor(sqrt(ncol(x)))),
        ntree = 50,
        nodesize = 25
    )
}
