# Required libraries
library(ggplot2)
library(cluster)      # for silhouette()
library(gridExtra)    # for plotting side-by-side

# Simulate data
set.seed(123)
n <- 100  # samples
p <- 50   # features

Y <- matrix(rnorm(n * p), nrow = n)
group <- rep(c(0, 1), each = n/2)
covariate <- rnorm(n)

# Design matrices
X1 <- model.matrix(~ factor(group))
X2 <- model.matrix(~ factor(group) + covariate)

# Residual extraction function
get_residuals <- function(Y, design) {
  beta_hat <- solve(t(design) %*% design) %*% t(design) %*% Y
  residuals <- Y - design %*% beta_hat
  return(residuals)
}

# Get residuals
resid_no_covar <- get_residuals(Y, X1)
resid_with_covar <- get_residuals(Y, X2)

# PCA
pca_no <- prcomp(resid_no_covar, scale. = TRUE)
pca_with <- prcomp(resid_with_covar, scale. = TRUE)

pc_scores_no <- pca_no$x[, 1:2]
pc_scores_with <- pca_with$x[, 1:2]

# K-means clustering
k <- 2
km_no <- kmeans(pc_scores_no, centers = k)
km_with <- kmeans(pc_scores_with, centers = k)

# Silhouette analysis
sil_no <- silhouette(km_no$cluster, dist(pc_scores_no))
sil_with <- silhouette(km_with$cluster, dist(pc_scores_with))

avg_sil_no <- mean(sil_no[, "sil_width"])
avg_sil_with <- mean(sil_with[, "sil_width"])

cat("Silhouette Score WITHOUT covariate: ", round(avg_sil_no, 3), "\n")
cat("Silhouette Score WITH covariate:    ", round(avg_sil_with, 3), "\n")

# Visualization using ggplot2
df_no <- data.frame(PC1 = pc_scores_no[,1], PC2 = pc_scores_no[,2], Cluster = factor(km_no$cluster))
df_with <- data.frame(PC1 = pc_scores_with[,1], PC2 = pc_scores_with[,2], Cluster = factor(km_with$cluster))

p1 <- ggplot(df_no, aes(PC1, PC2, color = Cluster)) + 
  geom_point(size = 2) + 
  ggtitle("Without Covariate") + 
  theme_minimal()

p2 <- ggplot(df_with, aes(PC1, PC2, color = Cluster)) + 
  geom_point(size = 2) + 
  ggtitle("With Covariate") + 
  theme_minimal()