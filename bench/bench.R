library(ggplot2)
bench <- read.table("2016-03-01-clang-700.1.81-x86_64-apple-darwin15.3.0.txt",
                    col.names = c("mode", "rows", "cols", "bandwidth",
                                  "X1", "X2", "X3", "X4", "X5"))
bench$rows <- factor(bench$rows)
bench$runtime <- rowMeans(subset(bench, select = c(X1, X2, X3, X4, X5)))
bench$runtime.stdev <- apply(subset(bench, select = c(X1, X2, X3, X4, X5)), 1, sd)

p <- ggplot(bench, aes(bandwidth, runtime, colour = mode)) +
  facet_grid(rows ~ ., scales = "free_y") +
  geom_line() +
  geom_point() +
  scale_x_log10(breaks = scales::pretty_breaks(n = 8)) +
  xlab("Bandwidth") +
  ylab("Runtime (s)") +
  ggtitle("Benchmarks partitioned by number of points")
ggsave("../img/2016-03-01-clang-700.1.81.x86_64-apple-darwin15.3.0.png", dpi = 100)
