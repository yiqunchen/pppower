# Labeled sample sizes from pppower for the ppi_py cross-check (check_ppi_py.py)
g <- expand.grid(rho = c(0.3, 0.5, 0.7, 0.9), N = c(200, 1000, 10000), delta = c(0.2, 0.3))
g$n_cl <- ceiling(1 / ((g$delta / (qnorm(0.975) + qnorm(0.8)))^2))
g$n_star <- mapply(function(r, N, d) ceiling(suppressWarnings(pppower::power_ppi_mean(
  delta = d, N = N, n = NULL, power = 0.8, sigma_y2 = 1, sigma_f2 = 1, cov_y_f = r))),
  g$rho, g$N, g$delta)
write.csv(g, file.path("simulation_studies", "external_checks", "pppower_nstar.csv"), row.names = FALSE)
