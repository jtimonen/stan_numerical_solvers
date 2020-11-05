library(tidyverse)
library(cmdstanr)
library(posterior)
library(loo)

lambda = 5.0
y = rpois(20, lambda)

model_ref = cmdstan_model("stan/poisson.stan")
ref_fit = model_ref$sample(data = list(N = length(y),
                                   y = y), parallel_chains = 4)

model = cmdstan_model("stan/approx_poisson.stan")

dxs = exp(seq(log(1e-2), log(1.0), length = 10))
dfs = list()
df = lapply(dxs, function(dx) {
  fit = model$sample(data = list(N = length(y),
                                 y = y,
                                 dx = dx), parallel_chains = 4)
  
  dfs[[i]] = fit$draws(c("log_lambda", "lambda", "log_ratio")) %>%
    as_draws_df %>%
    as_tibble() %>%
    mutate(which = "approximate",
           dx = dx)
}) %>%
  bind_rows() %>%
  select(log_lambda, lambda, log_ratio, which, dx)

df_ref = ref_fit$draws(c("log_lambda", "lambda")) %>%
  as_draws_df %>%
  as_tibble() %>%
  mutate(which = "approximate",
         dx = dx) %>%
  select(log_lambda, lambda, which, dx) %>%
  gather(parameter, value, log_lambda, lambda) %>%
  group_by(parameter) %>%
  summarize(lq = quantile(value, 0.2),
            m = quantile(value, 0.5),
            uq = quantile(value, 0.8))

df_actual = tibble(parameter = c("lambda", "log_lambda"),
                   value = c(lambda, log(lambda)))

df_resampled = df %>%
  group_by(dx) %>%
  summarize(pareto_k = psis(log_ratio, r_eff = NA)$diagnostics$pareto_k) %>%
  filter(pareto_k < 0.5) %>%
  left_join(df) %>%
  gather(parameter, value, log_lambda, lambda) %>%
  group_by(dx, parameter) %>%
  mutate(value = resample_draws(as_draws_matrix(as.matrix(value)),
                                weights = exp(log_ratio)) %>% as.vector) %>%
  mutate(resampled = TRUE)

df %>%
  gather(parameter, value, log_lambda, lambda) %>%
  mutate(resampled = FALSE) %>%
  bind_rows(df_resampled) %>%
  group_by(dx, parameter, resampled) %>%
  summarize(lq = quantile(value, 0.2),
            m = quantile(value, 0.5),
            uq = quantile(value, 0.8)) %>%
  ggplot() +
  geom_ribbon(aes(dx, ymin = lq, ymax = uq, group = resampled, fill = resampled), alpha = 0.5) +
  geom_errorbar(data = df_ref, aes(x = 0.01, ymin = lq, ymax = uq),
                size = 2, width = 0, color = "red") +
  geom_hline(data = df_actual, aes(yintercept = value), color = "black", linetype = "dashed") +
  scale_x_log10() +
  facet_grid(parameter ~ ., scales = "free_y")
