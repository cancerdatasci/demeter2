n_samples <- 10
seed_df <- data.frame(random_seed = round(runif(n_samples, min = 0, max = 100000)))
write.csv(seed_df, '~/CPDS/demeter2/seed_params.csv', row.names=F)
