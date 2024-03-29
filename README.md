# Experiments with numerical solvers in Stan (mostly ODE)

* `case_study_num` - a draft study with ODE and PDE solvers (see [draft](https://users.aalto.fi/~timonej3/case_study_num.html))
* `case_study_ode` - a case study with ODE solvers (using **cmdstanr** and the "new" ODE interface of Stan)
* `case_study_ode_rstan` - a case study with ODE solvers (using **rstan**)
* `benchmark_vec` - speed comparison between a vectorized and non-vectorized
way of computing likelihood of an ODE model (using **cmdstanr** and the "new" ODE interface of Stan)
* `max_num_steps_reached` - studying ODE solver and Stan sampler
behaviour when `max_num_steps` is reached
