parameters {
#include params.stan
}
transformed parameters {
#include tparams.stan
}
model {
#include prior.stan
}
