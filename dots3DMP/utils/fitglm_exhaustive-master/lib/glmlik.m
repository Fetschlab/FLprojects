function loglik = glmlik(X_dat, y_dat, y_pred0, distr)
% loglik = glmlik(X_dat, y_dat, y_pred, distr)

% 2016 (c) Yul Kang. hk2699 at columbia dot edu.

assert(strcmp(distr, 'binomial'), 'distr=%s is not implemented yet!', distr);
loglik = glmlik_binomial(X_dat, y_dat, y_pred0);
