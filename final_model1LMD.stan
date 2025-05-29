data{
  int<lower=1> n_trees;
  int<lower=1> n_xcovs;
  
  int<lower=1> n_stands;
  vector[n_stands] stand_elev;
  
  int<lower=1> n_trans;
  int<lower=1, upper=n_stands> stand_id[n_trans];
  vector[n_trans] trans_slope;
  
  int<lower=1> n_trans_year;
  matrix[n_trans_year, n_xcovs] xmat;
  
  int<lower=1, upper=n_trans> trans_id[n_trans_year];
  
  int<lower=1, upper=n_trans_year> trans_year_id[n_trees];
  vector[n_trees] tree_dbh;
  
  int<lower=1> n_tot;
  int<lower=1> n_obs;
  int<lower=1, upper=n_obs> obs_id[n_tot];
  int<lower=1, upper=n_stands> stand_id2[n_tot];
  int<lower=1, upper=n_trans> trans_id2[n_tot];
  int<lower=1, upper=n_trees> tree_id[n_tot];
  int<lower=1> n_visits[n_trees];
  
  matrix[n_tot, n_xcovs] vmat;
  
  vector[n_tot] exp_ind;
  vector[n_tot] hike_dist;
  int<lower=0, upper=1> dets[n_tot];
  
  vector[n_xcovs] mean_dbh;
}

transformed data{
  int<lower=0, upper=1> naive_inf[n_trees];
  int n_xcovs1m = n_xcovs - 1;
  matrix[n_trans_year, n_xcovs1m] xmat2 = xmat[1:n_trans_year, 2:n_xcovs];
  matrix[n_tot, n_xcovs1m] vmat2 = vmat[1:n_tot, 2:n_xcovs];
  
  {
    int pos1 = 1;
    for(n in 1:n_trees){
      naive_inf[n] = max(segment(dets, pos1, n_visits[n]));
      pos1 += n_visits[n];
    }
  }
}

parameters{
  //Infection parameters
  real theta0;
  real theta1;
  vector[n_xcovs1m] theta_year;
  
  vector[n_stands] alpha_stand_raw;
  real<lower=0> sigma_stand;
  
  vector[n_trans] alpha_trans_raw;
  real<lower=0> sigma_trans;
  
  vector[n_trans_year] alpha_trans_year_raw;
  real<lower=0> sigma_trans_year;
  
  real alpha_dbh;
  
  //Detection parameters
  real eta0;
  real eta1;
  vector[n_xcovs1m] eta_year;
  
  vector[n_obs] beta_obs_raw;
  real<lower=0> sigma_obs;
  
  real beta_exp;
  real beta_hike;
  real beta_dbh;
}

transformed parameters{
  //Prevalence
  vector[n_stands] mu_stand = theta0 + theta1 * stand_elev;
  vector[n_stands] alpha_stand = mu_stand + alpha_stand_raw * sigma_stand;
  vector[n_trans] alpha_trans = alpha_stand[stand_id] + 
                                alpha_trans_raw * sigma_trans;
  vector[n_trans_year] alpha_trans_year = alpha_trans[trans_id] + xmat2 * theta_year +
                          alpha_trans_year_raw * sigma_trans_year;
  vector[n_trees] logit_psi = alpha_dbh * tree_dbh + alpha_trans_year[trans_year_id];
  
  vector[n_trees] psi = inv_logit(logit_psi);
  vector[n_trees] log_psi = log(psi);
  vector[n_trees] log1m_psi = log1m(psi);
  
  //Detection
  vector[n_tot] logit_p = eta0 +
                          eta1 * trans_slope[trans_id2] +
                          vmat2 * eta_year + 
                          beta_exp * exp_ind +
                          beta_hike * hike_dist +
                          beta_dbh * tree_dbh[tree_id] +
                          beta_obs_raw[obs_id] * sigma_obs;
}

model{
  //Create variable to index over
  int pos2 = 1;
  
  //Stand and transect level effects on infection
  alpha_stand_raw ~ std_normal();
  alpha_trans_raw ~ std_normal();
  alpha_trans_year_raw ~ std_normal();
  
  //Observer effects on detection
  beta_obs_raw ~ std_normal();
  
  //Model infection
  for(n in 1:n_trees){
    int dets_temp[n_visits[n]] = segment(dets, pos2, n_visits[n]);
    vector[n_visits[n]] logit_p_temp = segment(logit_p, pos2, n_visits[n]);
    
    if(naive_inf[n] == 1){
      target += log_psi[n] +
                   bernoulli_logit_lpmf(dets_temp | logit_p_temp);
    } else {
      target += log_sum_exp(log_psi[n] +
                               bernoulli_logit_lpmf(dets_temp | logit_p_temp),
                            log1m_psi[n]);
    }
    
    pos2 += n_visits[n];
  }
  
  //Infection priors
  theta0 ~ normal(0, 5);
  theta1 ~ normal(0, 5);
  theta_year ~ normal(0, 5);
  
  alpha_dbh ~ normal(0, 5);

  sigma_stand ~ gamma(2, 0.25);
  sigma_trans ~ gamma(2, 0.25);
  sigma_trans_year ~ gamma(2, 0.25);
  
  //Detection priors
  eta0 ~ normal(0, 5);
  eta1 ~ normal(0, 5);
  eta_year ~ normal(0, 5);
  
  beta_exp ~ normal(0, 5);
  beta_hike ~ normal(0, 5);
  beta_dbh ~ normal(0, 5);
  
  sigma_obs ~ gamma(2, 0.25);
}

generated quantities{
  real prev_t1 = inv_logit(theta0 + alpha_dbh * mean_dbh[1]);
  real prev_t2 = inv_logit(theta0 + theta_year[1] + alpha_dbh * mean_dbh[2]);
  real prev_t3 = inv_logit(theta0 + theta_year[2] + alpha_dbh * mean_dbh[3]);
  real prev_t4 = inv_logit(theta0 + theta_year[3] + alpha_dbh * mean_dbh[4]);
  real prev_t5 = inv_logit(theta0 + theta_year[4] + alpha_dbh * mean_dbh[5]);
  
  real small_t1 = inv_logit(theta0 + alpha_dbh * -1.26);
  real small_t2 = inv_logit(theta0 + theta_year[1] + alpha_dbh * -1.26);
  real small_t3 = inv_logit(theta0 + theta_year[2] + alpha_dbh * -1.26);
  real small_t4 = inv_logit(theta0 + theta_year[3] + alpha_dbh * -1.26);
  real small_t5 = inv_logit(theta0 + theta_year[4] + alpha_dbh * -1.26);
  
  real med_t1 = inv_logit(theta0 + alpha_dbh * -0.08);
  real med_t2 = inv_logit(theta0 + theta_year[1] + alpha_dbh * -0.08);
  real med_t3 = inv_logit(theta0 + theta_year[2] + alpha_dbh * -0.08);
  real med_t4 = inv_logit(theta0 + theta_year[3] + alpha_dbh * -0.08);
  real med_t5 = inv_logit(theta0 + theta_year[4] + alpha_dbh * -0.08);
  
  real large_t1 = inv_logit(theta0 + alpha_dbh * 1);
  real large_t2 = inv_logit(theta0 + theta_year[1] + alpha_dbh * 1);
  real large_t3 = inv_logit(theta0 + theta_year[2] + alpha_dbh * 1);
  real large_t4 = inv_logit(theta0 + theta_year[3] + alpha_dbh * 1);
  real large_t5 = inv_logit(theta0 + theta_year[4] + alpha_dbh * 1);
  
  real xlarge_t1 = inv_logit(theta0 + alpha_dbh * 2.28);
  real xlarge_t2 = inv_logit(theta0 + theta_year[1] + alpha_dbh * 2.28);
  real xlarge_t3 = inv_logit(theta0 + theta_year[2] + alpha_dbh * 2.28);
  real xlarge_t4 = inv_logit(theta0 + theta_year[3] + alpha_dbh * 2.28);
  real xlarge_t5 = inv_logit(theta0 + theta_year[4] + alpha_dbh * 2.28);
}
