install.packages("readxl")   
library(readxl)


# --- 1. Loading the Data -----

# Canadian Pensioners’ Mortality 2014 – Tables 1: CPM 2014 Rates 
Data_raw <- read_excel(path = 'CPM2014 Mortality Table.xlsx')
Mortality_Table <- data.frame(
  age = as.numeric(Data_raw$Age),
  qx_male = as.numeric(Data_raw$Male),
  qx_female = as.numeric(Data_raw$Female) 
)
rm(Data_raw)
Mortality_Table <- Mortality_Table[!is.na(Mortality_Table$age), ]
str(Mortality_Table)


# --- 2. qx lookup function ----

get_qx <- function(ages,sexes){
  index <- match(ages,Mortality_Table$age)
  if (any(is.na(index))) warning("Some ages are outside the table")
  qx <- ifelse(sexes=='M',
                  Mortality_Table$qx_male[index],
                  Mortality_Table$qx_female[index])
  
  return(qx)
}



# --- 3. Generate a portfolio ---------------------------------
generate_portfolio <- function(n, age_min = 55, age_max = 90,
                               contrib=100, contrib_sd=0,
                               prob_male = 0.5, age_dist = "uniform",
                               age_mean = 70, age_sd = 8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sex <- ifelse(runif(n) < prob_male, "M", "F")
  ages_range <- age_min:age_max
  
  age <- switch(age_dist,
            uniform = sample(ages_range, n, replace = TRUE),
            normal  = round(pmax(age_min, pmin(age_max, rnorm(n, age_mean, age_sd)))),
            peaked  = sample(ages_range, n, replace = TRUE,
                                 prob = dnorm(ages_range, age_mean, age_sd)),
                stop("unknown age_dist: ", age_dist)
  )
  
  a <- if (contrib_sd > 0) {
    round(pmax(1, rnorm(n, contrib, contrib_sd)))
  } else {
    rep(contrib, n)
  }
  data.frame(id = 1:n, age = age, sex = sex,
             qx = get_qx(age, sex), contributions=a,
             stringsAsFactors = FALSE)
}


# --- 4. Distribution of the available mortality credits Sn ----

    # Convolution of two discrete distributions (each a data.frame of value, prob). 

convolution2 <- function(d1,d2){
  values <-  outer(d1$value, d2$value, "+")
  probs  <- outer(d1$prob,  d2$prob,  "*")
  merged <- tapply(probs, round(values, 8), sum)
  data.frame(value = as.numeric(names(merged)),
             prob  = as.numeric(merged))
}

    # Distribution of S_n over the members in `rows`.
    # Each member contributes X: value 0 w.p. 1-qx, value a w.p. qx.

pool_distribution <- function(tontine_data, rows = seq_len(nrow(tontine_data))) {
  one_member <- lapply(rows, function(j)
    data.frame(value = c(0, tontine_data$contributions[j]), 
               prob = c(1 - tontine_data$qx[j], tontine_data$qx[j])))
  Reduce(convolution2, one_member)
}




# --- 5. The conditional mean risk-sharing rule ----

# 
#  CMRS mortality credits (Denuit-Hieber-Robert, formula 2.2)
#  E[X_i | S_n = s] = a_i * q_i * P[S_n^{(-i)} = s - a_i] / P[S_n = s]
#  where X_i = a_i if member i dies (prob q_i), else 0.
#
#  Requires a portfolio data frame `pf` with columns:
#     id : member identifier
#     contributions  : accumulated contribution a_i (fixed)
#     qx : one-period death probability q_i

# CMRS credit for one member i, at every possible total s
cmrs_one <- function(pf, i, full = pool_distribution(pf))  # law of S_n 
{                                
  loo  <- pool_distribution(pf, setdiff(seq_len(nrow(pf)), i))   # law of S_n without i
  p_loo <- loo$prob[match(round(full$value - pf$contributions[i], 8), 
                          round(loo$value, 8))]
  p_loo[is.na(p_loo)] <- 0
  data.frame(id     = pf$id[i],
             s      = full$value,
             prob   = full$prob,
             credit = pf$contributions[i] * pf$qx[i] * p_loo / full$prob)
}



# --- 6. The proportional risk-sharing rule ----

prop_one <- function(pf, i) {
  full <- pool_distribution(pf)                       # law of S_n (states + probs)
  w_i  <- (pf$qx[i] * pf$contributions[i]) /
    sum(pf$qx * pf$contributions)               # proportional weight
  data.frame(id     = pf$id[i],
             s      = full$value,
             prob   = full$prob,
             credit = w_i * full$value)               # w_i * s in each state
}



# --- 7. The convex combination of CMRS and proportional sharing-rule ----
combo_one <- function(pf, i, delta) {
  cm <- cmrs_one(pf, i)
  pr <- prop_one(pf, i)
  data.frame(id     = pf$id[i],
             s      = cm$s,
             prob   = cm$prob,
             credit = delta * cm$credit + (1 - delta) * pr$credit)
}



# --- 8. Expected utility optimization ----

# Power (CRRA) utility with parameter gamma.
utility <- function(w, gamma) {
  if (abs(gamma - 1) < 1e-9) log(w) else (w^(1 - gamma)) / (1 - gamma)
}

# --- EXACT expected utility for member i under delta*CMRS + (1-delta)*proportional 
expected_utility_exact <- function(pf, i, delta, gamma, beta = 1) {
  a_i <- pf$contributions[i]
  q_i <- pf$qx[i]
  
  cc  <- combo_one(pf, i, delta)                               # credit at every pool state s
  loo <- pool_distribution(pf, setdiff(seq_len(nrow(pf)), i))  # law of S_n^(-i)
  
  # credit member i receives in the survive states (total = v) and die states (total = v + a_i)
  idx_survive <- match(round(loo$value,       8), round(cc$s, 8))
  idx_die     <- match(round(loo$value + a_i, 8), round(cc$s, 8))
  
  credit_survive <- ifelse(is.na(idx_survive), 0, cc$credit[idx_survive])
  credit_die     <- ifelse(is.na(idx_die),     0, cc$credit[idx_die])
  
  U_survive <- sum(utility(a_i + credit_survive, gamma) * loo$prob)
  U_die     <- sum(utility(credit_die,           gamma) * loo$prob)
  
  return((1 - q_i) * U_survive + beta * q_i * U_die)
}




# Plot expected utility vs delta for one member (fast: heavy parts computed once).
plot_eu_vs_delta <- function(pf, i, gamma, beta = 1, step = 0.05) {
  full <- pool_distribution(pf)
  cm   <- cmrs_one(pf, i, full)
  loo  <- pool_distribution(pf, setdiff(seq_len(nrow(pf)), i))
  a_i  <- pf$contributions[i]; q_i <- pf$qx[i]
  w_i  <- (q_i * a_i) / sum(pf$qx * pf$contributions)
  prop_credit <- w_i * cm$s

  idx_s <- match(round(loo$value,       8), round(cm$s, 8))
  idx_d <- match(round(loo$value + a_i, 8), round(cm$s, 8))

  deltas <- seq(0, 1, by = step)
  eu <- sapply(deltas, function(d) {
    credit <- d * cm$credit + (1 - d) * prop_credit
    cs <- ifelse(is.na(idx_s), 0, credit[idx_s])
    cd <- ifelse(is.na(idx_d), 0, credit[idx_d])
    (1 - q_i) * sum(utility(a_i + cs, gamma) * loo$prob) +
      beta * q_i * sum(utility(cd, gamma) * loo$prob)
  })

  plot(deltas, eu, type = "l", lwd = 2,
       xlab = expression(delta), ylab = "Expected utility",
       main = paste0("Member ", i, ": EU vs delta (gamma=", gamma, ")"))
  invisible(data.frame(delta = deltas, eu = eu))
}





# --- 9. Simulating a portfolio of 100 participants ----

pf <- generate_portfolio(n=100,age_min = 55, age_max = 85,
                                    contrib=50, contrib_sd=3,
                                    prob_male = 0.5,age_dist = 'peaked',
                                    age_mean = 70, age_sd = 8, seed = 125)
head(pf)
table(pf$sex)
table(pf$age)
summary(pf$age)
mean(pf$qx)
#hist(pf$age)




# --- 10. Optimal convex combination ----

plot_eu_vs_delta(pf, i = 1, gamma = 1.5)
plot_eu_vs_delta(pf, i = 1, gamma = 6.5)
plot_eu_vs_delta(pf, i = 1, gamma = 6.5, beta = 0.8)