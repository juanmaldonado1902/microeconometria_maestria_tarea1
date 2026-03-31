# Tarea 1, Econometría II
# Creada: Juan Pablo Maldonado
# Revisada y editada: Arturo Aguilar

# Clean environment and workspace
rm(list = ls())

# Librerias
pkgs <- c("readr", "dplyr", "tibble", "ggplot2", "sandwich", "lmtest",
          "stargazer","tidyr","car","tidyverse","broom","lmtest",
          "fixest","glmnet")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))


# Cargar y limpiar base
df <- read_csv("india_base_final.csv", show_col_types = FALSE) %>%
  mutate(
    across(c(treat, age, total_expenditure, total_expenditure_obs), as.numeric),
    total_expenditure_isobs = ifelse(is.na(total_expenditure_obs), 0, 1),
    Pole = as.factor(Pole),
    estrato = case_when(
      age <= 35 ~ "Joven",
      age <= 50 ~ "Mediana",
      TRUE ~ "Mayor"
    ) %>% factor(levels = c("Joven", "Mediana", "Mayor"))
  ) %>%
  filter(treat %in% c(0, 1))

# Definición de variables

covars <- intersect(c("gender", "age", "religion", "caste", "education", "homeBuilt"), names(df)) # generadas para la tabla de balance
controls_simple <- intersect(c("age", "birthplace", "gender", "caste"), names(df)) # generadas para el modelo de controles
controls <- intersect(c("gender", "birthplace", "age", "religion", "caste", "education", "homeBuilt"), names(df)) # variables exclusivas del pretratamiento que dejo para double lasso

# ====/// 1a: Tabla de balance \\\=====

balance_table <- map_dfr(covars, \(v) {
  d <- df %>% select(treat, all_of(v)) %>% drop_na()
  tibble(
    Variable = v,
    Media_Total = mean(d[[v]]),
    Media_T1 = mean(d[[v]][d$treat == 1]),
    Media_T0 = mean(d[[v]][d$treat == 0]),
    Diferencia = mean(d[[v]][d$treat == 1]) - mean(d[[v]][d$treat == 0]),
    p_value = t.test(d[[v]] ~ d$treat)$p.value
  )
}) %>%
  mutate(across(where(is.numeric), round, 4))

balance_table

# ====/// 1c: F-test balance \\\=====

model_balance <- lm(reformulate(covars, response = "treat"), data = df) # modelo para ver si explican treat
robust_se <- vcovHC(model_balance, type = "HC1")
coeftest(model_balance, vcov = robust_se)

coef_names <- setdiff(names(coef(model_balance)), "(Intercept)")
joint_test <- linearHypothesis(model_balance, paste0(coef_names, " = 0"), vcov = robust_se)
joint_test

# ====/// 2b: Neyman y MCO \\\=====

df_agg <- df %>% filter(!is.na(total_expenditure)) %>%
  group_by(treat) %>%
  summarise(Mean = mean(total_expenditure),
            Var = var(total_expenditure),
            num = n())

stargazer(df_agg,summary = F, 
          title = "Efectos de Tratamiento",
          out="tabla_ates.tex")

(tau_hat_neyman <- df_agg$Mean[df_agg$treat==1] -
                    df_agg$Mean[df_agg$treat==0])
(se_neyman <- sqrt((df_agg$Var[df_agg$treat==1] / df_agg$num[df_agg$treat==1]) + 
                     (df_agg$Var[df_agg$treat==0] / df_agg$num[df_agg$treat==0])))
(t <- tau_hat_neyman / se_neyman)
(p <- 2*(1-pnorm(abs(t))))

# Regresión OLS  Neyman (esta parte no es completamente necesaria porque ya puse el código para extraer los coeficientes)

reg_neyman <- lm(total_expenditure ~ treat, 
                 data = df %>% filter(!is.na(total_expenditure)))
reg_neyman_rob <- coeftest(reg_neyman, vcov = vcovHC(reg_neyman, type = "HC1"))
reg_neyman_rob

#  Neyman estratificado

df_strat <- df %>% filter(!is.na(total_expenditure), !is.na(age))

tabla_estratos <- df_strat %>%
  group_by(estrato, treat) %>%
  summarise(
    n = n(),
    media = mean(total_expenditure),
    varianza = var(total_expenditure),
    .groups = "drop"
  ) %>%
  mutate(grupo = ifelse(treat == 1, "T", "C")) %>%
  select(-treat) %>%
  pivot_wider(names_from = grupo, values_from = c(n, media, varianza), names_sep = "_") %>%
  mutate(
    n_h = n_T + n_C,
    peso = n_h / sum(n_h),
    tau_h = media_T - media_C,
    var_tau_h = varianza_T / n_T + varianza_C / n_C
  )

tabla_estratos

tau_strat <- sum(tabla_estratos$peso * tabla_estratos$tau_h)
var_strat <- sum((tabla_estratos$peso^2) * tabla_estratos$var_tau_h)
se_strat <- sqrt(var_strat)
t_strat <- tau_strat / se_strat

tau_strat; var_strat; se_strat; t_strat

mod_inter <- lm(total_expenditure ~ 0 + estrato + estrato:treat, data = df_strat) # modelo sin el intercepto para reconstruir el valor de neyman estratificado
summary(mod_inter)

coef_h <- coef(mod_inter)[grep(":treat", names(coef(mod_inter)))]
tau_strat_reg <- sum(tabla_estratos$peso * coef_h)
tau_strat_reg

# OLS + FE (controles)
run_feols <- \(y) {
  feols(
    as.formula(paste(y, "~ treat +", paste(controls_simple, collapse = " + "), "| Pole")),
    data = df, cluster = ~Pole
  )
}

ols_adult <- run_feols("adult_activity")
ols_child <- run_feols("child_activity")

etable(ols_adult, ols_child, keep = "treat", se = "cluster")

# Double LASSO
double_lasso_fe <- function(yname, data = df, xvars = controls, dname = "treat", fe = "Pole") {
  dsub <- data %>% select(all_of(c(yname, dname, fe, xvars))) %>% na.omit()
  fe_id <- dsub[[fe]]
  demean <- \(v) v - ave(v, fe_id, FUN = mean) # transformación within
  
  y <- demean(dsub[[yname]])
  d <- demean(as.numeric(dsub[[dname]]))
  X <- model.matrix(reformulate(xvars, response = NULL), data = dsub)[, -1, drop = FALSE]
  X <- apply(X, 2, demean)
  
  sel_y <- which(as.vector(coef(cv.glmnet(X, y, alpha = 1), s = "lambda.min"))[-1] != 0) #proceso de selección 
  sel_d <- which(as.vector(coef(cv.glmnet(X, d, alpha = 1), s = "lambda.min"))[-1] != 0)
  sel <- sort(unique(c(sel_y, sel_d)))
  X_sel <- if (length(sel) == 0) NULL else X[, sel, drop = FALSE]
  
  y_res <- if (is.null(X_sel)) y else resid(lm(y ~ X_sel)) # partialling-out
  d_res <- if (is.null(X_sel)) d else resid(lm(d ~ X_sel))
  
  final <- lm(y_res ~ d_res)
  V <- vcovCL(final, cluster = fe_id, type = "HC1")
  
  list(
    coeftest = coeftest(final, V),
    selected_cols = if (is.null(X_sel)) character(0) else colnames(X_sel)
  )
}

dl_adult <- double_lasso_fe("adult_activity")
dl_child <- double_lasso_fe("child_activity")

dl_adult$coeftest
dl_child$coeftest
dl_adult$selected_cols
dl_child$selected_cols

#  Atrición / Lee Bounds
attrition_tbl <- df %>%
  group_by(treat) %>%
  summarise(
    N = n(),
    tasa_observado = mean(total_expenditure_isobs),
    tasa_atricion = 1 - tasa_observado,
    .groups = "drop"
  )

attrition_tbl

df_obs <- df %>% filter(total_expenditure_isobs == 1)
tau_uncorrected <- with(df_obs, mean(total_expenditure_obs[treat == 1]) - mean(total_expenditure_obs[treat == 0]))
tau_uncorrected

p1 <- mean(df$total_expenditure_isobs[df$treat == 1])
p0 <- mean(df$total_expenditure_isobs[df$treat == 0])

trim_group <- function(x, share, drop = c("top", "bottom")) {
  drop <- match.arg(drop)
  x <- sort(x)
  n_trim <- floor(length(x) * share)
  if (n_trim <= 0) return(x)
  if (drop == "top") x[1:(length(x) - n_trim)] else x[(n_trim + 1):length(x)]
}

treated_y <- df %>% filter(treat == 1, total_expenditure_isobs == 1) %>% pull(total_expenditure_obs)
control_y <- df %>% filter(treat == 0, total_expenditure_isobs == 1) %>% pull(total_expenditure_obs)

if (abs(p1 - p0) < 1e-12) {
  lee_lower <- tau_uncorrected
  lee_upper <- tau_uncorrected
} else if (p1 > p0) {
  trim_share <- (p1 - p0) / p1
  lee_lower <- mean(trim_group(treated_y, trim_share, "top")) - mean(control_y)
  lee_upper <- mean(trim_group(treated_y, trim_share, "bottom")) - mean(control_y)
} else {
  trim_share <- (p0 - p1) / p0
  lee_lower <- mean(treated_y) - mean(trim_group(control_y, trim_share, "bottom"))
  lee_upper <- mean(treated_y) - mean(trim_group(control_y, trim_share, "top"))
}

lee_tbl <- tibble(
  Estimador = c("Sin corrección (obs)", "Lee - Límite inferior", "Lee - Límite superior"),
  Valor = round(c(tau_uncorrected, lee_lower, lee_upper), 4)
)

lee_tbl

# Poder estadístico
tau_hat <- tau_hat_neyman
sigma_pooled <- sqrt(
  ((length(y1) - 1) * var(y1) + (length(y0) - 1) * var(y0)) /
    (length(y1) + length(y0) - 2)
)
delta <- tau_hat / sigma_pooled

tau_hat; sigma_pooled; delta

power_diffmeans <- function(n, delta, alpha = 0.07, p = 0.5) {
  zcrit <- qnorm(1 - alpha / 2)
  delta_n <- delta / sqrt(1 / (p * n) + 1 / ((1 - p) * n))
  1 - (pnorm(zcrit - delta_n) - pnorm(-zcrit - delta_n))
}

alpha <- 0.07
p_treat <- 0.5
target_power <- 0.83

psi <- tibble(
  n = 20:5000,
  power = power_diffmeans(20:5000, delta = delta, alpha = alpha, p = p_treat)
)

psi

n_star <- psi %>% filter(power >= target_power) %>% slice(1) %>% pull(n)
n_star

plot(
  psi$n, psi$power, type = "l",
  xlab = "Tamaño muestral total (n)",
  ylab = "Poder estadístico ψ(n)",
  main = paste0("Curva de poder: α = ", alpha, ", delta = ", round(delta, 3), ", p = ", p_treat)
)
abline(h = target_power, lty = 2)
abline(v = n_star, lty = 2)
text(n_star, target_power, labels = paste0(" n* ≈ ", n_star, "\n power = ", target_power), pos = 4)

psi %>% filter(n %in% (n_star - 5):(n_star + 5))