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

# ====/// 3a: Neyman estratificado \\\=====

df_strat <- df %>% filter(!is.na(total_expenditure), !is.na(age)) %>%
  group_by(estrato,treat) %>%
  summarise(Mean = mean(total_expenditure),
            Var = var(total_expenditure),
            num = n(),
            .groups = "drop")

stargazer(df_strat,summary = F, 
          title = "Estadísticas por Estrato",
          out="tabla_desc_neyman.tex")



tabla_estratos <- df %>% filter(!is.na(total_expenditure), !is.na(age)) %>%
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
    var_tau_h = varianza_T / n_T + varianza_C / n_C) %>%
  select(estrato,tau_h,var_tau_h,n_h,peso)

stargazer(tabla_estratos,summary = F, 
          title = "Efectos de Tratamiento por Estrato",
          out="tabla_ates_neyman.tex")

(tau_strat <- sum(tabla_estratos$peso * tabla_estratos$tau_h))
var_strat <- sum((tabla_estratos$peso^2) * tabla_estratos$var_tau_h)
(se_strat <- sqrt(var_strat))
(t_strat <- tau_strat / se_strat)
(p_strat <- 2*(1-pnorm(abs(t_strat))))

# ====/// 3b: OLS estratificado \\\=====

mod_inter <- lm(total_expenditure ~ 0 + estrato + estrato:treat, data = df) 
summary(mod_inter)


# ====/// 4a: OLS estratificado \\\=====

# OLS + FE (controles)
controls_simple <- intersect(c("age", "birthplace", "gender", 
                               "caste"), names(df)) 

run_feols <- \(y) {
  feols(
    as.formula(paste(y, "~ treat +", paste(controls_simple, collapse = " + "), "| Pole")),
    data = df, cluster = ~Pole
  )
}

ols_adult <- feols(adult_activity ~ treat + age + birthplace + gender + caste | Pole,
                   data = df, cluster = ~Pole)
  
ols_child <- feols(child_activity ~ treat + age + birthplace + gender + caste | Pole,
                   data = df, cluster = ~Pole)

etable(ols_adult, ols_child,
       headers = c("Actividad adultos", "Actividad niños"),
       tex = TRUE,
       file = "Tabla_p4a.tex")

# ====/// 4b: Double lasso \\\=====

double_lasso_fe <- function(yname, data = df, xvars = controls, dname = "treat", fe = "Pole") {
  dsub <- data %>% select(all_of(c(yname, dname, fe, xvars))) %>% na.omit()
  fe_id <- dsub[[fe]]
  demean <- \(v) v - ave(v, fe_id, FUN = mean) 
  
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


# ====/// 5a: Atrición  \\\=====

attrition_tbl <- df %>%
  group_by(treat) %>%
  summarise(
    N = n(),
    tasa_atricion = 1-mean(total_expenditure_isobs))

stargazer(attrition_tbl,summary = F, 
          title = "Tasa de Atricion",
          out="atricion.tex")

# ====/// 5b: Lee Bounds  \\\=====

mean_control <- mean()

df_obs <- df %>% filter(total_expenditure_isobs == 1)

mean_control <- mean(df_obs$total_expenditure[df_obs$treat==0])
p_ar_treat <- (1-attrition_tbl$tasa_atricion[attrition_tbl$treat==0])/(1-attrition_tbl$tasa_atricion[attrition_tbl$treat==1])

df_lee <- df_obs %>% filter(treat==1) %>%
      arrange(total_expenditure) %>%
  summarise(lower_lim = mean(head(total_expenditure,floor(p_ar_treat*n()))),
            upper_lim = mean(tail(total_expenditure,floor(p_ar_treat*n()))))

lee_tbl <- tibble(
  Estimador = c("Lee - Límite inferior", "Lee - Límite superior"),
  Valor = round(c(df_lee$lower_lim - mean_control,
                  df_lee$upper_lim - mean_control), 2))

# ====/// 6: Poder estadistico  \\\=====

power_calc <- function(n, delta, alpha, p) {
  zcrit <- qnorm(1 - alpha / 2)
  pnorm(delta*sqrt(n*p*(1-p)) - zcrit)
}

delta_ney <- (tau_hat_neyman / sqrt(df_agg$Var[df_agg$treat==0]))
alpha <- 0.07
target_power <- 0.83
p_treat <- 0.5

psi <- tibble(
  n = 20:200,
  power = power_calc(20:200, delta = delta_ney, alpha = alpha, p = p_treat)
)

ggplot(psi, aes(x = n, y = power)) +
  geom_line() +
  geom_hline(yintercept = target_power, linetype = "dashed") +
  labs(
    x = "Tamaño muestral total (n)",
    y = "Poder estadístico ψ(n)",
    title = paste0("Curva de poder: α = ", alpha,
                   ", delta = ", round(delta_ney, 3),
                   ", p = ", p_treat)
  ) +
  theme_minimal()

ggsave("Graf_p6.png",  width = 5.54, height = 4.95)

(n_star <- psi %>% filter(power >= target_power) %>% slice(1) %>% pull(n))
