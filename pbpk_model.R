# =============================================================================
# PBPK Two-Compartment Model with DDI Simulation
# =============================================================================
# Author : Krima R. Patel, Pharm.D. | M.S. Pharmaceutical Sciences Candidate
#          University at Buffalo, SUNY
#
# Required packages:
#   install.packages(c("deSolve", "ggplot2", "dplyr", "tidyr", "patchwork"))
# =============================================================================

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

COL_BASE  <- "#1B3A6B"
COL_DDI   <- "#BA7517"
COL_DOSES <- c("#1B3A6B", "#0F6E56", "#9B2335")
THEME_CLN <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92"),
        axis.line = element_line(color = "grey60"),
        plot.title = element_text(size = 13, face = "plain"),
        plot.subtitle = element_text(size = 10, color = "grey45"),
        legend.position = "bottom",
        legend.text = element_text(size = 9))

# --------------------------------------------------------------------------
# Drug library
# --------------------------------------------------------------------------
drug_params <- list(
  midazolam = list(name="Midazolam", enzyme="CYP3A4", F=0.36, ka=3.5,
                   Vc=20, Vp=70, CL=28, Q=15, fm=0.95, mw=325.8,
                   doses=c(2.5,5,10), dose_unit="mg"),
  caffeine  = list(name="Caffeine",  enzyme="CYP1A2", F=1.00, ka=4.2,
                   Vc=35, Vp=15, CL=1.9, Q=5,  fm=0.95, mw=194.2,
                   doses=c(50,100,200), dose_unit="mg"),
  warfarin  = list(name="Warfarin",  enzyme="CYP2C9", F=0.99, ka=1.5,
                   Vc=8,  Vp=5,  CL=0.19,Q=1.5,fm=0.90, mw=308.3,
                   doses=c(2.5,5,10), dose_unit="mg")
)

# --------------------------------------------------------------------------
# DDI presets
# --------------------------------------------------------------------------
ddi_presets <- list(
  ketoconazole = list(label="Strong CYP3A4 inhibitor (ketoconazole)", enzyme="CYP3A4", R_inh=0.95),
  fluconazole  = list(label="Moderate CYP3A4 inhibitor (fluconazole)",enzyme="CYP3A4", R_inh=0.55),
  rifampicin   = list(label="CYP3A4 inducer (rifampicin)",            enzyme="CYP3A4", R_inh=-0.60),
  fluvoxamine  = list(label="Strong CYP1A2 inhibitor (fluvoxamine)",  enzyme="CYP1A2", R_inh=0.90),
  amiodarone   = list(label="CYP2C9 inhibitor (amiodarone)",          enzyme="CYP2C9", R_inh=0.75)
)

# --------------------------------------------------------------------------
# Effective clearance (DDI adjustment)
# CL_eff = CL * (1 - fm * R_inhibition)
# --------------------------------------------------------------------------
effective_cl <- function(drug, perpetrator = NULL) {
  if (is.null(perpetrator)) return(drug$CL)
  mod   <- ddi_presets[[perpetrator]]
  R_inh <- if (!is.null(mod) && mod$enzyme == drug$enzyme) mod$R_inh else 0
  max(drug$CL * (1 - drug$fm * R_inh), 0.001)
}

# --------------------------------------------------------------------------
# ODE system — two-compartment PBPK
# State: c(A_gut, A_central, A_periph)
#
# Oral:
#   dA_gut/dt = -ka * A_gut
#   dA_c/dt   = F*ka*A_gut - (CL_eff/Vc)*A_c - (Q/Vc)*A_c + (Q/Vp)*A_p
#   dA_p/dt   = (Q/Vc)*A_c - (Q/Vp)*A_p
#
# IV bolus: A_gut=0, A_c(0)=Dose
# IV infusion (1h): zero-order rate into central for t <= 1
# --------------------------------------------------------------------------
pbpk_ode <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    inf_rate <- if (route == "iv_infusion" && t <= 1) dose / 1.0 else 0
    dA_gut   <- if (route == "oral") -ka * A_gut else 0
    dA_c     <- (if (route == "oral") F_oral * ka * A_gut else 0) +
                inf_rate - (CL_eff/Vc)*A_central - (Q/Vc)*A_central + (Q/Vp)*A_periph
    dA_p     <- (Q/Vc)*A_central - (Q/Vp)*A_periph
    list(c(dA_gut=dA_gut, A_central=dA_c, A_periph=dA_p))
  })
}

# --------------------------------------------------------------------------
# Simulation
# --------------------------------------------------------------------------
simulate_pbpk <- function(drug, dose, t_end=24, route="oral",
                           perpetrator=NULL, n_points=500) {
  CL_eff <- effective_cl(drug, perpetrator)
  parms  <- list(ka=drug$ka, F_oral=drug$F, Vc=drug$Vc, Vp=drug$Vp,
                 CL_eff=CL_eff, Q=drug$Q, route=route, dose=dose)
  y0 <- switch(route,
    oral        = c(A_gut=dose, A_central=0, A_periph=0),
    iv_bolus    = c(A_gut=0,    A_central=dose, A_periph=0),
    iv_infusion = c(A_gut=0,    A_central=0,    A_periph=0)
  )
  times <- seq(0, t_end, length.out=n_points)
  sol   <- ode(y=y0, times=times, func=pbpk_ode, parms=parms, method="lsoda")
  C     <- pmax(sol[,"A_central"] / drug$Vc, 0) * 1000  # mg/L -> ng/mL

  AUC  <- sum(diff(times) * (C[-length(C)] + C[-1]) / 2)
  Cmax <- max(C);  Tmax <- times[which.max(C)]

  idx <- seq(floor(0.7*n_points), n_points)
  fit <- tryCatch({
    m <- lm(log(C[idx][C[idx]>0]) ~ times[idx][C[idx]>0])
    kel <- -coef(m)[2];  if (kel>0) log(2)/kel else NA_real_
  }, error=function(e) NA_real_)

  data.frame(time=times, C_ng_mL=C, CL_eff=CL_eff,
             AUC=AUC, Cmax=Cmax, Tmax=Tmax, t_half=fit)
}

pk_summary <- function(df, label) {
  tibble(Scenario=label,
         `AUC (ng·h/mL)`=round(df$AUC[1],1),
         `Cmax (ng/mL)` =round(df$Cmax[1],1),
         `Tmax (h)`     =round(df$Tmax[1],2),
         `t½ (h)`       =round(df$t_half[1],2))
}

# --------------------------------------------------------------------------
# Plot: DDI comparison
# --------------------------------------------------------------------------
plot_ddi <- function(drug_key, dose, perpetrator, route="oral",
                     t_end=24, save_png=FALSE) {
  drug  <- drug_params[[drug_key]]
  base  <- simulate_pbpk(drug, dose, t_end, route)
  ddi_r <- simulate_pbpk(drug, dose, t_end, route, perpetrator)
  auc_r <- ddi_r$AUC[1]/base$AUC[1]
  cmax_r<- ddi_r$Cmax[1]/base$Cmax[1]
  lbl   <- ddi_presets[[perpetrator]]$label

  df <- bind_rows(
    base  %>% mutate(scenario=paste0("Baseline  |  AUC=",round(base$AUC[1],1)," ng·h/mL")),
    ddi_r %>% mutate(scenario=paste0("+ DDI     |  AUC=",round(ddi_r$AUC[1],1)," ng·h/mL"))
  )
  scn_col <- setNames(c(COL_BASE, COL_DDI), unique(df$scenario))

  p1 <- ggplot(df, aes(time, C_ng_mL, color=scenario, linetype=scenario)) +
    geom_line(linewidth=1.1) +
    scale_color_manual(values=scn_col) +
    scale_linetype_manual(values=c("solid","dashed")) +
    labs(x="Time (h)", y="Plasma concentration (ng/mL)",
         title=paste0(drug$name,"  |  ",dose," mg  |  Route: ",route),
         subtitle=paste0("DDI scenario: ",lbl), color=NULL, linetype=NULL) +
    THEME_CLN

  mdf <- tibble(Metric=c("AUC ratio","Cmax ratio"), Value=c(auc_r,cmax_r))
  p2  <- ggplot(mdf, aes(Metric, Value, fill=Metric)) +
    geom_col(width=0.45, show.legend=FALSE) +
    geom_hline(yintercept=1, linetype="dashed", color="grey50") +
    geom_hline(yintercept=2, linetype="dotted", color="#cc3333", linewidth=0.8) +
    geom_hline(yintercept=0.5, linetype="dotted", color="#cc3333", linewidth=0.8) +
    geom_text(aes(label=paste0(round(Value,2),"x")), vjust=-0.4, fontface="bold", size=4.5) +
    scale_fill_manual(values=c("AUC ratio"=COL_DDI,"Cmax ratio"=COL_BASE)) +
    scale_y_continuous(limits=c(0, max(c(auc_r,cmax_r))*1.35+0.5)) +
    annotate("text",x=2.4,y=2.08,label="2x threshold",size=3,color="#cc3333",hjust=1) +
    labs(x=NULL, y="DDI ratio (with / without)",
         title="Regulatory DDI thresholds",
         subtitle="FDA: AUC ratio >2x = significant DDI") +
    THEME_CLN

  out <- p1 + p2 + plot_layout(widths=c(2,1))
  if (save_png) { fname <- paste0("PBPK_",drug_key,"_",perpetrator,"_DDI.png"); ggsave(fname,out,width=11,height=4.5,dpi=150); message("Saved: ",fname) }
  print(out)
  cat("\n── PK Summary ───────────────────────\n")
  print(bind_rows(pk_summary(base,"Baseline"), pk_summary(ddi_r,paste0("+ ",lbl))))
  cat(sprintf("\nAUC ratio = %.2fx  |  Cmax ratio = %.2fx\n\n", auc_r, cmax_r))
  invisible(list(base=base, ddi=ddi_r))
}

# --------------------------------------------------------------------------
# Plot: Dose escalation
# --------------------------------------------------------------------------
plot_dose_escalation <- function(drug_key, route="oral", t_end=24, save_png=FALSE) {
  drug <- drug_params[[drug_key]]
  df   <- lapply(seq_along(drug$doses), function(i) {
    d   <- drug$doses[i]
    res <- simulate_pbpk(drug, dose=d, t_end=t_end, route=route)
    res %>% mutate(dose_label=paste0(d," ",drug$dose_unit,
      "  |  AUC=",round(res$AUC[1],0),"  |  Cmax=",round(res$Cmax[1],1)," ng/mL"), dose_idx=i)
  }) %>% bind_rows()
  dcol <- setNames(COL_DOSES[seq_along(drug$doses)], unique(df$dose_label))
  p <- ggplot(df, aes(time, C_ng_mL, color=dose_label)) +
    geom_line(linewidth=1.1) + scale_color_manual(values=dcol) +
    labs(x="Time (h)", y="Plasma concentration (ng/mL)",
         title=paste0(drug$name," — Dose escalation  |  Route: ",route),
         subtitle="Two-compartment PBPK model", color=NULL) + THEME_CLN
  if (save_png) { fname <- paste0("PBPK_",drug_key,"_dose_escalation.png"); ggsave(fname,p,width=8,height=4.5,dpi=150); message("Saved: ",fname) }
  print(p); invisible(df)
}

# --------------------------------------------------------------------------
# Plot: Route comparison
# --------------------------------------------------------------------------
plot_route_comparison <- function(drug_key, dose=NULL, t_end=24, save_png=FALSE) {
  drug   <- drug_params[[drug_key]]
  dose   <- if (is.null(dose)) drug$doses[2] else dose
  routes <- c("oral","iv_bolus","iv_infusion")
  labels <- c("Oral","IV bolus","IV infusion (1 h)")
  df <- lapply(seq_along(routes), function(i) {
    res <- simulate_pbpk(drug, dose, t_end, route=routes[i])
    res %>% mutate(route_label=paste0(labels[i],"  |  Cmax=",round(res$Cmax[1],1),
                                       "  |  AUC=",round(res$AUC[1],0)))
  }) %>% bind_rows()
  rcol <- setNames(COL_DOSES[seq_along(routes)], unique(df$route_label))
  p <- ggplot(df, aes(time, C_ng_mL, color=route_label)) +
    geom_line(linewidth=1.1) + scale_color_manual(values=rcol) +
    labs(x="Time (h)", y="Plasma concentration (ng/mL)",
         title=paste0(drug$name,"  ",dose," ",drug$dose_unit," — Route comparison"),
         subtitle="Two-compartment PBPK model", color=NULL) + THEME_CLN
  if (save_png) { fname <- paste0("PBPK_",drug_key,"_routes.png"); ggsave(fname,p,width=8,height=4.5,dpi=150); message("Saved: ",fname) }
  print(p); invisible(df)
}

# =============================================================================
# Run all analyses
# =============================================================================
cat("=================================================================\n")
cat("  PBPK Two-Compartment Model  |  Krima R. Patel, Pharm.D.\n")
cat("  University at Buffalo, SUNY\n")
cat("=================================================================\n\n")

cat("── [1] Midazolam oral dose escalation ──\n")
plot_dose_escalation("midazolam", route="oral",  t_end=12, save_png=TRUE)

cat("\n── [2] Midazolam DDI: ketoconazole ──\n")
plot_ddi("midazolam", dose=5, perpetrator="ketoconazole", t_end=12, save_png=TRUE)

cat("\n── [3] Midazolam route comparison ──\n")
plot_route_comparison("midazolam", dose=5, t_end=12, save_png=TRUE)

cat("\n── [4] Warfarin DDI: amiodarone ──\n")
plot_ddi("warfarin", dose=5, perpetrator="amiodarone", t_end=72, save_png=TRUE)

cat("\n── [5] Caffeine DDI: fluvoxamine ──\n")
plot_ddi("caffeine", dose=100, perpetrator="fluvoxamine", t_end=24, save_png=TRUE)

cat("\n── All analyses complete. PNG files saved.\n")
