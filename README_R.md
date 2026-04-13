# PBPK Two-Compartment Model — R Analysis
**Author:** Krima R. Patel, Pharm.D. | M.S. Pharmaceutical Sciences Candidate  
**Institution:** University at Buffalo, SUNY  
**Domain:** DMPK · PBPK Modeling · Drug-Drug Interaction Assessment

---

## Overview

This R project implements a **two-compartment physiologically-based pharmacokinetic (PBPK) model** to simulate plasma drug concentration-time profiles and assess drug-drug interactions (DDI). The modeling approach directly mirrors workflows used in DMPK groups employing tools such as SimCYP and PK-Sim.

### Key features

| Feature | Detail |
|---|---|
| **ODE solver** | `deSolve::lsoda` (stiff ODE solver, same class as NONMEM/SimCYP numerical engines) |
| **Routes** | Oral (first-order absorption), IV bolus, IV infusion |
| **DDI mechanisms** | CYP3A4/1A2/2C9 inhibition and induction |
| **PK metrics** | AUC (trapezoidal), Cmax, Tmax, terminal t½ (log-linear fit) |
| **Visualization** | `ggplot2` + `patchwork` publication-quality figures |

---

## Model structure

### ODE system

```
Oral dosing:
  dA_gut/dt  = -ka × A_gut
  dA_c/dt    = F·ka·A_gut − (CL_eff/Vc)·A_c − (Q/Vc)·A_c + (Q/Vp)·A_p
  dA_p/dt    = (Q/Vc)·A_c − (Q/Vp)·A_p

C_plasma(t) = A_c(t) / Vc      [mg/L → ×1000 → ng/mL]
```

### DDI model (competitive inhibition)
```
CL_eff = CL × (1 − fm × R_inhibition)
```

### Drug library

| Drug | Enzyme | F | CL (L/h) | Vc (L) | t½ (approx) |
|---|---|---|---|---|---|
| Midazolam | CYP3A4 | 0.36 | 28.0 | 20 | ~1.8 h |
| Caffeine | CYP1A2 | 1.00 | 1.9 | 35 | ~4.5 h |
| Warfarin | CYP2C9 | 0.99 | 0.19 | 8 | ~40 h |

---

## Installation

```r
install.packages(c("deSolve", "ggplot2", "dplyr", "tidyr", "patchwork"))
```

Then run the full analysis:
```r
source("pbpk_model.R")
```

---

## What each analysis generates

| Analysis | Drug | Scenario | Output |
|---|---|---|---|
| 1 | Midazolam | Dose escalation (oral) | `PBPK_midazolam_dose_escalation.png` |
| 2 | Midazolam | DDI: ketoconazole (strong CYP3A4 inh.) | `PBPK_midazolam_ketoconazole_DDI.png` |
| 3 | Midazolam | Route comparison (oral vs IV) | `PBPK_midazolam_routes.png` |
| 4 | Warfarin | DDI: amiodarone (CYP2C9 inh.) | `PBPK_warfarin_amiodarone_DDI.png` |
| 5 | Caffeine | DDI: fluvoxamine (CYP1A2 inh.) | `PBPK_caffeine_fluvoxamine_DDI.png` |

---

## Example output (console)
```
── PK Summary ───────────────────────
  Scenario                           AUC (ng·h/mL)  Cmax (ng/mL)  Tmax (h)  t½ (h)
  Baseline                                    48.2          24.3      0.51    1.83
  + Strong CYP3A4 inhibitor (keto.)          471.1         163.7      0.88    8.24

AUC ratio = 9.77x  (FDA threshold: >2x = significant DDI)
```

---

## Relevance to DMPK / PBPK role

This project directly demonstrates model-informed drug development (MIDD) competencies:

- **Human dose prediction** — linear PK dose escalation across clinical dose range
- **DDI risk flagging** — AUC ratios benchmarked against FDA 2× threshold
- **Route optimization** — bioavailability comparison across administration routes
- **Scientific communication** — publication-quality ggplot2 figures with PK metric overlays

---

## References

1. Rowland M, Tozer TN. *Clinical Pharmacokinetics and Pharmacodynamics*, 4th ed. (2011)
2. Soetaert K et al. Solving differential equations in R. *R Journal* (2010)
3. FDA Guidance: In Vitro Drug Interaction Studies (2020)
4. Jamei M et al. The Simcyp population-based ADME simulator. *Expert Opin Drug Metab Toxicol* (2009)
