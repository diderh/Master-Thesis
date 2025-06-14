# Identifying the single and combined effects of ALAN and signal crayfish on benthic communities
---

## Table of Contents

- [Project Overview](#project-overview)
- [Features](#features)
- [Getting Started](#getting-started)
- [Workflow Steps](#workflow-steps)
- [Outputs](#outputs)
- [Conclusion](#conclusion)

---

# Project Overview
This R script analyzes the effects of environmental stressors (such as light pollution and crayfish presence) on macroinvertebrate functional diversity in aquatic ecosystems. The workflow processes benthic macroinvertebrate data, assigns trait information, calculates multiple functional diversity indices, and performs advanced statistical analyses to assess impacts of treatments across time and positions.

---

# Features
  - Data reshaping and cleaning for macroinvertebrate abundance tables.
  - Taxonomic and trait assignment using biomonitoR and reference data.
  - Calculation of multiple functional diversity indices: richness, dispersion, evenness, redundancy, Rao’s Q, and community weighted means (CWMs).
  - Multivariate analyses (NMDS, PERMANOVA, PERMDISP, SIMPER) to test for effects of treatment, time, and position.
  - Advanced statistical modeling using mixed-effects and Bayesian models for functional indices.
  - Extensive plotting and export of results for publication-quality figures.

---

# Getting Started
  - Required R packages: biomonitoR, vegan, ade4, StatMatch, ggplot2, tidyverse, readxl, brms, lme4, and others.
  - Input data: Macroinvertebrate abundance table in Excel (e.g., Benthos_final_updatedfile_Tanypodinae.xlsx).
  - Instructions are provided for package installation and data import.

---

# Workflow Steps
1. **Data Import & Reshaping:** Load macroinvertebrate abundance data, convert to long format, and create sample identifiers.
2. **Trait Assignment:** Merge sample data with reference trait tables, aggregate taxa, and assign trait values.
3. **Functional Indices Calculation:** Compute indices (richness, dispersion, evenness, redundancy, CWMs) using fuzzy-coded trait data.
4. **Statistical Analysis:**
   - NMDS ordination to visualize functional trait patterns.
   - PERMANOVA and PERMDISP to test for differences by treatment/time/position.
   - SIMPER to identify traits contributing most to group differences.
   - Bayesian/Mixed-effect models to analyze functional indices (richness, dispersion, evenness, Rao’s Q), with diagnostic and conditional effect plots.
5. **Visualization & Export:** Generate and save figures for NMDS, boxplots, posterior checks, and model parameter distributions.

---

# Outputs
  - Processed data tables for functional indices and trait assignments.
  - NMDS ordination plots, PERMDISP boxplots, density plots, and conditional effect visualizations.
  - Model summaries and posterior distributions for each functional diversity metric.
  - PNG images of all major plots for reporting.

---

# Conclusion
This study provides a comprehensive assessment of how environmental stressors, specifically artificial light at night (ALAN) and crayfish presence (alone and in combination), impact the functional diversity of aquatic macroinvertebrate communities over time and across stream positions (upstream, downstream).

## Findings:

   - **Functional Diversity Changes:** The script calculates multiple indices (richness, dispersion, evenness, redundancy, Rao’s Q) and finds measurable shifts in these metrics in response to different treatments (ALAN, crayfish, combined, control), at different times (before and after treatment), and positions (upstream, downstream).
   - **Statistical Significance:** Through PERMANOVA and Bayesian mixed-effects modeling, the analysis identifies statistically significant effects and interactions between treatment, time point, and position. The models indicate that both individual and combined stressors can alter the structure and function of macroinvertebrate communities.
   - **Trait Contributions:** SIMPER analyses pinpoint which functional traits contribute most to the observed differences among treatments, providing insight into which ecological roles or adaptations are most sensitive to environmental stressors.
   - **Community Shifts Visualized:** NMDS ordination and boxplots visually confirm community shifts in functional composition, with clear groupings and separations corresponding to experimental manipulations.
   - **Model Diagnostics:** Posterior predictive checks and trace plots show that the statistical models fit the observed data well and that findings are robust.

Overall, the results demonstrate that both artificial light and invasive crayfish (individually and together) can significantly alter the functional diversity and trait structure of stream macroinvertebrate communities. These changes are context-dependent, varying by stream location and over time, underscoring the complex ecological impacts of multiple, interacting anthropogenic stressors.

---

