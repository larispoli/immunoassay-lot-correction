# Immunoassay Lot-to-Lot Correction Pipeline

A diagnostic and correction tool for assessing and adjusting plate-to-plate and lot-to-lot variation in immunoassay results. Designed for retrospective application to reported concentrations in longitudinal datasets where raw sample-level signal data may not be available.


*NOTE from Louisa: I apologize for the long length of this document. Please read the entire thing at least once, as it should provide clarity about what happens with data when lot correction is utilized.*
**This is under beta testing at moment as need others to check my code and implementation of logic for handling assay variation**

## Contents

- [Overview](#overview)
- [When to Use This Pipeline](#when-to-use-this-pipeline)
- [When to Use ELISAtools Instead](#when-to-use-elisatools-instead)
- [How the Correction Works](#how-the-correction-works)
- [Correction Tiers](#correction-tiers)
- [Dataset Considerations](#dataset-considerations)
- [Selecting the Reference Lot](#selecting-the-reference-lot)
- [Reference Samples](#reference-samples)
- [Requirements](#requirements)
- [Repository Contents](#repository-contents)
- [Quick Start](#quick-start)
- [Interpreting Key Outputs](#interpreting-key-outputs)
- [Methods Summary for Manuscript Use](#methods-summary-for-manuscript-use)
- [References](#references)
                            
## **Overview**
                            
Immunoassays used in longitudinal studies are subject to variation from multiple sources, including manufacturing differences between reagent lots, batch variability within lots, and changes in assay components over time. This variation can introduce systematic bias into reported concentrations, affecting downstream analyses. While plate and lot can be included as covariates in mixed modeling analyses to help correct for the variation, for other analyses, such as cycle detection, ROC analysis, and time series modeling, the covariate option is not available, and a correction for the lot-to-lot variation must occur before the data is used.
                            
This pipeline provides a structured, reproducible workflow to assess the magnitude and nature of this variation, apply defensible corrections where appropriate, and validate the impact of those corrections on actual sample data. It produces an interactive R Notebook with diagnostic tables, figures, and guided decision points.
                            
## **When to Use This Pipeline**
                            
This pipeline is appropriate when:
                              
- You have reported concentrations from immunoassays run across multiple plates and/or reagent lots.
- Raw sample-level signal data (individual well ODs for each unknown sample) is not available or not practical to retrieve retrospectively.
- Standard curve data (concentration and signal values for each standard point) can be recovered from plate reports, instrument software, or laboratory records.
- You need corrected values on a consistent scale for analyses such as ROC curves, threshold-based detection, cross-analyte comparisons, etc.
                            
## **When to Use ELISAtools Instead**
                            
The ELISAtools package (Feng et al. 2019) is more appropriate when:
                              
- You have raw plate reader output files (96-well format with individual well ODs for all standards, controls, and unknowns).
- You want to recalculate sample concentrations from the adjusted standard curve rather than applying a multiplicative correction to already-reported values.
- Your data is prospective, and you can integrate ELISAtools into your routine data processing workflow.
  
This pipeline adapts the core methodology of Feng et al. (2019) -- simultaneous curve fitting with shared shape parameters -- but applies the resulting shift factor to reported concentrations rather than recalculating from raw signals. This makes it suitable for retrospective datasets where raw signals are unavailable.
                            
## **How the Correction Works**
                            
### **The Simultaneous Curve Fit**
                            
All standard curves in the dataset are fitted simultaneously using the `drc` package in R. The key constraint: all curves share the same shape parameters (upper asymptote, lower asymptote, slope, and for 5PL, asymmetry). Only the curve midpoint (ED50) is allowed to vary per plate.
                            
This enforces the assumption that variation between plates manifests as a horizontal shift on the concentration axis -- the curve slides left or right but retains the same shape. A sample producing a given signal would be assigned a higher or lower concentration depending on which plate it was run on, purely because of where that plate's curve sits on the x-axis.

### **The Shift Factor**

The shift factor for each plate is the ratio of the reference lot's median midpoint to that plate's individual midpoint:

    Shift Factor = Reference Lot Median Midpoint / Plate Midpoint

A shift factor of 1.0 means the plate's curve is in the same position as the reference lot's typical curve. A shift factor of 1.3 means the reference lot would report a concentration 1.3x higher than this plate for the same signal.

The corrected concentration is: Corrected Value = Reported Value x Shift Factor

Every plate receives its own shift factor, including plates in the reference lot. The reference lot's median midpoint defines the scale; individual plates within the reference lot are corrected toward their own lot's median just like plates in any other lot.

### **Outlier Detection**

Plates whose standard curves do not conform to the shared shape are identified using the interquartile range (IQR) method applied to per-plate root mean square error (RMSE) values. Plates with RMSE exceeding Q3 + 1.5 x IQR are classified as outliers (Tukey 1977).

Flagged plates are excluded from the calculation of the reference lot's median midpoint so they do not influence the baseline. They still receive their own per-plate correction factor, but their corrected values carry additional uncertainty. Users should identify which samples were run on flagged plates.
                            
#### *Single-Lot Datasets:* 

For datasets with only one reagent lot, the pipeline functions as a plate-to-plate correction tool. All plates are corrected toward the lot's median midpoint. The workflow is identical; there is simply no between-lot component.

## **Correction Tiers**

### **Tier 1 -- Standard Curve Shift Factor**

The primary correction method. Uses the simultaneous curve fit described above. Available for every plate in the dataset because every plate has standard curve data (a requirement for inclusion in the Assay Data table).

Tier 1 is based on the methodology described in Feng et al. (2019), adapted for application to reported concentrations rather than recalculation from raw signals.
  
  **When Tier 1 May Not Be Appropriate**

  *The simultaneous fit shares all curve shape parameters across plates. This works when lot-to-lot variation is a shift in assay sensitivity (the curve moves left or right) without a change in curve shape. When a fundamental assay change has occurred -- such as a new antibody batch with different binding characteristics, a change in cross-reactivity, or a switch in assay format -- the curves may differ in shape, not just position.*


### **Tier 2 -- Reference Sample Based Correction**

An independent, measurement-based correction. A biological reference sample (e.g., pooled sample from a relevant species) measured on multiple plates provides a direct correction factor:

    Tier 2 Factor = Reference Lot Mean Concentration / Plate Concentration (same sample)

This follows the bridging sample normalization approach described by Henson et al. (2022). Tier 2 does not depend on curve shape assumptions. It directly measures how much a plate's reported concentration for a known sample differs from the reference lot's value.

When both Tier 1 and Tier 2 are available, an agreement between them provides the strongest evidence that the correction is appropriate. Divergence suggests the lot effect involves more than a simple curve shift and warrants further investigation.

  **Tier 2 Implications**

  *If the selected reference lot has no biological reference sample data, Tier 2 will not be available. The pipeline alerts you to this after selection. If Tier 2 validation is important for your dataset, consider selecting a lot that has reference sample data.*

### **Choosing Between Tiers**

When both tiers agree, either is defensible; Tier 1 provides full coverage while Tier 2 validates it.

When Tier 1 shows increased variability in the reference sample consistency check (Phase 2.2): Tier 2 is the recommended correction method, as the shared-shape assumption does not adequately describe the variation.

When Tier 2 is unavailable (no reference samples or reference samples do not span all lots), Tier 1 is the only option. The Phase 2.2 diagnostics indicate how much confidence to place in it.


## **Dataset Considerations**

### The Five-parameter Logistic Curve Fitted (5PL) Datasets

The 5PL includes an asymmetry parameter (f) that describes how the upper and lower portions of the curve differ in steepness. In the simultaneous fit, this parameter is shared across all plates.

If the asymmetry genuinely differs across lots (e.g., because antibody binding kinetics changed between manufacturing batches), forcing a single shared value creates a compromise that does not accurately describe any lot. The midpoint can only shift the curve horizontally; it cannot compensate for a shape change. In this situation:

- The reference sample consistency check (Phase 2.2) will show that the simultaneous model increased variability compared to the original reported values.
- The correction factors from Tier 1 will be less reliable.
- Tier 2 (bridging sample correction) should be used as the primary correction method, as it does not depend on curve shape assumptions.

**Important distinction:** This applies when the antibody specificity remains the same but sensitivity or binding kinetics vary between batches. If the antibody was replaced with one that has a fundamentally different cross-reactivity profile (as in Wilson et al. 2021, where 11-deoxycortisol cross-reactivity changed from 15% to 58%), the assay is effectively measuring a different analyte. In that case, data from the two antibodies cannot be combined through correction — they represent different biological measurements and should be treated as separate datasets.

### The Four-parameter Logistic Curve Fitted (4PL) Datasets

For datasets where all plates use 4PL (four-parameter logistic), the shape is symmetric by definition. The simultaneous fit constrains slope and asymptotes to be shared. If these truly are consistent across lots (same antibody, same assay format), Tier 1 will work well. If the reference sample consistency check shows increased variability even with 4PL, the variation may involve changes in slope or asymptotes rather than just midpoint position.

### Mixed 4PL and 5PL Datasets

When some plates were originally fitted with 4PL and others with 5PL, the pipeline tests both shared models and selects the one with lower overall RMSE. In most cases, the 5PL will be selected because it can accommodate slight asymmetry in 4PL plates (the asymmetry parameter will be estimated near 1.0 for symmetric curves) while also fitting the genuinely asymmetric 5PL plates.

The pipeline reports the median difference between original reported concentrations and simultaneous model predictions for each model type separately. If the differences are small (under 5%), the original curve model choice did not substantially affect the reported values, and no action is needed. If the differences are large, and raw plate reader data is available, rerunning the affected plates with a consistent curve model would produce more accurate starting concentrations before applying the shift factor correction.

## Selecting the Reference Lot

The reference lot defines the baseline scale for all corrected values. Both Tier 1 and Tier 2 corrections are anchored to this lot.

### Items to Consider

- **Maintain consistency:** If a lot was the reference in previous analyses or publications, keep it. This ensures previously corrected values and any thresholds derived from them remain comparable.
- **Stability:** Prefer the lot with the most plates (more stable median midpoint) and the tightest clustering of midpoints (most internally consistent baseline).
- **Future direction:** If a particular lot represents the assay conditions you intend to use going forward, selecting it means future data should require minimal correction.
- **Fundamental assay changes:** If a complete antibody replacement or assay format change occurred, the new conditions may warrant a new reference lot. Values corrected to the new reference are not directly comparable to values corrected to the previous one. This effectively creates two separate analytical eras that should be documented (Wilson et al. 2021).
- **Scale, not relationships:** The choice of reference lot does not affect the relative correction between any two plates. It only sets the absolute scale that corrected values are expressed on.


## Reference Samples

### What Qualifies as a Reference Sample

A biological reference sample is a pooled biological sample (e.g., pooled serum or plasma from a relevant species) that is aliquoted, stored, and run on each plate alongside unknown samples. Its purpose is to track assay consistency across plates and lots.

### Bridging References

A reference sample that is measured across multiple lots provides a direct lot-to-lot comparison and enables Tier 2 correction. This is the most valuable type of reference data for longitudinal datasets.

### Cross-Reference Bridging

When transitioning between reference samples, running both the old and new reference on at least one plate establishes the ratio between them. The pipeline detects these dual-reference plates automatically and uses the ratio to extend Tier 2 coverage across lots that used different reference pools.

## Requirements

- R version 4.0 or later
- RStudio (for interactive notebook execution)
- R packages: readxl, dplyr, tidyr, ggplot2, drc, lme4, broom, patchwork, scales, knitr, kableExtra

Packages are installed automatically on first run if not already present.

## Repository Contents

| File | Description |
|------|-------------|
| `Lot_Correction_Pipeline.Rmd` | The R Notebook pipeline. Open in RStudio and run interactively. |
| `Assay_Plate_Data_Template.xlsx` | Template for plate-level metadata input. Contains four sheets: Instructions, AssayPlateData (blank template), Codebook (data dictionary), and Examples. |
| `README.md` | This file. |

## Quick Start

1. **Prepare your data.** Fill out the Assay Plate Data Template with one row per assay plate. See the Instructions and Examples sheets in the template for guidance.

2. **Open the pipeline.** Open `Lot_Correction_Pipeline.Rmd` in RStudio.

3. **Run chunks one at a time.** Prompts will appear in the Console at decision points.

4. **Review outputs.** Diagnostic tables and figures appear inline after each chunk. Review them before responding to prompts.

5. **Save.** After all chunks complete, save the file (Ctrl+S) to capture all outputs into the notebook (.nb.html).

6. **Rename.** Rename the .nb.html file to include the analyte name (e.g., `Lot_Correction_Pipeline_P4.nb.html`) before running a different analyte.

7. **Optional: Phase 6.** If you have longitudinal sample data, Phase 6 visualizes individual plots before and after correction to assess biological impact.

#### **Optional Input (for Phase 6 -- Apply & Visualize Correction):**

Sample data file (.csv or .xlsx) containing longitudinal sample concentrations linked to plate metadata. At minimum, this file must include:

| Column | Description |
|--------|-------------|
| Sample_ID | Unique identifier for the biological sample |
| Animal_ID | Individual animal identifier for longitudinal tracking |
| Plate_ID | Must match Plate_ID values in the Assay data file |
| Conc | Reported concentration; the correction factor is multiplicative, so order of operations with dilution factor does not matter |
| Coll_Date | Date of sample collection |

Additional columns (e.g., species, sex, reproductive status) are preserved but not required by the pipeline.

## Interpreting Key Outputs

### Phase 2.1 -- Simultaneous Curve Fit

**What to look for:** All plates fit the shared-shape model well (no outlier flags) = Tier 1 is appropriate. Flagged plates have standard curves that do not conform to the shared shape.

**If mixed curve models detected (4PL and 5PL):** The pipeline tests both and selects the better fit. Pay attention to the median difference by model type.

### Phase 2.2 -- Reference Sample Consistency Check

**This is the most important diagnostic.** It compares the CV of reference sample concentrations under the original individual fits versus the simultaneous shared model.

- **Predicted CV lower than Reported CV:** The simultaneous model improved consistency. Tier 1 correction is supported.
- **Predicted CV higher than Reported CV:** The simultaneous model made things worse. The shared-shape assumption does not hold. Tier 2 is the recommended correction method.

### Phase 3.3 -- Midpoint Plot

**What to look for:** Tight clustering within lots indicates consistent plates. Clear separation between lots indicates a systematic lot effect. Widely scattered midpoints within a lot indicate substantial plate-to-plate variation.

### Phase 4.3 -- Tier Comparison

**What to look for:** Points near the 1:1 line indicate agreement between correction methods. Systematic deviation suggests the curve-based and measurement-based corrections are capturing different aspects of variation.

### Phase 5.2 -- Post-Correction Validation

**What to look for:** Reference sample CVs should decrease after correction. If they increase, the correction may not be appropriate.

### Phase 6.4 -- Before & After Plots

**What to look for:** Where the orange (before) and blue (after) lines diverge, the correction changed the value. This should align with plate or lot transitions.

## Methods Summary for Manuscript Use

The following text can be adapted for the methods section of manuscripts using data corrected with this pipeline.

---

Plate-to-plate and lot-to-lot variation in immunoassay results was assessed and corrected using a simultaneous curve fitting approach adapted from Feng et al. (2019). Standard curve data (concentration and signal values for each standard point) were extracted from plate records for all assay plates contributing sample values to the dataset. Standard curves were fitted simultaneously using the drc package in R (Ritz et al. 2015), with shape parameters (upper asymptote, lower asymptote, slope) shared across all plates and only the curve midpoint (ED50) allowed to vary per plate. This enforces the assumption that variation between plates manifests as a horizontal shift on the concentration axis rather than a change in curve shape.

Plates with standard curves that did not conform to the shared shape were identified using the interquartile range (IQR) method applied to per-plate root mean square error (RMSE) values, with plates exceeding Q3 + 1.5 x IQR classified as outliers (Tukey 1977). Flagged plates were excluded from the calculation of the reference lot median midpoint but retained their own per-plate correction factor, with the understanding that corrected values from these plates carry additional uncertainty.

A reference lot was designated as the baseline, and the shift factor for each plate was calculated as the ratio of the reference lot's median midpoint to the plate's individual midpoint. The shift factor was applied as a multiplicative correction to the reported sample concentrations. Unlike the original ELISAtools implementation (Feng et al. 2019), which recalculates sample concentrations from raw signal data using the adjusted standard curve, this approach applies the correction to already-reported concentrations. This makes it suitable for retrospective datasets where raw sample-level signal data may not be available.

Where biological reference samples (pooled serum run alongside unknowns) were available, reference sample ratios were calculated as an independent correction method (Tier 2), following the bridging sample normalization approach described by Henson et al. (2022). When the same reference sample was measured across multiple lots, agreement between the curve-based (Tier 1) and reference-based (Tier 2) correction factors was evaluated as convergent evidence for the correction.

For analytes where the simultaneous curve fit did not improve reference sample consistency -- indicating that lot-to-lot variation involved changes in curve shape rather than position (e.g., due to antibody batch differences) -- the bridging sample correction (Tier 2) was used as the primary correction method, as it does not depend on curve shape assumptions (Wilson et al. 2021).

Correction effectiveness was evaluated by comparing reference sample concentration CVs before and after correction, and by visual inspection of individual longitudinal hormone plots with and without the correction applied.

For analyses using mixed-effects models, including plate as a random effect should be considered to absorb residual variation not captured by the standard curve correction (Whitcomb et al. 2010). For analyses using corrected values directly (e.g., threshold-based detection, time series), the corrected values represent the best available estimate of the true concentration on the reference lot's scale.
                            
All analyses were performed in R using the tidyverse suite for data manipulation and visualization.
                            
---
                              
## References
                              
- Feng F, Thompson MP, Thomas BE, et al. (2019). A computational solution to improve biomarker reproducibility during long-term projects. *PLoS ONE*, 14(4), e0209060. https://doi.org/10.1371/journal.pone.0209060
                            
- Henson RL, Volluz K, Saef BA, et al. (2022). A methodology for normalizing fluid biomarker concentrations across reagent lots. *Alzheimer's and Dementia*, 18(S6), e066912. https://doi.org/10.1002/alz.066912

- Ritz C, Baty F, Streibig JC, Gerhard D. (2015). Dose-response analysis using R. *PLoS ONE*, 10(12), e0146021. https://doi.org/10.1371/journal.pone.0146021

- Tukey JW. (1977). *Exploratory Data Analysis*. Reading, MA: Addison-Wesley Publishing Company. ISBN 978-0-201-07616-5.

- Whitcomb BW, Perkins NJ, Albert PS, Schisterman EF. (2010). Treatment of batch in the detection, calibration, and quantification of immunoassays in large-scale epidemiologic studies. *Epidemiology*, 21(Suppl 4), S44-S50. https://doi.org/10.1097/EDE.0b013e3181f37e32

- Wilson AE, Sergiel A, Selva N, et al. (2021). Correcting for enzyme immunoassay changes in long term monitoring studies. *MethodsX*, 8, 101212. https://doi.org/10.1016/j.mex.2021.101212
