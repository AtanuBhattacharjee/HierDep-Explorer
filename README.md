# HierDep Explorer

Companion Shiny application to *Longitudinal Mental Health and Neurodevelopmental Studies: Guide to Dealing with Non-independent Data Using Hierarchical Dependence* (Institute of Psychiatry, Psychology & Neuroscience, King's College London).

Pooled multi-site data carry three layers of dependence: participants within sites, measures within participants, and visits within participants. The app estimates each layer, shows what it does to inference, and gives the decision that follows at every step, from planning a study to testing its hypothesis.

## What the app does

| Tab | Content |
|---|---|
| Overview | The ABIDE I eye-status case and the main results by layer |
| Decision guide | Objective and hypothesis; Stage 1 planning (power, sites needed, decision path, analysis plan, script); Stage 2 analysis (decisions from the data, test of the hypothesis, planned against observed) |
| Data | Simulated ABIDE-format demo, the public ABIDE I phenotypic file, or an uploaded CSV |
| Site layer | Site ICC, design effect, effective sample size, site effects |
| Exposure | Level of each variable, site-level exposure test, real-data null experiment |
| Groups | Four analyses, random slope, prediction interval for a new site |
| Site make-up | Instrument crossed with site; institution against cohort |
| Measures | Within-site, between-site and pooled correlations |
| Time | Correlation between visits; longitudinal models |
| Planning, Few sites | Design effects for new studies; conversion of between-centre η² to the ICC |
| Report | Methods and results text and result tables |

Every analysis tab shows the inspection on the left and the decision it leads to on the right.

## Run locally

```r
install.packages(c("shiny", "bslib", "ggplot2", "DT", "lme4", "lmerTest", "nlme", "pbkrtest", "Matrix"))
shiny::runGitHub("HierDep-Explorer", "<github-user>")
```

## Data

The demo data are simulated in the format of the ABIDE I phenotypic file and contain no participant data. The ABIDE I option reads the public phenotypic file of the ABIDE Preprocessed release. Uploaded files stay in the session and are not stored.

## Licence

MIT
