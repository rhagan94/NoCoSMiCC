---
title: "NoCoSMiCC"
layout: splash
permalink: /
header:
  overlay_color: "#0b1020"
  overlay_filter: "0.55"
  actions:
    - label: "View code"
      url: "https://github.com/rhagan94/NoCoSMiCC"
    - label: "Methods"
      url: "/NoCoSMiCC/methods/"
excerpt: "Non-coding Somatic Mutations in Colorectal Cancer — CNV-aware background mutation-rate modeling and non-coding driver discovery."
feature_row:
  - image_path: /NoCoSMiCC/assets/img/icon-model.svg
    alt: "Modeling"
    title: "BMR modeling"
    excerpt: "Genome-wide background mutation-rate prediction over fixed bins with rigorous evaluation."
  - image_path: /NoCoSMiCC/assets/img/icon-control.svg
    alt: "Confounders"
    title: "Confounder control"
    excerpt: "Explicit tracking of CNV and other covariates to reduce artefactual driver calls."
  - image_path: /NoCoSMiCC/assets/img/icon-discovery.svg
    alt: "Discovery"
    title: "Driver discovery"
    excerpt: "Mutation-enrichment testing of candidate regulatory elements with prioritisation for follow-up."
---

## What this project does
- Builds genome-wide background mutation-rate models over fixed bins.
- Tests candidate regulatory regions for mutation enrichment beyond expectation.
- Tracks confounders (e.g., CNV) to avoid artefactual “drivers”.


