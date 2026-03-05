---
layout: archive
title: "CV"
permalink: /cv/
author_profile: true
redirect_from:
  - /resume
---

{% include base_path %}

Education
======
* **Ph.D. Candidate in Genomics**, University of Chinese Academy of Sciences (UCAS), 2023.09 - Present
  * Awards: "Three Good Student" of UCAS (2025)
* **B.S. in Bioinformatics**, Hunan Agricultural University (HNAU), 2019.09 - 2023.07
  * Awards: National Scholarship (2022), Outstanding Graduate of Hunan Province (2023), "Three Good Student" of HNAU (2020, 2021)

Work Experience
======
* **Joint Training Student**, BGI Research (Shenzhen) - Guangdong-Hong Kong-Macao Greater Bay Area Branch, 2022.07 - Present
  * Responsible for building data preprocessing and QC standardization workflows. Developed and maintained Python/R automation scripts. Used WDL and Snakemake to construct transferable pipelines, increasing analysis efficiency by 50% and supporting the data mining of 50,000+ NIPT samples.
  * Constructed multi-omics integration workflows for molecular subtyping and biomarker discovery in complex diseases.
* **Intern**, BGI Genomics (Shenzhen) - Maternal and Child Health Division, 2021.07 - 2021.09
  * Participated in NIPT product development, learned the biological characteristics and detection workflows of cfDNA, and assisted in sample data QC and experimental record organization.

Research Projects
======
* **Multi-omics Research on Human Stress Response in High-Altitude Environments** (Core Analyst), 2022.07 - 2024.03
  * Tested various cfRNA upstream analysis tools and built customized analysis pipelines. Conducted differential analysis and pathway enrichment to parse dynamic changes of human cell-free nucleic acids in high-altitude environments.
  * Led the construction of standard analysis workflows for 6 omics (transcriptomics, proteomics, phosphoproteomics, metabolomics, etc.), unifying data processing standards.
* **Placental Multi-omics** (Core Analyst), 2024.03 - Present
  * Used MOFA to integrate multi-omics data for subtyping pregnancy diseases and constructing molecular characteristic profiles.
  * Designed and implemented xQTL cross-omics integration methods to parse pathogenesis from the perspective of genetic variation.
* **Impact of Maternal Nutrition and Metabolism on Pregnancy Diseases and Outcomes** (Core Analyst), 2025.08 - Present
  * Constructed a standardized genomic analysis workflow for ultra-low depth (< 0.1x) cfDNA, integrating fragmentomics features to build disease risk prediction models.
* **Biomarker Screening Method Based on Cross-Omics Scoring Profile**, 2024.07 - Present
  * Constructed a novel multi-dimensional scoring profile integrating similarity, statistical significance, and functional chain-specific scores for cross-omics biomarker screening.

Skills
======
* **Programming Languages**: Python, R, Shell/Bash
* **Computing Environments**: Linux System Administration, Server Cluster Usage
* **Bioinformatics Analysis**: WDL/Snakemake Workflow Development & Deployment, High-throughput Data Processing, Multi-omics Integration, Biological Data Visualization
* **Languages**: English (CET-6, good reading, writing, and communication skills)

Academic Achievements
======
* **Publications**: 2 SCI Papers (1 First Author)
* **Software Copyrights**: 4
* **Invention Patents**: Participated in 2, 1 pending

Publications
======
  <ul>{% for post in site.publications reversed %}
    {% include archive-single-cv.html %}
  {% endfor %}</ul>
  
Service and Leadership
======
* Disciplinary Inspection Committee Member, 516 Party Branch, UCAS (2023.09 - 2024.03)
* Deputy Director, Office of the Youth Volunteer Association, Youth League Committee of HNAU (2020.09 - 2021.09)
* Secretary of the Youth League Branch, Bio-Sci 19-1 & Bio-Info 19-3 (2019.09 - 2021.09)