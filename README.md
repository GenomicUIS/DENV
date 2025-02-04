# DENV

The R code **“DENGUE-serotyopes”** aims to download, process and clean up the complete genome sequences of Dengue Virus (*Orthoflavivirus denguei*) of the four serotypes (DENV1, DENV2, DENV3, DENV4) from the NCBI nucleotide database. It uses several functions to perform secure HTTP requests, search for genome IDs, download sequences in FASTA format, validate the downloaded sequences, and combine them into a single file. In addition, the code includes a process to remove gaps and undetermined nucleotides from the combined sequences. This process is crucial for obtaining high-quality genomic data, which is essential for research studies, treatment development, viral identification techniques such as RT-PCR and RT-LAMP, and dengue vaccines. The automation of these tasks allows handling large volumes of data efficiently and accurately with the use of High Performance Computing (HPC), facilitating progress in understanding and combating this disease in endemic countries such as Colombia.

The Python code **“DENGUE_Filter”** in a second verification code that uses the Biopython library to filter sequences in a FASTA file, removing those containing a fraction of ambiguous (non-ACGT) bases above a specified threshold. Thus, the is_valid_sequence function checks whether a sequence is valid according to the proportion of ambiguous bases, while the filter_fasta function processes the FASTA file, writing only valid sequences to a new file. This process is crucial in genomic and bioinformatics studies, as it guarantees the quality of the sequences used in subsequent analyses, avoiding errors and inaccurate results due to ambiguous or low-quality data. The automation of this filtering improves efficiency and accuracy in the preparation of genomic data for scientific research.

The R code **“DENGUE-maps”** performs a descriptive analysis of dengue cases in Colombia, using data from the epidemiological surveillance system Sistema Nacional de Vigilancia en Salud Pública (SIVIGILA) the National Institute of Health of Colombia. Link: https://www.ins.gov.co/Direcciones/Vigilancia/Paginas/SIVIGILA.aspx. The analysis includes the temporal trend of annual cases, distribution by gender and distribution by age group. In addition, the code generates maps of geographic distribution of cases by department, both nationally and individually for each year between 2007 and 2023. This analysis is crucial to identify patterns and trends in the incidence of Dengue Virus, which can help public health authorities design more effective prevention and control strategies. The visualization of the data through graphs and maps facilitates the interpretation and communication of the findings, allowing for better evidence-based decision making.

**References**

[1] E. A. Martínez Álvarez, “SIVIGILA, an infrastructure mobilizing diseases, policies and practices in public health surveillance,” Revista Colombiana de Sociología, vol. 39, no. 2, Jul. 2016, doi: 10.15446/rcs.v39n2.58977.

[2] E. W. Sayers et al., “GenBank 2024 Update,” Nucleic Acids Res, vol. 52, no. D1, pp. D134–D137, Jan. 2024, doi: 10.1093/nar/gkad903.

[3] R Core Team, “R: A language and environment for statistical computing,” 2021, R Foundation for Statistical Computing, Vienna, Austria.: 4.4.2. Accessed: Feb. 01, 2025. [Online]. Available: https://www.R-project.org/

[4]	RStudio Team, “RStudio: Integrated Development for R,” 2020, RStudio, PBC, Boston, MA: 2024.12.0.

[5]	Python Software Foundation, “Python Language Reference,” 3.12.2. Accessed: Feb. 01, 2025. [Online]. Available: https://docs.python.org/3/

[6]	P. J. A. Cock et al., “Biopython: freely available Python tools for computational molecular biology and bioinformatics,” Bioinformatics, vol. 25, no. 11, pp. 1422–1423, Jun. 2009, doi: 10.1093/bioinformatics/btp163.

[7]	H. Wickham, ggplot2: Elegant Graphics for Data Analysis, Second Edition. in Use R! Cham: Springer International Publishing, 2016. doi: 10.1007/978-3-319-24277-4.

[8]	H. Wickham, R. François, L. Henry, K. Müller, and D. Vaughan, “dplyr: A Grammar of Data Manipulation,” 2023, 1.1.4. Accessed: Feb. 01, 2025. [Online]. Available: https://github.com/tidyverse/dplyr

[9]	E. Pebesma, “Simple Features for R: Standardized Support for Spatial Vector Data,” R J, vol. 10, no. 1, pp. 439–446, 2018.

[10]	D. Dunnington, “ggspatial: Spatial Data Framework for ggplot2,” 2023. Accessed: Feb. 01, 2025. [Online]. Available: https://github.com/paleolimbot/ggspatial

[11]	H. Wickham, D. Vaughan, and M. Girlich, “tidyr: Tidy Messy Data,” 2024, 1.3.1. Accessed: Feb. 01, 2025. [Online]. Available: https://tidyr.tidyverse.org


