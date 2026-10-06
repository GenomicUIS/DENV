# SIVIGILA Dengue Pipeline

R workflow for harmonizing, auditing, deduplicating and visualizing historical Colombian SIVIGILA notifications for **dengue (210)**, **severe dengue (220)** and **dengue mortality (580)**, covering **2007–2025**.

Developed at the **Laboratorio de Genómica Celular Aplicada (LGCA), Universidad Industrial de Santander (UIS), Bucaramanga, Colombia**, in the context of a Master's thesis in Biology. Source code is maintained in the [GenomicUIS/DENV repository](https://github.com/GenomicUIS/DENV/tree/main/code_SIVIGILA) [1].

The workflow connects historical data preparation with traceable epidemiological summaries and high-resolution figures.

## Contents

- [Features](#features)
- [Data sources](#data-sources)
- [Requirements](#requirements)
- [Project structure](#project-structure)
- [Workflow](#workflow)
- [Script reference](#script-reference)
- [Configuration](#configuration)
- [How to run](#how-to-run)
- [Outputs](#outputs)
- [Data quality and interpretation](#data-quality-and-interpretation)
- [Loading processed data](#loading-processed-data)
- [Reproducibility and software citations](#reproducibility-and-software-citations)
- [How to cite this project](#how-to-cite-this-project)
- [License and attribution](#license-and-attribution)
- [References](#references)

## Features

- Historical harmonization of annual column names, data types and selected SIVIGILA codes.
- Source traceability through event, filename, file year and Excel row fields.
- Annual canonical and analytical partitions for event 210.
- Auditing of missing values, invalid dates, age units, hospitalization information and chronological consistency.
- Within-event exclusion of analytical redundancy, with retained representatives and excluded records documented in audit tables.
- Temporal, demographic, social, geographical and care-related summaries.
- Event-specific figures exported as 600 dpi TIFF files, with PDF versions, PNG previews, figure indexes and saved R plot objects.
- Comparative social and care figures across events 210, 220 and 580.
- DANE MGN 2025 cartography for **32 departments and Bogotá D. C.**, including an enlarged representation of San Andrés, Providencia and Santa Catalina [6].

## Data sources

Download annual Excel files manually from the [INS SIVIGILA microdata portal](https://portalsivigila.ins.gov.co/buscador) [2]. The current configuration requires one `.xls` or `.xlsx` file for each year from 2007 to 2025 in each of the three event folders.

The data dictionaries and notification forms referenced by the scripts are listed in references [3–5]. A current dictionary is a documentary aid: historical codes and case definitions should also be checked against the documentation applicable to each year.

Departmental geometries are retrieved from the **DANE MGN 2025 Departamento layer, ID 319** [6]. They are downloaded and cached when the required local cartography is absent. Archive files alone are not the cache read by `departmental_maps.R`; extracted files must be placed under the configured data root using the filenames expected by that script.

Record the download date and the exact input files used in each analysis.

## Requirements

Use R [7] with a compatible installation of the packages below. Record the R version and package versions from the successful analysis environment.

| Package | Role | Reference |
| --- | --- | --- |
| `dplyr` | Data manipulation and grouped summaries | [8] |
| `stringr` | String handling | [9] |
| `fs` | Filesystem operations | [10] |
| `purrr` | Iteration and functional programming | [11] |
| `readxl` | Excel import | [12] |
| `janitor` | Column-name cleaning | [13] |
| `tidyr` | Reshaping tables | [14] |
| `ggplot2` | Statistical graphics | [15,23] |
| `scales` | Axis labels, number formats and visual scales | [16] |
| `patchwork` | Multi-panel plot composition | [17] |
| `sf` | Spatial vector data and coordinate transformations | [18,24,25] |
| `lwgeom` | Geodetic support used by the spatial workflow | [19] |
| `lubridate` | Date parsing | [20,26] |
| `readr` | Delimited-table exports | [21] |
| `tibble` | Explicit creation of tabular objects | [22] |

`tibble` is called directly by the scripts even though it is also installed as a dependency of other packages. Packages distributed with R, including `stats`, `utils`, `grid` and `grDevices`, are described in [7].

Install the external R packages once:

```r
install.packages(c(
  "dplyr", "stringr", "fs", "purrr", "readxl", "janitor",
  "tidyr", "ggplot2", "scales", "patchwork", "sf", "lwgeom",
  "lubridate", "readr", "tibble"
))
```

The code uses `dplyr::pick()` and inside-legend options in `ggplot2`. These require `dplyr` 1.1.0 or later and `ggplot2` 3.5.0 or later, respectively.

When compiling spatial packages from source, additional system libraries may be needed, including GDAL, GEOS, PROJ, SQLite and the dependencies of `units`. See the package documentation [18,19]. The figure exports also require working TIFF/PNG graphics devices and Cairo PDF support:

```r
capabilities(c("tiff", "png", "cairo"))
```

**Memory:** `09_prepare_dengue.R` and the partition-writing stage of `11_deduplicate_dengue.R` process annual tables. However, `02_initial_diagnosis.R` loads all three historical event databases, global duplicate checks collect keys across years, and mortality/severe-dengue preparation consolidates historical data. RAM needs depend on the input size, the number of columns and these intermediate objects. This project used a computer with 12 GB of RAM, which was sufficient to run all the analyses. Please check your computer and its performance before running the program.

## Project structure

Keep the scripts together in `code_SIVIGILA/`, matching the repository. With the default configuration, the data root is its parent directory. An alternative data root can be supplied through `SIVIGILA_PATH`.

| Location relative to the data root | Content |
| --- | --- |
| `DENGUE/` | Event 210: one annual Excel file for 2007–2025 |
| `DENGUE GRAVE/` | Event 220: one annual Excel file for 2007–2025 |
| `MORTALIDAD DENGUE/` | Event 580: one annual Excel file for 2007–2025 |
| `code_SIVIGILA/` | Pipeline scripts and this README, in the default layout |
| `datos_externos/cartografia/` | Downloaded geometries and cartography caches |
| `resultados/` | Event-specific data, reports and figures; comparative outputs |
| `salidas/diagnostico_inicial/` | Initial audit tables |
| `salidas/auditoria_integral/` | Saved figure evidence objects |

`run_all.R` writes logs to `../salidas/log_run_all/` relative to the script directory. This location is independent of a separately configured `SIVIGILA_PATH`.

Example input filename: `DENGUE/2018.xlsx`. The filename must contain a four-digit year. The inventory validation rejects missing years, additional years and multiple files for the same year.

## Workflow

1. Resolve paths, validate dependencies and inventory the annual files.
2. Define import functions and produce the initial diagnostic reports.
3. Prepare, deduplicate and plot mortality notifications (scripts 03–05).
4. Prepare, deduplicate and plot severe-dengue notifications (scripts 06–08).
5. Prepare annual dengue partitions, document their variables, deduplicate and plot them (scripts 09–12).
6. Generate comparative social and care figures (scripts 13–14).

## Script reference

### Orchestration and shared utilities

| Script | Purpose |
| --- | --- |
| `000_run.R` | Sources `run_all.R` and calls `correr_pipeline()`. |
| `run_all.R` | Runs scripts 00–14 sequentially in a shared environment; supports `desde`, `hasta` and `detener_en_error`. Writes a timestamped log and, after completing the selected loop, `last_run.csv`. |
| `00_config.R` | Resolves the data root, checks required folders/packages and the 2007–2025 inventory, defines audit helpers and sources `99_utils.R`. |
| `01_load_databases.R` | Defines `cargar_evento()` and `cargar_todas_las_bases()`; imports Excel as text, normalizes names and adds provenance. Defining these functions does not itself load every database. |
| `02_initial_diagnosis.R` | Loads the three historical databases and audits row counts, exact duplicates, missing values and dates. |
| `99_utils.R` | Shared helpers for missingness, dates, codes, record keys, duplicate marking and age groups. |
| `figures_quality.R` | Shared numerical formatting and `guardar_evidencia_figura()`, which saves PDF versions and R plot objects. |
| `departmental_maps.R` | Downloads/caches DANE MGN 2025 geometries and draws departmental maps with scale bars, a north indicator and enlarged island representations. |
| `load_processed_dengue.R` | Optional interactive helper; defines `cargar_dengue_procesado()`. It allows you to load and run only the dengue database, without having to re-execute the entire code. It is outside the execution order in `run_all.R`. |

### Event processing

| Script | Purpose |
| --- | --- |
| `03_prepare_dengue_mortality.R` | Imports event 580, parses dates, derives age/geographical/care variables and adds quality flags. |
| `04_deduplicate_dengue_mortality.R` | Excludes redundant mortality records using the compared analytical fields; saves a deduplicated base and a plotting subset with `con_fin == "2"`. |
| `05_generate_dengue_mortality_graphics.R` | Defines 11 mortality figures covering time, demography, geography, care, social characteristics, cause codes, vulnerability and quality. |
| `06_prepare_severe_dengue.R` | Harmonizes event 220, derives variables, marks duplicate candidates and exports canonical/analytical tables, dictionaries and quality reports. |
| `07_deduplicate_severe_dengue.R` | Excludes records matching analytical content, retains a representative and exports the exclusion audit. |
| `08_generate_severe_dengue_graphics.R` | Defines 11 severe-dengue figures, including hospitalization and the proportion with a fatal recorded outcome. |
| `09_prepare_dengue.R` | Processes event 210 by year; aligns schemas, derives variables, saves keys for global duplicate detection and writes annual canonical/analytical partitions. |
| `10_document_dengue.R` | Consolidates variable profiles, quality reports, distributions, the data dictionary and `GUIA_ANALISIS_TESIS.md`. |
| `11_deduplicate_dengue.R` | Selects representatives from globally marked analytical-content groups; exports deduplicated annual partitions and an exclusion audit. |
| `12_generate_dengue_graphics.R` | Defines 12 dengue figures, including temporal patterns, geography, case classification, care, social characteristics and comparative event trends. |
| `13_integrate_dengue_events.R` | Compares area, insurance regime and ethnic belonging across the events. Uses event-210 summary tables and event-220/580 RDS inputs, with code-to-label fallback where needed. |
| `14_integrated_care_figure.R` | Generates a six-panel comparison of hospitalization and recorded care-related intervals from deduplicated analytical data. |

These are source-code descriptions, not a guarantee that every figure is produced for every new input release. Review execution logs and output indexes.

## Configuration

Run the entry point from the directory containing the scripts. `00_config.R` checks that `01_load_databases.R` exists in the current working directory.

For a separate data root:

```r
Sys.setenv(SIVIGILA_PATH = "/absolute/path/to/SIVIGILA_data")
setwd("/absolute/path/to/DENV/code_SIVIGILA")
```

The environment variable can also be defined in `.Renviron`. If it is absent, the default root is the parent of the script directory. Setting the variable does not move the scripts or change the independent log path in `run_all.R`.

If a local package library exists under `R_library/<R-major.minor>` in the data root, the configuration adds it to `.libPaths()`. Its package-availability check occurs before this addition, so required packages must already be discoverable when configuration begins.

## How to run

From `code_SIVIGILA/`:

```r
source("000_run.R", encoding = "UTF-8")
```

Alternatively:

```r
source("run_all.R", encoding = "UTF-8")
correr_pipeline()
```

Partial runs select a range of scripts; they do not automatically reconstruct omitted upstream outputs:

```r
source("run_all.R", encoding = "UTF-8")
correr_pipeline(hasta = "10")  # scripts 00–10

# Requires previously generated event data and the event-210 figure tables.
correr_pipeline(desde = "13") # scripts 13–14
```

`detener_en_error = TRUE` is the default. Existing output files can be overwritten by reruns. Preserve a copy of a completed analysis before rerunning it with changed inputs or dependencies.

## Outputs

Output filenames and organization differ between events. In particular, event 210 is saved as annual partitions rather than the single historical RDS files used for event 220.

| Event/output | Main location or filename |
| --- | --- |
| Dengue canonical partitions | `resultados/dengue/datos_procesados/canonica_por_anio/dengue_canonica_<year>.rds` |
| Dengue analytical partitions | `resultados/dengue/datos_procesados/analitica_por_anio/dengue_analitica_<year>.rds` |
| Deduplicated dengue partitions | `canonica_sin_duplicados_por_anio/` and `analitica_sin_duplicados_por_anio/` under the same data directory |
| Severe-dengue deduplicated tables | `resultados/dengue_grave/datos_procesados/dengue_grave_canonica_sin_duplicados.rds` and `dengue_grave_analitica_sin_duplicados.rds` |
| Mortality deduplicated table | `resultados/mortalidad_dengue/datos_procesados/mortalidad_dengue_canonica_sin_duplicados.rds` |
| Mortality plotting subset | `resultados/mortalidad_dengue/datos_procesados/mortalidad_dengue_muestra_graficas.rds` |
| Exclusion audits | `resultados/<event>/reportes/deduplicacion/registros_excluidos_duplicados.csv` |
| Event figures | `resultados/<event>/figuras/`, with topic folders and `_previsualizaciones/` |
| Figure inventories and package versions | `INDICE_FIGURAS.csv` and `VERSIONES_PAQUETES.csv` in the event figure directories |
| Comparative social figure | `resultados/integracion_eventos/figuras/01_perfil_social_comparativo.tiff` |
| Comparative care figure | `resultados/integracion_atencion/Figura_2_atencion_tres_eventos.tiff` |
| Saved event plot objects | `salidas/auditoria_integral/objetos_figuras/` |

Scripts 05, 08 and 12 export TIFF, PDF and PNG figure versions. Scripts 13 and 14 currently export TIFF figures only.

The analytical tables generally select variables for analysis; this does not imply that all rows with quality flags have been removed. The mortality plotting subset additionally restricts the recorded final condition to deceased.

## Data quality and interpretation

- **Textual missingness:** `es_vacio()` recognizes empty values, selected textual sentinels and strings beginning with `*`. It does not universally treat `"000"` as missing. Zero-only municipal codes are handled separately by `unir_divipola_notificacion()`.
- **Dates:** `parsear_fecha()` handles Excel serial dates and multiple textual formats, stripping appended time components. Ambiguous day/month strings still require examination in the context of the source files.
- **Duplicate column names:** the event-210 importer keeps the first occurrence of a repeated column name. The existence of `flag_cod_evento_duplicado_inconsistente` alone does not prove that every discarded duplicate column was compared before removal.
- **Quality flags:** these describe data consistency. They can be visualized as audit results, but should not be interpreted as clinical outcomes.
- **Deduplication:** the exclusion criterion is matching content in the fields compared by each event script, after that script's normalization and exclusions. It is not validated patient-level or episode-level linkage. Repeated `CONSECUTIVE` identifiers and strict epidemiological keys are flagged rather than used alone to exclude records in the event-210/220 workflows.
- **Representative selection:** duplicate groups prioritize the most recent `fec_aju`, then `fec_arc_xl`; `consecutive` is used as a further ordering criterion. Retain the exclusion audit alongside the analysis.
- **Across-event comparison:** the workflow does not link individuals across events 210, 220 and 580. Adding the three counts does not establish a total number of unique dengue episodes. A fatal final condition within event 220 and a mortality notification in event 580 are distinct analytical sources.
- **Denominators:** maps display counts, not incidence rates. Population denominators must be incorporated separately. The severe-dengue proportion uses event-210 plus event-220 notifications as its stated denominator; it is not population incidence.
- **Geography:** occurrence, residence and notification codes refer to different locations. Use the geographical definition appropriate to each indicator.
- **Cartographic display:** island geometries are enlarged fourfold in linear scale and their separation is adjusted for readability. Their displayed separation does not represent real distance.
- **Historical comparisons:** changes in definitions, completeness, ascertainment and coding can affect trends. The present-day dictionary does not replace documentation of historical changes.

## Loading processed data

From the script directory:

```r
source("load_processed_dengue.R", encoding = "UTF-8")

# Explicitly request the deduplicated analytical partitions.
dengue_2025 <- cargar_dengue_procesado(
  2025, tipo = "analitica", sin_duplicados = TRUE
)

dengue_2020_2025 <- cargar_dengue_procesado(
  2020:2025, tipo = "analitica", sin_duplicados = TRUE
)
```

The current default is `sin_duplicados = FALSE`; use `TRUE` for analyses intended to use the deduplicated partitions.

`combinar = FALSE` returns a named list, but **still loads every requested partition into memory**. It avoids building the combined table, not the loading of all selected years. For lower memory use, request one year at a time:

```r
for (anio in SIVIGILA_ANIOS) {
  datos_anuales <- cargar_dengue_procesado(
    anio, tipo = "analitica", sin_duplicados = TRUE,
    combinar = TRUE, verbose = FALSE
  )
  # Calculate and save the required annual summary here.
  rm(datos_anuales)
  invisible(gc())
}
```

## Reproducibility and software citations

Keep the input inventory, download dates, source commit, logs, exclusion audits, output indexes, R version and dependency versions with each completed analysis.

The references below document the software and official sources. CRAN pages describe the packages consulted on 6 October 2026. Obtain the recommended references from the actual analysis environment with `citation()` [28]:

```r
citation()
citation("ggplot2")
citation("sf")
citation("lubridate")
sessionInfo()
sf::sf_extSoftVersion()
```

The accompanying optional `cite_packages.R` exports the recommended citations for R and all 15 external packages, BibTeX records, installed versions, session information and spatial-library versions. It does not rerun the analysis or install packages. If an environment lockfile is added later, create it from the validated analysis environment rather than assuming the latest available package versions reproduce a previous run.

## How to cite this project

Suggested institutional reference:

> Cadena-Caballero, C.E., Vera-Cala, L.M., Barrios Hernández, C.J., Martinez-Perez, F. (2026). SIVIGILA Dengue Pipeline: historical dengue surveillance data processing in R. Laboratory of Applied Cellular Genomics, Industrial University of Santander, Bucaramanga, Santander, Colombia. Version 0.1. Available from: https://github.com/GenomicUIS/DENV/tree/main/code_SIVIGILA

For reproducibility, include the release tag or full commit used. The version inspected while preparing this documentation was commit `aafa9fd770564615037ba556b365719327e4322a`, dated 6 October 2026. This is a source identifier, not a software DOI or proof of runtime validation. Cite the thesis separately. Cite INS data [2] and DANE cartography [6] separately from the code.

## License

This project is licensed under the **MIT License** [27]. It permits use, modification and redistribution, including commercial use, provided that the copyright and permission notices are retained. The proposal applies to original pipeline code and associated original documentation within `code_SIVIGILA/`.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

INS data, DANE cartography, R and third-party packages retain their own rights and conditions. They are not relicensed by the pipeline's MIT license. Preserve applicable notices when distributing third-party material or dependencies. The proposal does not automatically apply to the other projects or scripts at the root of the DENV repository.

Academic users are encouraged to cite this software, its developers, the data sources and the relevant packages. Academic citation is a request for recognition; it is not an additional restriction inserted into the MIT license.

## References

1. Cadena-Caballero, C.E., Vera-Cala, L.M., Barrios Hernández, C.J., Martinez-Perez, F. (2026). SIVIGILA Dengue Pipeline: historical dengue surveillance data processing in R. Laboratory of Applied Cellular Genomics, Industrial University of Santander, Bucaramanga, Santander, Colombia. Available from: https://github.com/GenomicUIS/DENV/tree/main/code_SIVIGILA
2. Instituto Nacional de Salud. Portal Sivigila web 4.0: búsqueda de microdatos [Internet]. Bogotá: INS; [cited 2026 Oct 6]. Available from: https://portalsivigila.ins.gov.co/buscador
3. Instituto Nacional de Salud. Diccionario de datos Sivigila 2026 [Internet]. Bogotá: INS; 2026 [cited 2026 Oct 6]. Available from: https://www.ins.gov.co/BibliotecaDigital/diccionario-datos-sivigila-2026-VF.pdf
4. Instituto Nacional de Salud. Manual del usuario Sivigila 4.0 [Internet]. Bogotá: INS; [cited 2026 Oct 6]. Available from: https://www.ins.gov.co/BibliotecaDigital/manual-del-usuario-sivigila-4-0.pdf
5. Instituto Nacional de Salud. Ficha de notificación: dengue (210), dengue grave (220), mortalidad por dengue (580). Versión 04 [Internet]. Bogotá: INS; 2015 [cited 2026 Oct 6]. Available from: https://www.ins.gov.co/Direcciones/Vigilancia/sivigila/FichasdeNotificacion/DENGUE%20F210-220-580.pdf
6. Departamento Administrativo Nacional de Estadística. Marco Geoestadístico Nacional 2025: capa Departamento (ID 319) [dataset on the Internet]. Bogotá: DANE; 2025 [cited 2026 Oct 6]. Available from: https://geoportal.dane.gov.co/mparcgis/rest/services/MGN2025/Serv_CapasMGN_2025/FeatureServer/319
7. R Core Team. R: a language and environment for statistical computing [computer program on the Internet]. Vienna: R Foundation for Statistical Computing; [cited 2026 Oct 6]. Available from: https://www.r-project.org/
8. Wickham H, François R, Henry L, Müller K, Vaughan D. dplyr: a grammar of data manipulation [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=dplyr. doi:10.32614/CRAN.package.dplyr.
9. Wickham H. stringr: simple, consistent wrappers for common string operations [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=stringr. doi:10.32614/CRAN.package.stringr.
10. Hester J, Wickham H, Csárdi G. fs: cross-platform file system operations based on libuv [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=fs. doi:10.32614/CRAN.package.fs.
11. Wickham H, Henry L. purrr: functional programming tools [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=purrr. doi:10.32614/CRAN.package.purrr.
12. Wickham H, Bryan J. readxl: read Excel files [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=readxl. doi:10.32614/CRAN.package.readxl.
13. Firke S. janitor: simple tools for examining and cleaning dirty data [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=janitor. doi:10.32614/CRAN.package.janitor.
14. Wickham H, Vaughan D, Girlich M. tidyr: tidy messy data [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=tidyr. doi:10.32614/CRAN.package.tidyr.
15. Wickham H, Chang W, Henry L, Pedersen TL, Takahashi K, Wilke C, et al. ggplot2: create elegant data visualisations using the grammar of graphics [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=ggplot2. doi:10.32614/CRAN.package.ggplot2.
16. Wickham H, Pedersen TL, Seidel D. scales: scale functions for visualization [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=scales. doi:10.32614/CRAN.package.scales.
17. Pedersen TL. patchwork: the composer of plots [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=patchwork. doi:10.32614/CRAN.package.patchwork.
18. Pebesma E. sf: simple features for R [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=sf. doi:10.32614/CRAN.package.sf.
19. Pebesma E. lwgeom: bindings to selected liblwgeom functions for simple features [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=lwgeom. doi:10.32614/CRAN.package.lwgeom.
20. Spinu V, Grolemund G, Wickham H. lubridate: make dealing with dates a little easier [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=lubridate. doi:10.32614/CRAN.package.lubridate.
21. Wickham H, Hester J, Bryan J. readr: read rectangular text data [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=readr. doi:10.32614/CRAN.package.readr.
22. Müller K, Wickham H. tibble: simple data frames [computer program on the Internet]. CRAN; [cited 2026 Oct 6]. Available from: https://cran.r-project.org/package=tibble. doi:10.32614/CRAN.package.tibble.
23. Wickham H. ggplot2: elegant graphics for data analysis. 2nd ed. Cham: Springer; 2016. doi:10.1007/978-3-319-24277-4.
24. Pebesma E. Simple features for R: standardized support for spatial vector data. R J. 2018;10(1):439–46. doi:10.32614/RJ-2018-009.
25. Pebesma E, Bivand R. Spatial data science: with applications in R. Boca Raton: Chapman and Hall/CRC; 2023. doi:10.1201/9780429459016.
26. Grolemund G, Wickham H. Dates and times made easy with lubridate. J Stat Softw. 2011;40(3):1–25. doi:10.18637/jss.v040.i03.
27. Open Source Initiative. The MIT License [Internet]. OSI; [cited 2026 Oct 6]. Available from: https://opensource.org/license/mit
28. R Core Team. Citing R and R packages in publications: citation [Internet]. R documentation; [cited 2026 Oct 6]. Available from: https://stat.ethz.ch/R-manual/R-devel/library/utils/html/citation.html
