# Linfocitos T en piel de vitiligo, reanálisis de scRNA-seq (GSE298871)

## Resumen

Se reanalizaron 17 biopsias de piel (6 sanas, 5 no lesionales y 6 lesionales de pacientes con vitiligo) con secuenciación de ARN de célula única, para preguntar si los linfocitos T de la piel no lesional ya están alterados antes de que aparezca la lesión. Los linfocitos T de piel lesional mostraron una respuesta a interferón gamma, el único resultado que superó la corrección por comparaciones múltiples, aunque depende de dos muestras. En piel no lesional no se encontró un cambio robusto, y las tendencias observadas quedan como hipótesis para validar en otra cohorte.

![Respuesta a interferón por condición](results/figuras/firma_IFN_principal.png)

*Cada punto es una muestra. Los puntos blancos son 486_LS y BIO023_LS, las dos muestras con mayor respuesta.*

## Pregunta y antecedentes

En el vitiligo, linfocitos T CD8 destruyen a los melanocitos, las células que producen el pigmento de la piel. La piel lesional (LS) es la que ya perdió el pigmento, y la no lesional (NL) es piel pigmentada de los mismos pacientes.

Gellatly et al. (2021) describieron activación inmune subclínica en piel no lesional, contenida en parte por linfocitos T reguladores (Treg). Brunner et al. (2025), autores del dataset, encontraron con RNA-seq bulk de 15 pacientes que la piel no lesional ya tiene marcadores inflamatorios de los ejes T1, T2 y T17. Este proyecto cambia la estrategia: aísla los linfocitos T y los compara entre condiciones con la persona como unidad de análisis, para ver si esa alteración del tejido también aparece dentro de los linfocitos T.

## Datos

GSE298871 en GEO, generado con la química 5' de 10x Genomics. GEO solo ofrece las matrices de conteo de Cell Ranger 7.2, sin lecturas originales ni datos clínicos por paciente. El sexo se dedujo de genes del cromosoma Y.

| Condición | Muestras | Muestras con 20 o más linfocitos T |
|---|---|---|
| Sano | 6 | 5 |
| NL | 5 | 4 |
| LS | 6 | 5 |

## Métodos

Cada paso está en un notebook, con la explicación de los conceptos en sus celdas de texto.

| Paso | Herramienta | Decisión principal | Motivo |
|---|---|---|---|
| Control de calidad (`02`) | Scanpy, Scrublet | Umbrales por muestra a 5 MAD en genes y moléculas, y a 3 MAD o 20% en ARN mitocondrial | La mediana mitocondrial varió entre 1.9% y 8.9% según la muestra, y un umbral fijo habría castigado más a unas muestras que a otras |
| Integración (`03`) | Harmony | Corregir por muestra en la PCA, sin tocar los conteos | Las muestras tenían diferencias técnicas, y los conteos crudos se necesitan para la expresión diferencial |
| Anotación (`03`) | Leiden, Wilcoxon | Marcadores conocidos, y genes más distintivos para los grupos sin marcador claro | Así se identificaron Schwann y glándula sudorípara sin suponerlas de antemano |
| Subtipos de T (`04`) | Leiden, Wilcoxon | Excluir los segmentos V y J del receptor de T al agrupar | Evita que una clona expandida forme su propio grupo |
| Unidad de análisis (`04` a `09`) | | Resumir los linfocitos T de cada biopsia en un solo valor antes de comparar | Dos muestras contienen el 83% de los T lesionales, y comparar célula por célula produce falsos positivos |
| Composición (`04`) | Mann-Whitney, Spearman | Porcentaje de cada subtipo por muestra | Pruebas por rangos, adecuadas para 4 o 5 valores por grupo |
| Expresión diferencial (`05`, `06`) | DESeq2, apeglm | Pseudobulk: suma de conteos de los T de cada muestra, 3,521 genes evaluados | Es el método recomendado para comparar condiciones con datos de célula única |
| Firmas (`07`) | `score_genes` | Seis programas: residencia (TRM), Th1, Th2, Th17, citotoxicidad e interferón | Th2 y Th17 porque Brunner et al. los reportaron en NL, el resto por su papel conocido en vitiligo |
| Quimiocinas (`08`) | Scanpy | Porcentaje de células con CXCL9, CXCL10 o CXCL11 por tipo y muestra | LIANA se descartó porque su filtro eliminaba las tres quimiocinas |
| Sensibilidad (`07`, `09`) | DESeq2 | Repetir sin muestras atípicas, con sexo en el modelo y LS contra NL del mismo paciente | Comprobar que los resultados no dependen de una muestra o de un desbalance |

Todos los valores p se corrigieron con Benjamini-Hochberg. El valor corregido se reporta como q, y el umbral de significancia es q < 0.05.

Como control, los linfocitos T por muestra coinciden de cerca con los que reportaron Brunner et al. en su tabla suplementaria S6 (por ejemplo 458 contra 433 en DER115_LS). Los melanocitos estuvieron prácticamente ausentes en las seis muestras lesionales (0 a 0.1%), lo esperado en vitiligo.

![Tipos celulares](results/figuras/umap_tipos_celulares.png)

## Resultados

### Subtipos de linfocitos T

Se obtuvieron 1,445 linfocitos T en seis subtipos. Al revisar los genes más distintivos de cada uno, el grupo nombrado al principio "CD4 residente" resultó estar definido por genes de estrés por disociación (FOS, JUNB, HSPA1B), que se activan cuando la biopsia se digiere con enzimas (van den Brink et al., 2017). Se renombró "CD4 con estrés". Los demás subtipos se confirmaron, por ejemplo Treg con FOXP3, IKZF2 y TIGIT.

![Marcadores de los subtipos](results/figuras/dotplot_T_subtipos.png)

### La piel lesional muestra respuesta a interferón gamma

La firma de respuesta a interferón fue el único resultado significativo (tendencia de Spearman, rho = 0.69, q = 0.037). No es un gradiente: NL se parece a Sano y LS está más alto. Las otras cinco firmas no cambiaron ([figura](results/figuras/firmas_T_por_condicion.png)).

| Escenario | p LS contra Sano | p tendencia |
|---|---|---|
| Todas las muestras | 0.032 | 0.006 |
| Sin la muestra sana atípica CON194 | 0.063 | 0.016 |
| Sin DER114_LS y DER115_LS | 0.036 | 0.006 |
| Sin 486_LS y BIO023_LS | 0.143 | 0.055 |

El resultado no depende de las muestras con más linfocitos T, pero sí de 486_LS y BIO023_LS, que podrían ser lesiones activas. El heatmap muestra que los genes que suben son STAT1, IRF1, GBP2 y GBP5, que responden sobre todo al interferón gamma, mientras que ISG15 e IFI6, más propios del interferón tipo I, no siguen la condición. Esto coincide con el mecanismo del vitiligo: los CD8 producen interferón gamma, que activa la vía JAK y STAT1 e induce CXCL9 y CXCL10, que atraen más linfocitos T (Rashighi et al., 2014). Junto con la ausencia de melanocitos, funciona como control positivo del análisis.

![Genes de interferón y residencia por muestra](results/figuras/heatmap_T_IFN_TRM.png)

### En piel no lesional no se detectó un cambio robusto

La composición mostró tendencias hacia más CD8 GZMK+ y menos Treg en NL (p = 0.063 a 0.11), sin superar la corrección (q mínimo de 0.19). La tendencia hacia menos CD4 con estrés probablemente es técnica.

![Subtipos por condición](results/figuras/subtipos_T_por_condicion.png)

En la expresión diferencial, solo FXYD2 fue significativo en NL contra Sano (q = 0.029), seguido de ZNF683 (q = 0.20), que codifica Hobit, un factor asociado a residencia en tejido (Mackay et al., 2016) y a CD8 efectores en humanos (Vieira Braga et al., 2015). Ninguno de los dos es un resultado confiable.

| Motivo | Detalle |
|---|---|
| Pocos conteos | Alrededor de 10 moléculas por muestra |
| Sensibles al modelo | FXYD2 sobrevive al quitar muestras (q de 0.022 a 0.024), pero pierde significancia con el sexo en el modelo (q = 0.11) |
| No son independientes | Son los dos marcadores principales del subtipo CD8 / gamma delta |
| Sin respaldo | La firma TRM completa no cambió (p = 0.73), y en el heatmap ZNF683 sube sin CD69, ITGAE ni RGS1 |
| Posible confusión | Las cuatro muestras NL vienen de donantes con prefijo DER y ninguna sana, y GEO no explica qué significa |

![Genes candidatos](results/figuras/genes_T_NL_vs_Sano.png)

La comparación de LS contra NL del mismo paciente (3 pares) no dio genes significativos, y los genes con menor p fueron de estrés por disociación. Las quimiocinas CXCL9/10/11 se detectaron en menos del 2% de las células, casi nunca en queratinocitos, aunque Brunner et al. sí detectaron CXCL9 con RNA-seq bulk. Esto apunta a un límite de detección de la técnica y no permite evaluar ese eje ([figura](results/figuras/quimiocinas_por_condicion.png)).

## Conclusión

Los linfocitos T de piel lesional muestran una respuesta a interferón gamma consistente con la biología conocida del vitiligo, sostenida principalmente por dos muestras. En piel no lesional no se encontró una alteración robusta, lo que con cuatro muestras no descarta que exista. Las tendencias de más CD8 GZMK+ y menos Treg en NL serían compatibles con la activación subclínica que describieron Gellatly et al. (2021), y son la hipótesis a validar en una cohorte independiente.

## Limitaciones

| Limitación | Efecto | Posible solución |
|---|---|---|
| Pocas muestras (5, 4 y 5) | Solo se detectan diferencias grandes. Con 4 contra 5 muestras, el p mínimo en Mann-Whitney es 0.016 | Validar en otra cohorte o combinar datasets |
| NL y LS del mismo paciente | No son independientes, y la comparación pareada solo tiene 3 pares | Más pacientes con ambas biopsias |
| Sexo desbalanceado | FXYD2 pierde significancia al ajustarlo | Cohorte con sexo balanceado |
| Origen de las muestras | Todas las NL tienen prefijo DER y ninguna sana | Consultar a los autores |
| Sin datos clínicos | No se sabe si 486 y BIO023 tenían enfermedad activa (VIDA) | Pedir los datos a los autores |
| Estrés por disociación | Define un subtipo completo y domina la comparación LS contra NL | Modelar una firma de estrés |
| Detección de quimiocinas | No se pudo evaluar el eje CXCL9/10 y CXCR3 | RNA-seq bulk, PCR o transcriptómica espacial |

## Reproducir el análisis

Se necesitan conda o mamba. Hay dos entornos, uno para Python y otro para R.

```bash
mamba env create -f envs/environment.yml
mamba env create -f envs/environment-r.yml
```

Para usar R en JupyterLab hay que registrar su kernel una vez, con los dos entornos activos:

```bash
mamba activate vitiligo-t
mamba activate --stack vitiligo-t-r
R -e 'IRkernel::installspec(name = "vitiligo-t-r", displayname = "R (vitiligo-t-r)")'
```

Después se descargan los datos desde la raíz del repositorio y se corren los notebooks en orden:

```bash
bash scripts/01_descargar.sh
```

| Notebook | Kernel | Contenido |
|---|---|---|
| `01_carga` | Python | Carga de las 17 muestras |
| `02_qc` | Python | Control de calidad y dobletes |
| `03_tipos_celulares` | Python | Integración, agrupamiento y anotación |
| `04_linfocitos_T` | Python | Subtipos de T y composición |
| `05_pseudobulk` | Python | Suma de conteos por muestra |
| `06_deseq2` | R | Expresión diferencial y heatmap |
| `07_firmas` | Python | Puntajes de firmas y sensibilidad |
| `08_quimiocinas` | Python | Expresión de CXCL9/10/11 |
| `09_sensibilidad` | R | Expresión diferencial sin muestras atípicas y LS contra NL pareado |

Las figuras se guardan en `results/figuras` y las tablas en `results/tablas`.

## Referencias

**Biología**

Brunner PM, David E, Del Duca E, et al. Transcriptomic profiling of vitiligo patients shows polar immune dysregulation in involved and uninvolved skin. *J Allergy Clin Immunol*. 2025, 156(4):993-1007. [doi:10.1016/j.jaci.2025.06.002](https://doi.org/10.1016/j.jaci.2025.06.002). Datos: [GSE298871](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298871).

Gellatly KJ, Strassner JP, Essien K, et al. scRNA-seq of human vitiligo reveals complex networks of subclinical immune activation and a role for CCR5 in Treg function. *Sci Transl Med*. 2021, 13(610):eabd8995. [doi:10.1126/scitranslmed.abd8995](https://doi.org/10.1126/scitranslmed.abd8995)

Rashighi M, Agarwal P, Richmond JM, et al. CXCL10 is critical for the progression and maintenance of depigmentation in a mouse model of vitiligo. *Sci Transl Med*. 2014, 6(223):223ra23. [doi:10.1126/scitranslmed.3007811](https://doi.org/10.1126/scitranslmed.3007811)

Shiu J, Zhang L, Lentsch G, et al. Multimodal analyses of vitiligo skin identify tissue characteristics of stable disease. *JCI Insight*. 2022, 7(13):e154585. [doi:10.1172/jci.insight.154585](https://doi.org/10.1172/jci.insight.154585)

Mackay LK, Minnich M, Kragten NAM, et al. Hobit and Blimp1 instruct a universal transcriptional program of tissue residency in lymphocytes. *Science*. 2016, 352(6284):459-463. [doi:10.1126/science.aad2035](https://doi.org/10.1126/science.aad2035)

Vieira Braga FA, Hertoghs KML, Kragten NAM, et al. Blimp-1 homolog Hobit identifies effector-type lymphocytes in humans. *Eur J Immunol*. 2015, 45(10):2945-2958. [doi:10.1002/eji.201545650](https://doi.org/10.1002/eji.201545650)

van den Brink SC, Sage F, Vértesy Á, et al. Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations. *Nat Methods*. 2017, 14(10):935-936. [doi:10.1038/nmeth.4437](https://doi.org/10.1038/nmeth.4437)

**Software y métodos**

Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018, 19:15. [doi:10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)

Wolock SL, Lopez R, Klein AM. Scrublet: computational identification of cell doublets in single-cell transcriptomic data. *Cell Syst*. 2019, 8(4):281-291. [doi:10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005)

Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods*. 2019, 16:1289-1296. [doi:10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0)

Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. *Sci Rep*. 2019, 9:5233. [doi:10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)

McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for dimension reduction. arXiv. 2018. [arXiv:1802.03426](https://arxiv.org/abs/1802.03426)

Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol*. 2014, 15:550. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics*. 2019, 35(12):2084-2092. [doi:10.1093/bioinformatics/bty895](https://doi.org/10.1093/bioinformatics/bty895)

Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B*. 1995, 57(1):289-300. [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)
