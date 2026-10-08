# PAAD — revisión de los datasets descargados

Recomendación de Claude tras verificar los registros de muestra de cada serie.
**No borra nada**: la decisión es tuya. Fecha: 2026-10-06.

Veredicto del pipeline: `1021 encontrados -> 251 pasaron puertas 1-2 -> 17 confirmados + 13 ambiguos`

Es el cáncer más productivo con diferencia. **18 datasets utilizables, ~861 muestras**,
frente a las 183 de TCGA-PAAD.

---

## El hallazgo principal

### GSE93326 — 204 muestras, epitelio y estroma separados por microdisección

```
tissue: pdac
compartment: epithelium
patient id: cumc_001
library: nugen
```

*"Gene expression from laser capture microdissected pancreatic cancer epithelium and
stroma"*. Usa **microdisección por captura láser** para separar los dos compartimentos.

Esto importa de forma directa para el proyecto: **CLDN1 es una proteína de unión estrecha
epitelial**. En bulk RNA-seq convencional su señal queda diluida por el estroma, que en
PDAC es masivo (puede ser el 80% del tumor). Este dataset te permite medir CLDN1 en el
compartimento donde realmente se expresa.

Es el dataset más valioso de todo el esfuerzo en GEO, en los cuatro cánceres procesados.

---

## Conservar — 18 datasets

### Confirmados por el pipeline (15)

| Dataset | Muestras | Qué es |
|---|---|---|
| **GSE93326** | 204 | LCM epitelio/estroma (ver arriba) |
| **GSE124230** | 97 | QuantSeq |
| **GSE254877** | 84 | Cohorte |
| **GSE162689** | 59 | **Páncreas normal**: islotes + tejido exocrino, donantes no diabéticos |
| **GSE225767** | 55 | Biorepositorio |
| **GSE179351** | 54 | Cohorte |
| **GSE124231** | 48 | Cohorte |
| **GSE151580** | 33 | Cohorte |
| **GSE306366** | 24 | Cohorte |
| **GSE169321** | 20 | Cohorte |
| **GSE135686** | 18 | Cohorte |
| **GSE336518** | 17 | PDAC |
| **GSE161208** | 10 | PDAC |
| **GSE214894** | 8 | Cohorte |
| **GSE235244** | 8 | hPDAC |

**Nota sobre GSE162689:** declara `tissue: islet tissue` y `tissue: exocrine tissue
adjacent to islets`, con `disease state: non-diabetic`. Es **páncreas sano**, no cáncer.
Solapa con la función de GTEx (362 muestras de páncreas), así que su valor marginal es
limitado — pero con TCGA-PAAD teniendo solo **4 normales**, cualquier tejido sano extra
cuenta. Decídelo tú.

### Ambiguos que sí valen (3)

| Dataset | Muestras | Lo que declara | Por qué conservarlo |
|---|---|---|---|
| **GSE242915** | 79 | `tissue: tumor`, `tumor-draining lymph nodes`, `chemo-naive` | Tejido de paciente sin tratar. Ojo: mezcla tumor y ganglios |
| **GSE319924** | 36 | `tissue: stroma`, `id paired: s001`, `timepoint: pre` | Estroma pareado pre/post quimio |
| **GSE226307** | 7 | `tissue type: poorly/moderately differentiated pdac` | PDAC genuino con grado histológico |

---

## Borrar — 12 datasets

```bash
cd ~/master-project/dataset/GEO/PAAD
rm -rf GSE199102 GSE226829 GSE309679 GSE185190 GSE160306 GSE272019 \
       GSE330683 GSE64018 GSE295683 GSE265950 GSE313827 GSE63124
```

### Transcriptómica espacial, no bulk RNA-seq (2)

**GSE199102** — 608 "muestras", 161 MB. Confirmado por el pipeline como tejido, y lo es
(`tissue: pancreas`), pero el título es *"**Spatial** transcriptome profiling of
pancreatic cancer"*. Sus archivos lo delatan: `BioProbeCountMatrix`, `.pkc`,
`segments_annotation` — es **GeoMx de NanoString**. Las 608 entradas son regiones de
interés de unos pocos cortes, no muestras independientes.

**GSE226829** — 277 "muestras", **768 MB**. El mismo caso, y además descargó cuatro
imágenes PNG de 141-212 MB cada una (711 MB de los 768). Es el único dataset que ha
gastado disco de forma seria en todo el proceso.

> **Limitación detectada.** El pipeline no distingue GeoMx/espacial de bulk RNA-seq. El
> tipo de dataset de GEO es el mismo (*expression profiling by high throughput
> sequencing*) y las muestras declaran tejido legítimamente. Los indicios están en los
> nombres de archivo (`BioProbeCountMatrix`, `.pkc`, `segments_annotation`) y en el
> título. Es material para la sección de limitaciones.

### Especie mixta (1)

**GSE309679** — 13 muestras. Declara `tissue: human normal pancreas` **y**
`murin bat: brown adipose tissue` (tejido adiposo marrón de ratón). Sus archivos son
`hsa-counts.txt.gz` y `mmu-counts.txt.gz`. Es un **paper metodológico** sobre inferencia
de redes ligando-receptor que usa varios tejidos y especies.

### No es tejido pancreático (4)

| Dataset | Muestras | Lo que declara |
|---|---|---|
| **GSE185190** | 80 | `peripheral blood`, `treatment: placebo` — ensayo de diabetes tipo 1 |
| **GSE272019** | 48 | `cell type: cd14+ monocytes` — monocitos en cultivo |
| **GSE160306** | 79 | `tissue: retina macula` — retinopatía diabética. **Ya apareció en CHOL** |
| **GSE64018** | 24 | `countlevel_12asd_12ctl` — 12 autismo vs 12 control, tejido cerebral |

### Cultivo (2)

**GSE295683** — 66 muestras, archivos `KRAS_Mutant_HPNE_Raw_counts.xlsx`. **HPNE** es una
línea celular pancreática inmortalizada. Quedó ambiguo solo porque el `.xlsx` no se puede
verificar sin descargarlo entero.

**GSE330683** — 18 muestras, `cell type: cafs` (fibroblastos asociados a cáncer, que se
cultivan). Además su matriz tiene **478 columnas numéricas frente a 18 muestras
declaradas** — la comprobación de consistencia saltó correctamente.

### Sin verificar por mí (3)

Estos tres los dejo en "borrar" por prudencia, pero **no los he comprobado**. Si quieres
recuperarlos, mira su `summary.md` y su página de GEO:

| Dataset | Muestras | Motivo de la marca |
|---|---|---|
| GSE265950 | 12 | El resumen menciona línea celular, organoide y no humano |
| GSE313827 | 18 | El resumen no nombra cáncer; identificadores no reconocidos |
| GSE63124 | 20 | La matriz tiene 8 columnas frente a 20 muestras declaradas |

---

## Precisión de esta ejecución

| | Resultado |
|---|---|
| Confirmados correctos | **15 de 17** (88%) |
| Fallos | 2, los dos por transcriptómica espacial |
| Ambiguos que sí valen | 3 de 13 |

Los dos fallos son del mismo tipo nuevo: **GeoMx espacial que declara tejido
correctamente**. No es un error de la puerta 4 — las muestras *son* tejido. Es que el
pipeline no comprueba la **modalidad** del experimento.

---

## Datos para la memoria

De los 1,021 candidatos de páncreas, los motivos de descarte dominantes:

```
   361  solo un _RAW.tar por muestra, sin matriz de counts     (35%)
   120  ningún archivo reconocible como matriz de counts       (12%)
   107  solo valores normalizados (FPKM/TPM/CPM)               (10%)
    51  los counts llamados "raw" contienen decimales           (5%)
    39  marca de tipo de muestra en la cabecera de la matriz
    32  sin archivos suplementarios: los datos solo están en SRA
    16  solo counts de circRNA/isoformas/miRNA
```

**El 62% se descarta por motivos puramente técnicos de publicación**, no por el contenido
biológico. Es un argumento concreto y cuantificado sobre las limitaciones de reutilizar
datos públicos.

Y entre los que llegaron a la puerta 4, los nombres de línea celular más repetidos fueron
**PANC-1** (16 series) y **MIA PaCa-2** (7). La investigación en PDAC en GEO está
dominada por esas dos líneas.
