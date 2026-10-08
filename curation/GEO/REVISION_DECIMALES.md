# Revisión de los rechazos por decimales

Fecha: 2026-10-07. Motivo: comprobar si la regla `decimales ⇒ no son counts`
descartó datasets recuperables.

**Resultado: 36 datasets con 2 527 muestras de tejido humano fueron rechazados
por error.** Es más de lo que aportan CHOL, ESCA, PAAD y STAD juntos (2 198).

---

## 1. Por qué se revisó

La regla está implementada como **hecho**, no como sospecha: rechaza sin
apelación y el dataset no se descarga.

```python
if whole is False:
    hard.append("the 'raw' counts contain decimals, so they are not counts")
```

Es correcta para FPKM, TPM y CPM, donde la profundidad de secuenciación ya se
eliminó y la información es irrecuperable. Pero **no distingue esos de los
*expected counts* de RSEM y salmon**, que son decimales legítimos.

RSEM reparte las lecturas ambiguas —las que encajan en varias isoformas de un
mismo gen— de forma proporcional a su probabilidad. Si una lectura tiene un 70 %
de probabilidad de venir del transcrito A, suma 0,7 a ese gen. El resultado
conserva la escala del conteo y la relación con la profundidad, que es lo que
DESeq2 necesita. Redondearlos es el procedimiento estándar: lo hace `tximport`,
el paquete de Bioconductor para esto.

## 2. Método

De los seis `decisions.csv`, los datasets cuyo **único** motivo de rechazo era
el decimal: **54 datasets, 3 481 muestras**. (Hubo 237 con ese motivo, pero 183
tenían además otro motivo duro, así que seguirían rechazados.)

Para cada uno se pidió la cabecera y las primeras 60 filas del archivo de counts
por rango HTTP, reutilizando `peek_file()` del propio pipeline, y se midió la
distribución de valores.

**El discriminador es la magnitud, no la cantidad de decimales:**

| | Mediana de genes expresados | Máximo |
|---|---|---|
| counts (enteros o de RSEM) | decenas a miles | decenas de miles o más |
| log2, VST, FPKM, TPM | < 30 | < 300 |

Un log2 o un VST **comprimen** el rango: todo cabe por debajo de ~20. Un count
no, porque conserva la escala absoluta. Los 8 casos que quedaron en el límite se
abrieron a mano.

## 3. Resultado

| Veredicto | Datasets | Muestras |
|---|---|---|
| **Recuperables, con tejido declarado** | **36** | **2 527** |
| Recuperables, muestras sin declarar | 2 | 111 |
| No recuperables (log2 / VST / TPM) | 15 | 486 |
| Rechazar por otro motivo | 1 | 357 |

Detalle por dataset en `curation/GEO/decimal_rejections.csv`.

Los 36 recuperables **ya pasaron la puerta 4**: los motivos se acumulan, así que
la puerta 4 corrió igualmente y confirmó `declares: tissue` en todos. No son
líneas celulares: son biopsias humanas con counts utilizables.

### Los no recuperables son inequívocos

| Dataset | n | Mediana | Máximo | Qué es |
|---|---|---|---|---|
| GSE179443 | 137 | 0.92 | 12.1 | log2 |
| GSE133979 | 68 | 6.23 | 13.4 | el nombre lo dice: `variance_stabilised` |
| GSE210351 | 41 | 7.77 | 15.2 | el nombre lo dice: `vsdCounts` |
| GSE317207 | 33 | 6.53 | 13.2 | log2 |
| GSE316796 | 2 | 0.50 | 7.2 | `CorrectedCountsMatrix` — transformado |

Todo el rango cabe entre 0 y 20. No hay nada que redondear.

### GSE140182 sigue rechazado, por otro motivo

357 muestras, y su archivo se llama `M1M2_processed_expected_count_matrix`:
counts de RSEM, enteros de hecho. Pero las columnas son `Donor1.M1_01`,
`Donor1.M1_07`… — **macrófagos M1/M2 sorteados por donante**, no tejido. El
título lo confirma: *"Single cell RNA sequencing of advanced gastric cancer"*.
El rechazo era correcto; el motivo registrado, no.

## 4. Dos modos de fallo distintos, no uno

### 4.1. RSEM y salmon (34 de los 36)

El caso esperado. Muchos lo dicen en el nombre del archivo y aun así se
rechazaron:

```
GSE314812_rsem.merged.gene_counts.tsv.gz
GSE255122_salmon.merged.gene_counts.tsv.gz
GSE288907_salmon.merged.gene_counts.tsv.gz
GSE202853_allSample_RSEM_genes_ExpReadCount.txt.gz
GSE162960_RSEM_counts.csv.gz
GSE243188_RSEM_counts.tsv.gz
GSE198697_rsem_gene_counts.tsv.gz
GSE233421_3Samples_RSEM_genes_ExpReadCount.txt.gz
```

**Ocho nombres de archivo con `rsem`, `salmon` o `ExpReadCount` dentro.** Una
sola comprobación del nombre antes de aplicar la regla habría evitado estos
ocho, y es la mejora más obvia que queda.

### 4.2. Columnas de counts y de FPKM intercaladas (2 casos)

Este no lo esperaba. GSE146889 (176 muestras) publica las dos cosas en el mismo
archivo, columna a columna:

```
GeneId  GeneName  MSI_MLH1G_normal_10_count  MSI_MLH1G_normal_10_rpkm  ...
ENSG00000223972  DDX11L1  12  0.05619978143906497  0  0.0  5  0.03586146...
```

Los counts son **enteros perfectos**. Los decimales vienen de las columnas
`_rpkm` de al lado. El detector tokeniza la fila y encuentra decimales, así que
rechaza — sin ver que la mitad de las columnas son counts íntegros.

GSE92945 (21 muestras) hace lo mismo con `N1, N1_fpkm, N2, N2_fpkm, …`.

**Son 197 muestras rechazadas por una regla que miraba las columnas
equivocadas.** Al unificar hay que seleccionar solo las columnas `_count`.

## 5. El hallazgo que lo justifica todo

**GSE143584 — 401 muestras.**

```
Gene expression from laser capture microdissected epithelium
```

Microdisección láser de **epitelio** en adenocarcinoma ductal de páncreas, con
counts de precisión completa:

```
GENE   S1    S2    ...
A1CF   4548.694301644136   190.4097599778177   ...
```

Supera a GSE93326 (204 muestras), que hasta ahora era el hallazgo más valioso de
toda la curación, y por el mismo motivo: **CLDN1 es una proteína epitelial de
unión estrecha, y en PDAC el estroma puede ser el 80 % del tumor.** En bulk
convencional la señal se diluye; aquí se mide donde se expresa.

Con GSE119968 (20, también microdisección en colorrectal) y GSE93326 (204), el
dataset tendría **625 muestras con epitelio separado del estroma**. Eso deja de
ser un hallazgo suelto y pasa a ser un eje de análisis propio.

Otros recuperables grandes:

| Dataset | n | Cáncer | Qué aporta |
|---|---|---|---|
| GSE205154 | 289 | PAAD | primarios **y metastásicos** |
| GSE276114 | 177 | LIHC | proteo-transcriptómica de HCC avanzado |
| GSE319878 | 171 | COAD | privación socioeconómica y expresión |
| GSE310957 | 129 | PAAD | estados de fibroblastos (DeCAF) |
| GSE183984 | 113 | COAD | **longitudinal clínico-genómico** |
| GSE330465 | 104 | COAD | antes/después de FOLFOX neoadyuvante |
| GSE254660 | 68 | ESCA | con recidiva, edad y supervivencia en la cabecera |

## 6. Cómo cambia el balance

| | Antes | Después | Cambio |
|---|---|---|---|
| Datasets conservados | 112 | 148 | +36 |
| Muestras | ~5 623 | ~8 150 | **+2 527 (+45 %)** |

Por cáncer, el cambio no es uniforme:

| Cáncer | Muestras antes | Recuperadas | Después |
|---|---|---|---|
| PAAD | ~861 | **+1 054** | ~1 915 |
| COAD | ~1 543 | +826 | ~2 369 |
| LIHC | ~1 882 | +454 | ~2 336 |
| ESCA | 152 | **+96 (+63 %)** | 248 |
| STAD | 853 | +67 | 920 |
| CHOL | 332 | 0 | 332 |

**PAAD más que dobla**, y pasa a ser el segundo cáncer mejor cubierto. ESCA, que
era el más pobre, crece un 63 %.

## 7. Qué hacer con esto

Estos 36 datasets **no están descargados**: se rechazaron antes de la fase 4.
Para incorporarlos hay dos caminos.

**Opción A — descargar solo estos 36.** Un script corto que lea
`decimal_rejections.csv`, filtre `veredicto == RECUPERABLE` y `declara ==
tissue`, y baje esos archivos. No toca el pipeline ni el caché. Es lo que
recomiendo: son ~2 h de descarga y el resultado es verificable dataset a dataset.

**Opción B — descongelar el código.** Cambiar la regla a: si el nombre del
archivo contiene `rsem`, `salmon` o `expectedcount`, los decimales son sospecha
en vez de hecho; y al tokenizar, ignorar las columnas cuyo nombre acabe en
`_fpkm`, `_tpm`, `_rpkm` o `_cpm`. Después reinspeccionar los seis cánceres
(~7 h) con el caché borrado.

La opción A da el mismo resultado en un tercio del tiempo y sin riesgo de
cambiar decisiones ya revisadas. La B es la correcta si el pipeline va a volver
a usarse para otros cánceres.

## 8. Para la memoria

La cifra honesta cambia, y conviene decirla con el matiz:

> El filtro automático descartó 36 datasets con 2 527 muestras de tejido humano
> por una regla demasiado estricta: trataba cualquier valor decimal como señal de
> normalización, sin distinguir los *expected counts* de RSEM y salmon, que son
> decimales legítimos y redondeables. Se detectó al revisar manualmente los
> rechazos, y se recuperaron. El episodio ilustra el límite de cualquier filtro
> automático sobre datos de terceros: **la regla era correcta en el 85 % de los
> casos y aun así perdía el 31 % de las muestras utilizables.**

De los 54 rechazos por decimales, 15 eran correctos y 1 lo era por otro motivo.
La regla acertaba en 16 de 54 **datasets**, pero el coste no se mide en datasets
sino en muestras: las 2 527 perdidas frente a las 486 correctamente excluidas.
