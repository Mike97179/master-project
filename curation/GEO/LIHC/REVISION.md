# LIHC — revisión manual de los datasets descargados

Fecha de la ejecución: 2026-10-07 12:01 → 14:26 (2 h 25 m)
Embudo: **2036 candidatos → 363 pasan puertas 1-2 → 27 confirmados + 24 ambiguos = 51 descargados (867 MB)**

Verificación hecha a mano sobre `sample_metadata.tsv`, sobre la cabecera real de
cada matriz y, en los `.xlsx`, abriendo el archivo.

**Resultado: conservar 34 datasets (~1882 muestras), borrar 17 (661 MB de 867).**

LIHC es el cáncer con mejor cosecha clínica: cuatro cohortes traen supervivencia,
METAVIR o respuesta a inmunoterapia.

---

## A. BORRAR — 17 datasets, 661 MB

### A.1 Single-cell y single-nuclei disfrazados de bulk (3) — 597 MB, el 69 % del disco

| GSE | n declarado | Tamaño | Qué es |
|---|---|---|---|
| GSE98638 | 6 | **340 MB** | células T infiltrantes, Smart-seq2, **5063 columnas**. Trae un `HCC.bulk.S5.count` de 5 muestras, pero no justifica 340 MB |
| GSE292298 | 6 | **218 MB** | snRNA-seq; las columnas son códigos de barras `AAACCTGAGACGCACA.1_1` |
| GSE278324 | 6 | 39 MB | snRNA-seq; columnas `AAACCCAGTCTCGGAC-1_h799` |

Los tres los marcó la puerta 3 por desajuste de columnas o por identificadores
raros, y los tres son exactamente lo que sospechaba. Entre ellos, el 69 % del
disco de LIHC.

### A.2 Las matrices `.xlsx` que no resistieron la apertura (4)

Abrir los `.xlsx` a mano es lo que más información nueva dio en LIHC:

| GSE | n | Qué apareció al abrirlo |
|---|---|---|
| GSE237697 | 15 | la primera columna **no es un gen**, es una muestra (`ATP1B1-PRKACA-1`), con decimales (`5450.066`); y son **hepatocitos primarios transducidos** con el gen de fusión, no tejido |
| GSE144494 | 135 | las muestras son **CTC de pacientes `Brx71`, `Brx50`** — líneas de células tumorales circulantes de **cáncer de mama**, no hígado |
| GSE120021 | 8 | `cell type: Hep3B cells` / `SK-Hep1 cells`, xenoinjertos de líneas celulares |
| GSE82177 | 27 | la matriz **no usa símbolos de gen estándar**: las filas son `14q(0)`, `14q(I-1)`, `ACA58`, `ACA5b`. **Cero coincidencias con `CLDN*`.** Inutilizable |

GSE82177 es el segundo caso, después de GSE230524 en COAD, en que la sospecha
blanda *"the first column does not look like gene identifiers"* acierta de lleno.

### A.3 Sangre y células inmunes aisladas (5)

CLDN1 es epitelial; en estas muestras no hay epitelio.

| GSE | n | Qué es |
|---|---|---|
| GSE317443 | 20 | `tissue: blood`, HCC vs sanos |
| GSE120663 | 7 | **PBMC** de sangre periférica (estaba entre los *confirmados*) |
| GSE245905 | 9 | células innatas **sorteadas** de hígado (también confirmado) |
| GSE276062 | 4 | **monocitos** de sangre tratados con sobrenadante tumoral |
| GSE285415 | 6 | los mismos monocitos, otro experimento del mismo grupo |
| GSE289187 | 6 | células **HUVEC** (endotelio de cordón umbilical) tratadas con EVs |

### A.4 Otros (3)

| GSE | n | Motivo |
|---|---|---|
| GSE254461 | 57 | **gliomas** — el mismo falso positivo que ya salió en COAD |
| GSE92432 | 10 | método de aislamiento de células de **colon** — también repetido de COAD |
| GSE221589 | 27 | *"Human hepatocytes can give rise to intrahepatic cholangiocarcinoma"*: los títulos son `phHep-derived iCCA_12`, **iCCA generado en laboratorio** a partir de hepatocitos primarios. Solo 2 de 27 son `Patient-derived` |
| GSE255163 | 30 | **duplicado**: ya se conserva en COAD (metástasis hepática de colorrectal). Se queda allí |

---

## B. CONSERVAR — 34 datasets, ~1882 muestras

### B.1 Las cohortes con datos clínicos — lo mejor de LIHC

**GSE263786 (243) — la mayor de todas.** `tissue: liver; disease state: Normal / Disease`,
HTSeq con Ensembl versionado. Progresión de enfermedad hepática.

**GSE141198 (148) + GSE141200 (58) = 206 muestras, HCC de Taiwán, con supervivencia:**

```
etiology: HBV; ctnnb1 status: mut; efs days: 2922; efs event: 1; os days: 4206; os event: 1
```

Etiología, mutación de CTNNB1, supervivencia libre de evento y global, con
censura. Permite análisis de supervivencia, no solo expresión diferencial.

**GSE144269 (140) — HCC de Mongolia, 70 pares tumor / no-tumor** del mismo paciente.

**GSE244591 (136)** — tejido tumoral hepático. Ojo: declara
`cell line origin: Human Hepatocarcinoma`, que es un **uso incorrecto del campo** —
los títulos son `HCC_001`…`HCC_136` y es tejido. Estuvo a punto de caer por esto.

**GSE237330 (56) + GSE237331 (50) = 106 muestras de biopsias precancerosas**
con `virus: HBV/HCV`, `metavir: 0_1 / 3_4`, `hcc: yes/no` y `censored time`.
Cirrosis que progresa (o no) a HCC: el grupo de riesgo intermedio que falta en TCGA.

**GSE148355 (128)** — serie de progresión completa en un solo dataset:

```
Normal tissue, liver | Fibrosis-low | Fibrosis-high | Cirrhosis tissue | (+3 más)
```

Con edad, sexo y hepatectomía por paciente. Ensembl + símbolo, enteros, CLDN1 presente.
**Nota sobre su marca:** el motivo que arrastraba era *"the text points at LIHC, not CHOL"* —
una sospecha **heredada del caché de CHOL**, donde sí era correcta. En LIHC apunta
al cáncer correcto. Ver la sección D.

**GSE214846 (130)** — `tissue: hepatocellular carcinoma` / `paracancerous normal tissues`,
símbolos de gen, enteros, CLDN1 presente. Solo estaba marcado porque la cabecera
lleva el nombre del gen sin etiqueta de columna.

### B.2 Inmunoterapia, antes y después

| GSE | n | Qué |
|---|---|---|
| GSE285963 | 84 | HCC **pre y post atezolizumab + bevacizumab**, con edad y sexo |
| GSE302495 | 70 | HCC **pre y post nivolumab + ipilimumab** (potencialmente resecable) |
| GSE181946 | 17 | HCC de pacientes con **progresión vs respuesta parcial a anti-PD1** |
| GSE287319 | 12 | HCC **pre y post radioterapia con iones de carbono** |

### B.3 Hepatoblastoma y tumores pediátricos

| GSE | n | Qué |
|---|---|---|
| GSE133039 | 66 | **hepatoblastoma** tumor + hígado normal, con sexo |
| GSE306477 | 21 | hepatoblastoma / tumor hepático transicional / **HCC pediátrico** / hígado normal, con genotipo de CTNNB1 y edad en meses |
| GSE301049 | 12 | hepatoblastoma, tumor y adyacente emparejados |

No es HCC del adulto, es otra entidad. Conviene analizarlo aparte, pero es
precisamente la clase de contraste que TCGA-LIHC no tiene.

### B.4 Pares tumor / no-tumor y resto

GSE242315 (71), GSE198946 (48, HBV+ con TP53 y CTNNB1), GSE136711 (41,
**multirregión del mismo tumor** — mide heterogeneidad intratumoral), GSE184733
(34), GSE251942 (25, HBV-HCC: tumor + no-tumor + hígado normal de control),
GSE185700 (24, FFPE: VU_HCC / HCV_HCC / Normal), GSE193567 (18),
GSE146719 (12), GSE70089 (12, hígado **y pulmón** — usar solo las 6 de hígado),
GSE290907 (30), GSE272510 (6), GSE200809 (6, VHC sin tratar vs curado con AAD),
GSE185799 (6), GSE110345 (4), GSE220670 (5).

### B.5 Tres con reservas, pero vale la pena quedárselos

**GSE334651 (120) — RSEM, decimales.** Al abrir el `.xlsx`:

```
gene_id            R1      R10     R101    R102
ENSG00000000003    23.0    13.82   16.32   24.31
```

Son *expected counts* de RSEM, no enteros. Es la excepción a la regla de
"decimales ⇒ no son counts": el camino correcto es redondearlos (lo que hace
`tximport`), no descartarlos. 120 muestras con grupos MASLD-HCC / MASLD / CTRL
valen el trabajo extra. **Hay que documentar el redondeo en la memoria.**

**GSE112705 (40) — la mitad no es RNA-seq.** Declara
`molecule subtype: ribosome-protected fragments` en 20 de las 40 muestras: es
*ribosome profiling*. Usar solo las columnas `totalRNA` (tumor y normal adyacente
de 10 pacientes).

**GSE192862 (9 declaradas, 4 columnas) — microdisección láser.** `liver section
(laser capture)` / `tumor section (laser capture)`, HCV+ y HCV−. Es el equivalente
hepático de GSE93326 en PAAD, pero solo publicaron 4 columnas de las 9 muestras.
Pocas, pero separa epitelio del estroma, que es justo lo que CLDN1 necesita.

---

## C. Comandos de borrado

No he borrado nada. Para aplicarlo:

```bash
cd /home/mike/master-project/dataset/GEO/LIHC
rm -rf GSE98638 GSE292298 GSE278324 GSE237697 GSE276062 GSE82177 \
       GSE120021 GSE254461 GSE144494 GSE317443 GSE255163 GSE245905 \
       GSE92432 GSE289187 GSE285415 GSE221589 GSE120663
```

Libera 661 MB y deja LIHC en ~206 MB.

---

## D. Dos observaciones sobre el pipeline

**1. El caché comparte las sospechas entre cánceres.** GSE148355 llegó con la
marca *"the text points at LIHC, not CHOL"*. Ese motivo se escribió durante la
ejecución de CHOL y quedó guardado en `curation/GEO/inspected.csv`, que es común a
los seis cánceres. Al correr LIHC, la fila se reutilizó tal cual, con un motivo
que ya no aplicaba. **No es un error de decisión** (el dataset se descargó y se
revisó, que es lo que queremos), pero sí hace que algunos motivos del log sean
válidos para otro cáncer. Conviene mencionarlo en la guía antes de citar motivos
textuales en la memoria.

**2. El `.xlsx` es el punto ciego real del filtro.** Seis datasets llegaron
marcados solo por estar en Excel. Al abrirlos: dos eran válidos (GSE181946,
GSE193567), uno tenía decimales recuperables (GSE334651) y **tres eran basura
que ninguna otra puerta habría detenido** (GSE237697 sin columna de genes,
GSE144494 de mama, GSE120021 líneas celulares). La tasa de acierto de esa marca
es del 50 %, la más alta de todas las sospechas blandas.

---

## E. El embudo

```
 972  solo un _RAW.tar por muestra, sin matriz            (48 %)
 251  solo valores normalizados (FPKM/TPM/CPM)            (12 %)
 195  ningún archivo reconocible como matriz              (10 %)
  94  counts "raw" con decimales
  50  los datos solo están en SRA
```

**El 70 % se descarta por cómo se publicó el dato** — frente al 63 % de COAD y el
62 % de PAAD. LIHC es el peor de los tres: casi la mitad de los candidatos
publican solo un `.tar` de archivos por muestra.

Líneas celulares dominantes: **HepG2** (21 series), **Huh7** (19), **Hep3B** (8).
En COAD eran HCT116 y SW620; en PAAD, PANC-1 y MIA PaCa-2.
