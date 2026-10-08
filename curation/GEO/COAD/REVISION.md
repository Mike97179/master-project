# COAD — revisión manual de los datasets descargados

Fecha de la ejecución: 2026-10-06 23:59 → 2026-10-07 02:31 (2 h 31 m)
Embudo: **2080 candidatos → 457 pasan puertas 1-2 → 36 confirmados + 35 ambiguos = 71 descargados (749 MB)**

Verificación hecha a mano sobre `sample_metadata.tsv` (lo que el autor declaró en
cada GSM) y sobre la cabecera real de cada matriz.

**Resultado: conservar 43 datasets (~1543 muestras), borrar 28 (578 MB de 749).**

---

## A. BORRAR — 28 datasets, 578 MB

### A.1 Líneas celulares y organoides que pasaron la puerta 4 (7)

Estas son las fallas reales del filtro en COAD. Todas declaran `tissue: colon`
o similar **y a la vez** una línea celular que el detector no reconoció.

| GSE | n | Qué es en realidad | Por qué pasó |
|---|---|---|---|
| GSE328387 | 6 | `tissue: colon; cell line: RKO` | **RKO no tiene dígitos** — el detector exige ≤2 palabras *más un dígito* |
| GSE289480 | 40 | `cell line: RKO; cell type: human colorectal carcinoma` | igual que el anterior |
| GSE304081 | 12 | `cell line: SC; cell type: Colorectal Cancer` + tratamientos con fármacos | **"SC" tampoco tiene dígitos** |
| GSE140984 | 12 | `cell line: HT-29 (ATCC HTB-38); tumor type: colorectal` | el valor tiene **3 tokens** por el catálogo ATCC, el detector pedía ≤2 |
| GSE178698 | 12 | organoides de colon dañados crónicamente (el título lo dice) | los GSM solo dicen `region: transverse colon`, sin la palabra organoide |
| GSE330032 | 6 | `KRAS G13C mutant CRC organoids` (en el *título* de la muestra) | los `characteristics` solo dicen `tissue: CRC` |
| GSE169720 | 27 | fibroblastos CAF en cultivo (`CAF1_1`, `CAF1_2`…) | declara `cell type: colon cancer associated fibroblasts` |

El patrón es claro y es nuevo: **el nombre de la línea celular sin dígito**
(RKO, SC) y **el nombre con el código de catálogo** (HT-29 (ATCC HTB-38)).
Era la limitación que ya estaba documentada en `docs/03_geo_guide.md`; COAD la
expone cuatro veces porque RKO y SW/HT son las líneas estándar del campo.

### A.2 Tejido equivocado (6)

| GSE | n | Qué es |
|---|---|---|
| GSE160306 | 79 | **retina** (`tissue: Retina Macula`) — tercera vez que aparece, ya salió en CHOL y PAAD |
| GSE254461 | 57 | **gliomas** (oligodendroglioma, glioblastoma, meningioma) |
| GSE253560 | 70 | lesiones escamosas **anales** (HGSIL / LGSIL / ASCC), no colorrectales |
| GSE162945 | 20 | metástasis de **cérvix y pulmón** irradiadas |
| GSE186466 | 6 | **tejido adiposo** subcutáneo y visceral de pacientes con caquexia |
| GSE288482 | 52 | pan-cáncer (suprarrenal, mama, urotelial…), solo una parte es colorrectal |

### A.3 Sangre y células inmunes aisladas (6)

CLDN1 es una proteína epitelial de unión estrecha. En estas muestras no hay epitelio.

| GSE | n | Qué es |
|---|---|---|
| GSE289582 | 115 | CD4 y CD8 purificados de **sangre periférica** |
| GSE302921 | 50 | **evRNA de plasma** (es el par de GSE302922, que sí es tejido) |
| GSE310282 | 11 | **Treg** sorteados de colon (epitelio / lámina propia) |
| GSE247742 | 14 | **ILC3** sorteados de colon |
| GSE165814 | 8 | **ILC3** sorteados de colon |
| GSE265950 | 12 | **macrófagos** derivados de PBMC de donantes sanos |

### A.4 Single-cell disfrazado de bulk (4) — 507 MB, el 68 % del disco

| GSE | n | Tamaño | Qué es |
|---|---|---|---|
| GSE144735 | 18 | **311 MB** | 10x single-cell de 6 pacientes belgas (cohorte KUL3) |
| GSE150050 | 28 | 98 MB | pools de 384 células de sangre, pulmón, hígado |
| GSE160686 | 18 | 78 MB | tumor disociado tratado *ex vivo* 24 h con anti-TGFβ |
| GSE92432 | 10 | 1 MB | paper de método para aislar células individuales |

GSE144735 es, por sí solo, el 42 % de todo lo que COAD ocupa en disco.

### A.5 Otros (5)

| GSE | n | Motivo |
|---|---|---|
| GSE230524 | 226 | **la matriz no tiene identificadores de gen**: las filas son `1, 2, 3…` y la cabecera solo trae 35 nombres de muestra para 226 GSM declarados. Inutilizable sin la anotación, que el autor no publicó |
| GSE280311 | 12 | 10x Flex multiplexado; la matriz real estaba en el bundle de **522 MB que el pipeline saltó bien**, así que en disco solo quedan el README y la referencia de features |
| GSE288086 | 9 | línea celular **BJAB** (linfoma B), paper de inhibidores de RelB |
| GSE192400 | 6 | macrófagos **THP-1** co-cultivados con *F. nucleatum* |
| GSE147971 | 17 | **organoides** derivados de paciente (cribado de fármacos) |

GSE230524 merece una nota: es el único caso en los cinco cánceres en que la
sospecha blanda *"the first column does not look like gene identifiers"* acertó
de lleno. Lo marcó como ambiguo, se descargó, y al abrirlo no hay genes.
Justifica la decisión de descargar los dudosos en lugar de descartarlos a ciegas.

---

## B. CONSERVAR — 43 datasets, ~1543 muestras

### B.1 Los que valen más

**GSE262671 (113) + GSE262672 (33) + GSE262670 (25) — 171 muestras, cáncer colorrectal primario.**
Los tres son el proyecto *NanoCMSer* (clasificador de subtipos moleculares
consensus, CMS). `tissue: Primary colorectal cancer`, fresh frozen, con sexo y
fijación declarados. Los tres quedaron ambiguos **solo porque los counts están
en `.xlsx`**, que no se puede leer por rango de bytes. Los abrí a mano:

```
gene             CO129  CO157  CO123 ...
ENSG00000223972    ...    ...    ...
```

Identificadores Ensembl, enteros, y el número de columnas cuadra exactamente con
las muestras declaradas (26-1=25, 114-1=113, 34-1=33). **Son válidos.** Es el
mayor hallazgo de COAD y, en número de muestras, el segundo de todo GEO después
de GSE93326 en PAAD.

**GSE290622 (102) — tumor colorrectal bulk.** Títulos `JM209B,tumor`,
`tissue: tumor; cell type: CRC`. 102 tumores.

**GSE166925 (96) — tejido intestinal completo de pacientes con EII.**
UC y Crohn, inflamado y no inflamado, intestino grueso y delgado, con sexo y
paciente. No es cáncer, pero es mucosa colónica humana bien anotada: sirve como
referencia de inflamación frente a tumor.

**GSE183202 (90) — metástasis peritoneal y tumor primario** emparejados.

**GSE183150 (60) — 60 muestras de mucosa colorrectal normal** de un ensayo de
ejercicio físico. TCGA-COAD tiene pocos normales; esto los complementa.

**GSE171468 (59) — tejido normal adyacente** al tumor, con sexo, edad, IMC y
toxicidad al tratamiento.

**GSE179351 (54)** — ensayo fase II, `tissue: metastatic tumor` con
`cancer type: MSS CRC` **y** `cancer type: PDAC` en el mismo dataset. Cuenta
doble: sirve para COAD y para PAAD.

### B.2 Gradiente normal → adenoma → carcinoma

Lo más útil para medir CLDN1 a lo largo de la progresión:

| GSE | n | Qué separa |
|---|---|---|
| GSE164541 | 15 | epitelio colónico normal / adenoma / CRC primario |
| GSE95132 | 31 | normal adyacente / tumor primario / cripta normal / **foco de cripta aberrante** |
| GSE110536 | 13 | colon / pólipo serrado sésil / adenoma serrado tradicional / PAF / carcinoma |
| GSE222298 | 36 | adenoma / microadenoma / mucosa sana (quimioprevención en PAF) |
| GSE200130 | 37 | no neoplásico / peritumoral / tumor |

### B.3 Pares tumor-normal y metástasis

GSE196006 (42, CRC de inicio precoz, emparejado), GSE152430 (49, FFPE con
estadio y localización), GSE251845 (44, por segmento del colon), GSE242676 (48,
metástasis peritoneal), GSE264227 (42, normal/tumor + tratado con inhibidor de
SMYD3), GSE256268 (42), GSE274551 (35, recto/hígado/pulmón/ganglio con CMS y
RECIST), GSE66207 (33, biopsias de colon en Crohn), GSE255163 (30, primario +
metástasis hepática), GSE211831 (30), GSE303444 (30, CRC con MMR deficiente y
tipo histológico), GSE151165 (30, tumor/normal en metástasis hepática),
GSE213331 (27, recto localmente avanzado pre/pCR/npCR), GSE145429 (22),
GSE200407 (22, sensibilidad a FOLFOX), GSE233517 (22, pre/post radioterapia),
GSE280384 (17), GSE292184 (17, **metástasis cerebral** y hepática — poco común),
GSE138202 (16), GSE242786 (15, recto, respuesta completa vs incompleta),
GSE289583 (14), GSE121842 (12), GSE302922 (12), GSE235850 (10), GSE179979 (7),
GSE156172 (4), GSE200427 (4).

---

## C. Comandos de borrado

No he borrado nada. Para aplicarlo:

```bash
cd /home/mike/master-project/dataset/GEO/COAD
rm -rf GSE144735 GSE150050 GSE160686 GSE288086 GSE160306 GSE330032 \
       GSE289582 GSE186466 GSE265950 GSE254461 GSE288482 GSE253560 \
       GSE192400 GSE328387 GSE289480 GSE162945 GSE92432 GSE310282 \
       GSE304081 GSE302921 GSE280311 GSE247742 GSE178698 GSE169720 \
       GSE165814 GSE147971 GSE140984 GSE230524
```

Libera 578 MB y deja COAD en ~171 MB.

---

## D. Nota sobre el embudo

De los 2080 candidatos, los motivos de descarte más frecuentes:

```
 820  solo un _RAW.tar por muestra, sin matriz de counts     (39 %)
 251  solo valores normalizados (FPKM/TPM/CPM)               (12 %)
 243  ningún archivo reconocible como matriz                 (12 %)
  76  counts "raw" con decimales
  55  marca de tipo de muestra en el nombre del archivo: CELL_LINE
  54  los datos solo están en SRA
```

**El 63 % se descarta por cómo se publicó el dato, no por su contenido.**
En PAAD fue el 62 %. La cifra se repite en los dos cánceres con más candidatos,
lo que la convierte en un dato defendible para la memoria.

De los que llegaron a la puerta 4, la línea celular más repetida fue **HCT116**
(43 series) seguida de **SW620** (14) y **HT-29** (9). En PAAD eran PANC-1 y
MIA PaCa-2. Cada cáncer tiene sus dos o tres líneas dominantes.

También hubo **un error 502 de NCBI** (GSE154429, candidato 1561). Queda sin
inspeccionar; el caché es incremental, así que una nueva ejecución lo retomaría.
