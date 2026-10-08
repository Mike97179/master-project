# STAD — revisión manual de los datasets descargados

Fecha de la ejecución: 2026-10-07 14:33 → 15:08 (34 m 35 s)
Embudo: **485 candidatos → 95 pasan puertas 1-2 → 6 confirmados + 5 ambiguos = 11 descargados (129 MB)**

**Resultado: conservar 5 datasets (**853 muestras**, suma verificada), borrar 6 (53 MB de 129).**

STAD es el cáncer con menos candidatos (485 frente a ~2000 de COAD y LIHC) y con
la tasa de descarte más alta por líneas celulares. Pero los cinco que sobreviven
son excelentes: **853 muestras, todas tejido gástrico humano con estadio o diagnóstico.**

---

## A. CONSERVAR — 5 datasets, 853 muestras

Los cinco los verifiqué abriendo la matriz. En los cinco el número de columnas
cuadra **exactamente** con las muestras declaradas, los counts son enteros y
**CLDN1 está presente**:

| GSE | n | Columnas | Identificadores | Qué es |
|---|---|---|---|---|
| **GSE184336** | **461** | 462−1 ✓ | Ensembl, 58 736 genes | tumor gástrico + **normal emparejado**, con estadio (ⅠA–ⅢB) y sexo |
| **GSE275647** | **216** | 217−1 ✓ | Ensembl versionado | estómago, **cáncer gástrico avanzado (AGC)**, con edad, sexo y diagnóstico |
| **GSE158631** | 94 | 95−1 ✓ | símbolos, 21 196 genes | cáncer gástrico metastásico: **tumor y ganglio linfático** del mismo paciente |
| **GSE179252** | 76 | 77−1 ✓ | Ensembl | tumor gástrico + normal emparejado, con estadio, edad y sexo |
| **GSE277340** | 9 | 9 + anotación ✓ | Ensembl + símbolo | **estómago normal** / enfermedad de Ménétrier / poliposis juvenil |

### GSE184336 es el mayor hallazgo de todo el esfuerzo en GEO

**461 muestras.** Supera a GSE263786 (243, LIHC) y a las 206 de las dos cohortes
de Taiwán juntas. Tumor gástrico con normal emparejado y estadio declarado por
muestra. Para un cáncer con solo 485 candidatos en total, es un resultado
desproporcionado.

GSE179252 es del **mismo grupo y mismo título** ("Next Generation Sequencing
Facilitates Quantitative Analysis…"): es la misma cohorte publicada en dos tandas.

**Comprobado cruzando los nombres de columna de las dos matrices:**

```
GSE179252:  B1   … B76    (76 columnas)
GSE184336:  B77  … B557   (461 columnas)

columnas en común: 0
```

Cero solapamiento y rangos consecutivos: uno acaba en `B76`, el otro empieza en `B77`.
**Las 537 muestras suman correctamente.** (El rango de GSE184336 llega a `B557` con 461
columnas, así que hay 96 huecos en la numeración — muestras que no pasaron el control de
calidad y no se publicaron. No afecta a nada.)

### GSE277340: pocas muestras, mucho valor

Solo 9, pero son **estómago normal** (`condition: Normal`) frente a dos
enfermedades gástricas hiperproliferativas no tumorales. TCGA-STAD tiene 36
normales; esto añade controles de otra procedencia y, sobre todo, un contraste
proliferativo sin cáncer.

---

## B. BORRAR — 6 datasets, 53 MB

### B.1 Xenoinjertos en ratón (2)

| GSE | n | Qué delata |
|---|---|---|
| GSE164961 | 53 | `host strain: NOD.Cg-Prkdcscid Il2rgtm1Wjl/SzJ` y títulos `SNU-JAX-G001 (PDX)` — **PDX en ratón NSG**. La puerta 3 ya lo había marcado: 16 columnas para 53 muestras declaradas |
| GSE262056 | 6 | el título lo dice entero: *"gastric cancer **PDX mice model**"*. Estaba entre los **confirmados**, porque los GSM solo declaran `tissue: gastric cancer` |

GSE262056 es la falla del filtro en STAD: la palabra "PDX" aparece en el título
de la serie pero no en los registros de muestra, y la puerta 1 trató la mención
como sospecha blanda que la puerta 4 luego levantó al confirmar "tejido".

### B.2 Tejido equivocado (3)

| GSE | n | Qué es |
|---|---|---|
| GSE289743 | 83 | **carcinoma escamoso cutáneo** (piel), anti-PD-1. Su marca decía *"the text points at PAAD, not STAD"* — la marca también se equivocaba, no es páncreas |
| GSE212521 | 34 | trazador FAP en humanos: vejiga, hígado, pulmón, ganglio. Metástasis de varios orígenes, **ninguna muestra de estómago primario**, y la mitad con `disease: na` |
| GSE288482 | 52 | pan-cáncer (suprarrenal, mama, urotelial…). Ya recomendé borrarlo en COAD |

### B.3 Duplicado de otro cáncer (1)

**GSE148355 (128, 34 MB)** — es el dataset **hepático** de progresión
normal→fibrosis→cirrosis→HCC que recomiendo conservar en LIHC. Aquí llegó por
arrastre del caché y no tiene nada que ver con estómago. Se queda en LIHC.

---

## C. Comandos de borrado

No he borrado nada. Para aplicarlo:

```bash
cd /home/mike/master-project/dataset/GEO/STAD
rm -rf GSE148355 GSE289743 GSE212521 GSE164961 GSE288482 GSE262056
```

Libera 53 MB y deja STAD en ~76 MB.

---

## D. Dos cosas que conviene anotar

**1. El caché cruzado vuelve a aparecer, y ahora dos veces.** GSE148355 llegó con
*"the text points at LIHC, not CHOL"* y GSE289743 con *"the text points at PAAD,
not STAD"*. Ambos motivos se escribieron en ejecuciones anteriores y quedaron en
`curation/GEO/inspected.csv`, que es común a los seis cánceres. El segundo además
era incorrecto en origen: GSE289743 es piel, no páncreas.

La consecuencia práctica: **los motivos textuales del log no son citables sin
comprobarlos**, porque pueden venir de otro cáncer. Las decisiones siguen siendo
correctas (se descargan y se revisan), pero la trazabilidad del *por qué* se
degrada. Vale la pena mencionarlo en `docs/03_geo_guide.md`.

**2. GSE179252 se bajó un GTF de 42 MB.** `GSE179252_Homo_sapiens_Ensemble_94.gtf.gz`
es la anotación del genoma, no datos. Son 42 de los 45 MB del dataset. No es un
error — el autor lo publicó como suplementario y de hecho sirve para mapear los
Ensembl a símbolos — pero es el único caso en los seis cánceres en que la mayoría
del peso de un dataset es un archivo de referencia.

---

## E. El embudo, y por qué STAD es distinto

```
159  solo un _RAW.tar por muestra, sin matriz        (33 %)
 74  solo valores normalizados (FPKM/TPM/CPM)        (15 %)
 73  ningún archivo reconocible como matriz          (15 %)
 36  solo counts de circRNA/isoform/miRNA            (7 %)
 15  los datos solo están en SRA
```

**El 70 % se descarta por cómo se publicó el dato**, igual que LIHC.

Lo llamativo es la puerta 4: de los 95 que llegaron a ella, **la inmensa mayoría
declara una línea celular**. Las dominantes son **AGS** (7 series), **HGC-27**
(6 entre `hgc-27` y `hgc27`) y **MKN45** (4). Aparecen también SGC-7901, MGC-803,
SNU-668, SNU-719, NUGC-4, KATO III, NCI-N87 y YCCEL1.

Dicho en claro: **de 485 candidatos en cáncer gástrico, solo 5 publican bulk
RNA-seq de tejido humano con counts enteros.** Es la proporción más baja de los
seis cánceres (1 %) y el argumento más contundente que tengo para la memoria
sobre los límites de reutilizar datos públicos.

El contraste con el resultado también es llamativo: esos 5 datasets aportan
**853 muestras**, frente a las 407 de TCGA-STAD.
