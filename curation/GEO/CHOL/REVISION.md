# CHOL — revisión de los datasets descargados

Recomendación de Claude tras verificar los registros de muestra de cada serie.
**No borra nada**: la decisión es tuya. Fecha: 2026-10-06.

Veredicto del pipeline: `204 encontrados -> 57 pasaron puertas 1-2 -> 6 confirmados + 4 ambiguos`

---

## Conservar (5 datasets, 332 muestras)

| Dataset | Muestras | Lo que declaran las muestras | Por qué |
|---|---|---|---|
| **GSE244807** | 246 | `surgical specimen` | La cohorte más grande. Ensembl, 246 columnas |
| **GSE107943** | 57 | `tumor (intrahepatic cholangiocarcinoma)` | Tumor + normal adyacente pareados. `Sex`, `age`, `dsfree(mo)` |
| **GSE162396** | 12 | `cholangiocarcinoma tissue` | Columnas `CCC-1..CCC-12` |
| **GSE308395** | 10 | `tissue: gallbladder` | Recuperado al pasar las marcas del abstract a sospecha |
| **GSE244331** | 7 | `tumor tissue` | Columnas `LLHM-T` / `LLHM-N`, pareado |

---

## Borrar — 5 datasets

```bash
cd ~/master-project/dataset/GEO/CHOL
rm -rf GSE306865 GSE261998 GSE221589 GSE160306 GSE148355
```

### GSE306865 — FALSO POSITIVO del filtro

El pipeline lo marcó `tissue confirmed`. **No lo es.** Sus registros:

```
tissue: tumor samples from korean patients with icc     <- describe el ORIGEN
cell line: sck
cell line: sck-r (cisplatin-resistant sck)
cell type: human intrahepatic cholangiocarcinoma cells
```

Serie **mixta**: el campo `tissue:` habla del origen de la línea celular, no de las
muestras. Sus 4 muestras son SCK vs SCK-R (resistente a cisplatino), como confirma su
archivo `WT_vs_R.tsv.gz`.

**Causa raíz, para la memoria:** la regla que distingue el nombre de una línea de una
frase exige que haya un dígito. `sck` no lo tiene, así que no se reconoció como nombre y
ganó el campo `tissue:`. Las líneas sin dígito (SCK, HeLa, AGS) son un punto ciego
conocido; `RBE` y `NOZ` se salvan solo porque están en la lista de nombres.

### GSE261998 — xenoinjerto, 10 muestras

```
tissue: hucc-a1 cdx tumors
treatment: cldn1 mab
```

`CDX` = *cell-line-derived xenograft*. No es biopsia.

> **Pero guarda la referencia.** Es un estudio de **anticuerpo monoclonal anti-CLDN1** en
> colangiocarcinoma: exactamente la proteína del proyecto, en contexto terapéutico. Los
> datos no sirven para la matriz; el paper es muy relevante. Coméntaselo a Abrar.

### GSE221589 — modelo de ingeniería genética, 27 muestras

```
tissue: liver tumor
genotype: nras, akt3, arid1a
genotype: sv40tt, tert
```

Título: *"Human hepatocytes can give rise to intrahepatic cholangiocarcinoma"*. Son
hepatocitos transformados con oncogenes, no tumores de paciente.

### GSE160306 — retina, 79 muestras

```
tissue: retina macula
```

Retinopatía diabética. No tiene relación con vías biliares. Apareció porque su resumen
coincidía con palabras de la búsqueda.

**Limitación detectada:** la lista `OTHER_ORGANS` no incluye ojo ni retina. El sistema lo
marcó ambiguo (no lo declaró limpio), así que se comportó como debía.

### GSE148355 — pertenece a LIHC, 128 muestras

Título: *"Preoperative immune landscape... in hepatocellular carcinoma patients with liver
transplantation"*. Sus muestras incluyen `normal tissue, liver`.

**No lo pierdes:** aparecerá correctamente al ejecutar LIHC, donde sí encaja. Bórralo solo
de la carpeta de CHOL para no duplicarlo.

---

## Rechazos que verifiqué y están bien

Por si al revisar `decisions.csv` alguno te extraña:

| Dataset | Mensaje confuso | Realidad |
|---|---|---|
| GSE179599 | `declare culture (ORGANOID): cca tumor tissue` | Tiene **ambos**: `cca tumor tissue` y `cholangiocarcinoma tumor organoids` + `cultured in expansion medium`. Organoides, correctamente fuera |
| GSE244424 | `declare culture (CELL_LINE): human biliary tissues` | `cell line: hibepic` con `ppdpf overexpressed`. Células epiteliales biliares con sobreexpresión. Correcto |

El mensaje muestra el primer valor de los registros, que puede contradecir el veredicto.
Es confuso al revisar, no un error de clasificación.

---

## Balance del cambio de filtro

Esta ejecución usó la versión donde las marcas del abstract son sospecha y no rechazo.

| | Antes | Ahora |
|---|---|---|
| Confirmados correctos | 4 | **5** |
| Falsos positivos | 0 | 1 (GSE306865) |
| Precisión | 100% | 83% |

Perdí precisión y gané un dataset (GSE308395, 10 muestras de vesícula biliar).

Mi criterio: buen intercambio. Un falso positivo entre 6 se detecta en dos minutos
abriendo su ficha; un dataset perdido en silencio no se detecta nunca.
