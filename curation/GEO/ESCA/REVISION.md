# ESCA — revisión de los datasets descargados

Recomendación de Claude tras verificar los registros de muestra de cada serie.
**No borra nada**: la decisión es tuya. Fecha: 2026-10-06.

Veredicto del pipeline: `317 encontrados -> 56 pasaron puertas 1-2 -> 5 confirmados + 3 ambiguos`

> **Contexto importante.** La ejecución anterior de ESCA daba **0 confirmados**. El cambio
> que pasó las marcas del abstract de rechazo a sospecha recuperó 5 datasets, y los he
> verificado uno a uno: **los 5 son tejido genuino**. Sin ese cambio, ESCA se habría
> quedado sin ningún dato de GEO.

---

## Conservar (7 datasets, 152 muestras)

| Dataset | Muestras | Lo que declaran las muestras | Valor |
|---|---|---|---|
| **GSE234304** | 44 | `tissue: tumor`, `non-tumor`, **`barretts esophagus`** | El más interesante: incluye la **lesión precursora** |
| **GSE273848** | 40 | `tissue: tumor`, `treatment: yes neoadj`, `disease stage: relapse` | Adenocarcinoma con anotación clínica y de tratamiento |
| **GSE235537** | 22 | `tissue: esophageal squamous cell carcinoma`, `patient: patient 1..` | **11 pares tumor/normal** (`counts-11pairs`) |
| **GSE305720** | 20 | `tissue: esophageal squamous cell carcinoma` | Tumor primario T2 con y sin metástasis ganglionar |
| **GSE194116** | 12 | `tissue subtype: escc tumor` / `adjacent normal tissue` | Pareado tumor/normal |
| **GSE119436** | 8 | `tissue: normal esophagus` + `age`, `gender`, `stage` | Con datos clínicos |
| **GSE157373** | 6 | `tissue: escc tissue` + `adjacent normal tissues` | Pareado, pero ver la advertencia |

### Advertencia sobre GSE157373

Sus muestras son tejido pareado legítimo y sus identificadores son símbolos de gen
reales (verificados: 34 de 39 coinciden). Pero el estudio es sobre **redes ceRNA de
circRNA y lncRNA**.

Eso importa por una razón técnica: los estudios de RNA no codificante suelen usar **RNA
total con depleción de rRNA**, no selección por poly-A. Cambia sistemáticamente el perfil
de cuantificación a nivel de gen. Al fusionarlo con TCGA (poly-A) tendrás un efecto de
lote adicional sobre este dataset concreto.

Con 6 muestras, mi opinión es que aporta poco frente al riesgo. Lo pongo en "conservar"
porque es técnicamente válido, pero sería el primero que yo descartaría si quieres
simplificar.

---

## Borrar — 1 dataset

```bash
cd ~/master-project/dataset/GEO/ESCA
rm -rf GSE291098
```

### GSE291098 — panel multiorgánico, 10 muestras

```
tissue: esophagus tumor
tissue: bladder
tissue: bone
```

Título: *"RNA Transcripts Serve as a Template for Double-Strand Break Repair"*. No es un
estudio de cáncer de esófago: es un trabajo de **reparación del ADN** que usa un panel de
tejidos de varios órganos, del cual el esófago es solo una parte.

Apareció porque la búsqueda coincidió con palabras del resumen. El sistema lo marcó
ambiguo (`the abstract names no cancer type clearly`), que es lo correcto.

**Nota técnica:** la regla de órganos ajenos no lo rechazó en firme porque **el esófago
también aparece** — es el comportamiento que se diseñó para no perder las metástasis de
nuestro cáncer. Aquí el efecto secundario es que deja pasar un panel multiorgánico como
ambiguo. Funcionó como debía; el caso es el que hay que decidir a mano.

---

## Precisión de esta ejecución

| | Resultado |
|---|---|
| Confirmados correctos | **5 de 5** |
| Ambiguos que sí valen | 2 de 3 |
| Falsos positivos | 0 |

Mejor que en CHOL, donde hubo 1 falso positivo (GSE306865). La diferencia: en ESCA
ninguna serie tenía un nombre de línea celular sin dígito, que es el punto ciego conocido.

---

## Comparación con la ejecución anterior

| | Antes del cambio | Ahora |
|---|---|---|
| Confirmados | 0 | 5 |
| Ambiguos | 1 | 3 |
| Muestras de tejido | 0 | **152** |

Los 5 recuperados estaban siendo rechazados porque su resumen mencionaba líneas
celulares. Verificado que ninguno lo es.

Dato para la memoria: de los 317 candidatos de esófago, **104 se descartaron por publicar
solo un `_RAW.tar` por muestra** sin matriz de counts, y **41 por publicar solo valores
normalizados**. Esos dos motivos técnicos, no biológicos, eliminan el 46% del total.
