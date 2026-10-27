# Data for the paper

High-throughput screening data behind the revised paper *Accelerated Pd-catalysed cross-coupling optimisation with compressed sensing*. The data were generated at F. Hoffmann-La Roche Ltd and exported from its screening database. The scripts in [`reproduction/paper`](../../reproduction/paper) use exactly these files; [`reproduction/README.md`](../../reproduction/README.md) explains how to run them.

## Files

| File | Content |
|---|---|
| `pooled_plates/ELN032036-334_plate-23_Buchwald-Hartwig_raw.csv` | Pooled CatMix plate of the sterically hindered Buchwald--Hartwig amination (3 sampling points) |
| `pooled_plates/ELN032036-424_plate-13_Meldrum_raw.csv` | Pooled plate of the arylation of Meldrum's acid (2 sampling points) |
| `pooled_plates/ELN032036-425_plate-11_Barbituric_raw.csv` | Pooled plate of the arylation of barbituric acid (2 sampling points) |
| `pooled_plates/ELN032036-810_plate-1_Barluenga_raw.csv` | Pooled plate of the three-component Barluenga coupling (3 sampling points) |
| `pooled_plates/ELN032036-821_plate-1_Suzuki_raw.csv` | Pooled plate of the late-stage Suzuki--Miyaura coupling (3 sampling points) |
| `pooled_plates/manifest.json` | SHA-256 of each pooled-plate file, and of the complete export it was taken from |
| `reference_334_single_catalyst_plates_15-18.csv` | Single-catalyst reference plates of the amination: 43 precatalysts at 8 mol%, one per well, under eight solvent/base combinations (ELN032036-334, plates 15 to 18) |
| `titration_334_plate_26.csv` | XantPhos titration plate of the amination (ELN032036-334, plate 26, 3 sampling points) |
| `validation_yields.csv` | Externally calibrated LC--MS assay yields (%) of the catalysts validated individually in each campaign; conditions A and B for campaigns 424 and 425 are those defined in the paper |

The numbers in the file names are the electronic-notebook references used in the paper (ELN032036-334 is campaign 334, and so on).

## Columns

| Column | Meaning |
|---|---|
| `ELN_ID` | Electronic-notebook reference of the campaign |
| `SAMPLE_ID` | Identifier of the analysed sample (one well at one sampling point) |
| `PLATENUMBER` | Plate number within the notebook |
| `SAMPLENUMBER` | Sampling point (1, 2, 3, ...) |
| `PLATEROW`, `PLATECOLUMN` (`PlateRow`, `PlateColumn`) | Well position |
| `COMPOUND` | Classification of the LC--MS peak: `Product`, `lim SM` (limiting starting material), `other SM`, `unknown`, and so on |
| `AREA_TOTAL_REDUCED` | Reduced total area percentage of the peak (SI section S1.2) |
| `ColumnLabel`, `RowLabel` | Plate labels: the solvent/base combination of a column (pooled plates) or the precatalyst of a column and the solvent/base combination of a row (reference plates) |

A pooled-plate file contains every row of the export (every peak, classification, sampling point and well) but only the nine columns the analysis reads. The complete exports also contain instrument and inventory fields that the analysis does not use; their SHA-256 hashes are listed in `manifest.json` as `full_export_sha256`. The reference and titration files contain the rows of the stated plates and the columns shown. In the export the reference rows come from, the `ELN_ID` column held a copy of the area values; it has been replaced by the notebook reference.
