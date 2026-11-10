# Selection language

Display rows ([representations](representations.md)) pick their atoms with a short text selection, similar to
VMD's. Selections are typed in the row editor. While you type:
- the atom count updates;
- a list suggests keywords (**↑ / ↓** to choose, **Enter** or **Tab** to accept);
- a mistake is underlined with a message, for example *Unknown keyword "protien" — did you mean "protein"?*.

Keywords, atom names and residue names are case-insensitive. Chain ids are case-sensitive. Quote a name that
contains special characters, for example `name "O5'"`.

## Classes

| Keyword (aliases) | Atoms |
|---|---|
| `all` (`everything`), `none` | Every atom, no atom |
| `protein` | Amino-acid residues |
| `backbone` | Protein N, CA, C, O and their hydrogens; the nucleic-acid phosphate–sugar backbone |
| `sidechain` | Protein atoms that are not backbone |
| `nucleic` (`dna`, `rna`) | Nucleotides |
| `water` (`waters`, `solvent`) | Water molecules (HOH, WAT, SOL, TIP3…) |
| `ions` (`ion`) | Monatomic ions (NA, CL, SOD, CLA, POT, CAL…) |
| `lipid` (`lipids`) | Lipids and sterols (POPC, CHL1, PSM, Lipid21 PA/PC/OL…) |
| `carbohydrate` (`carbohydrates`, `sugar`, `sugars`, `glycan`, `glycans`) | Monosaccharide residues |
| `ligand` (`ligands`) | Everything else: hetero groups, small molecules, cofactors |
| `hydrogen` (`hydrogens`) | Hydrogen atoms |
| `heavy` | Every atom except hydrogen |
| `metal` (`metals`) | Metal atoms (everything that is not a nonmetal or metalloid) |
| `helix`, `sheet` (`strand`), `coil` | Protein residues by secondary structure |

How residues are put into classes is described in [Representations](representations.md#how-residues-are-recognised).

## QM/MM layers and flags

| Keyword (aliases) | Atoms |
|---|---|
| `layer H`, `qm`, `high` | The high layer or QM region |
| `layer M`, `qm2`, `medium` | The medium layer or QM2 region |
| `layer L`, `mm`, `low` | The low layer or MM region |
| `layer H M` | Several layers at once |
| `link` | QM/MM link atoms |
| `frozen`, `active` | Frozen atoms and active (not frozen) atoms |
| `selected` (`selection`) | The atoms selected in the 3D view |

## Properties with values

| Keyword (aliases) | Values | Example |
|---|---|---|
| `chain` (`segid`, `segname`) | Chain or segment ids | `chain A B` |
| `resname` (`resn`) | Residue names; wildcards `*` and `?` | `resname POPC CHL1`, `resname PO*` |
| `resid` (`resi`, `resnum`) | Residue numbers and ranges | `resid 10-20 35`, `resid 10 to 20` |
| `name` | Atom names; wildcards | `name CA N C O`, `name H*`, `name C1?` |
| `element` (`elem`) | Element symbols | `element C N`, `element Pt` |
| `index` | 0-based atom indices and ranges | `index 0-99` |
| `serial` | 1-based atom numbers and ranges | `serial 1-100` |

Ranges can be written as `10-20`, `10:20`, `10..20` or `10 to 20`. Several values mean "any of these".

## Operators

| Syntax | Meaning |
|---|---|
| `A and B` | Atoms in both |
| `A or B` | Atoms in either |
| `not A` | Atoms not in A |
| `( … )` | Grouping |
| `within D of A` | Atoms within D Å of any atom of A (A itself included) |
| `same residue as A` | Whole residues that contain at least one atom of A |

`not` binds most tightly, then `and`, then `or`. So `protein and not hydrogen or ligand` means
`(protein and (not hydrogen)) or ligand`. When in doubt, use parentheses. `within` and `same residue as` apply to
the selection that follows them. Write `within 5 of (ligand or ions)` to give them a combination.

## Examples

```text
protein and not hydrogen                       heavy atoms of the protein
same residue as within 5 of ligand             the binding site, whole residues
water and within 3.5 of ligand                 water molecules close to the ligand
resname POPC and name P                        phosphorus atoms of POPC
lipid and not (name H* or element H)           lipids without hydrogens
chain A and resid 100-150 and backbone         a stretch of backbone
carbohydrate or (protein and within 6 of carbohydrate)
                                               sugars and the protein around them
layer H or link                                QM region and its link atoms
layer L and within 8 of layer H                MM atoms near the QM region
same residue as selected                       the residues of the atoms you clicked
element Pt or (within 2.6 of element Pt)       a Pt centre and the atoms bonded to it
```
