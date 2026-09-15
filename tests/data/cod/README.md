# Crystal structures from the Crystallography Open Database

Test fixtures for the CIF route, downloaded unmodified from the
[Crystallography Open Database](https://www.crystallography.net/cod/), whose
data are dedicated to the public domain (CC0).

| file | compound | why it is here |
|---|---|---|
| `7033930.cif` | ferrocene, Fe(C5H5)2 | Fe on an inversion centre (half a molecule in the file), eta-5 ligands, a disordered Cp whose major component is group 2 |
| `1538403.cif` | trans-[PtCl2(NH3)2] | three atoms in the asymmetric unit, Pt on an inversion centre, no hydrogens located, an element beyond Kr |
| `2106540.cif` | [Co(NH3)6]Cl3 | four independent Co on special positions, counter-ions, no hydrogens located |
| `2010630.cif` | Λ-[Co(en)3]Cl3 double salt | helicity reference, Λ by the authors, resolved crystal (Flack 0.07) |
| `2006181.cif` | Λ-[Ru(phen)3](PF6)2 | helicity reference, Λ; also holds disordered solvents |
| `2016683.cif` | Δ-[Ru(phen)3]2+ salt | helicity reference, Δ by the authors |
| `2017108.cif` | rac-[Co(en)2(oxamato)]+ | a racemic crystal (P-1): both enantiomers present |
| `2219943.cif` | (R,Sp)-1-PPh2-2-(1-ethoxyethyl)ferrocene | planar-chirality reference, Sp by the authors (Schlögl) |

The helicity and planar-chirality sign conventions in `m2i/crystal` were fixed
against these labelled structures (and more of their kind; see the docstrings),
not from memory: implemented from memory, the Δ/Λ sign came out inverted.

