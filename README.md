# BSub-aa
Functionality for calculating the amino acid mass fractions of Bacillus subtilis strain 168 from a given protein mass fraction

# Usage
Basic usage in Python

 ```python
  prot_frac = 0.5   # protein mass fraction (g/gDW)
  aa_fracs = get_aa_fractions(prot_frac, mass=False)    # mmol/gDW amino acid, True for g/gDW amino acid
 ```
