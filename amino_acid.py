import os
import requests
from requests.exceptions import HTTPError
from Bio import SeqIO

AMINO_ACIDS_THREE = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 
                     'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 
                     'Thr', 'Trp', 'Tyr', 'Val']
AMINO_ACIDS_ONE = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 
                   'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
MOLAR_MASSES = [89.09, 175.21, 132.12, 132.09, 121.16, 146.12, 146.14, 75.07,
              155.15, 131.17, 131.17, 147.19, 149.21, 165.19, 115.13, 105.09,
              119.12, 204.22, 181.19, 117.15]


def get_aa_fractions(protein, mass=True):
    """
    Given the protein biomass mass fraction (g/gDW), calculate the corresponding
    amino acid mass or molar fractions. The relative amino acid distribution is
    calculated from the protein-coding genes of the genome of Bacillus subtilis
    strain 168 (UniProt Proteome ID UP000001570).

    Parameters
    ----------
    protein : float
        Protein biomass mass fraction (g/gDW)
    mass : bool
        True for amino acid mass fractions (g/gDW, default), False for molar 
        fractions (mmol/gDW)
    
    Returns
    -------
    aa_fracs : dict[str, float]
        Dictionary with three letter amino acid codes as keys, and amino acid
        mass/molar fractions as values
    """
    # Download B. subtilis proteome
    filename = 'bsu.fasta'
    if not os.path.isfile(filename):
        url = ('https://rest.uniprot.org/uniprotkb/stream?format=fasta&query'
               '=%28%28proteome%3AUP000001570%29%29')
        try:
            response = requests.get(url)
            response.raise_for_status()
        except HTTPError as err:
            print(f'HTTP error occured (status code: {err.response.status_code})')
        with open(filename, 'w') as handle:
            handle.write(response.text)
    
    # Calculate relative amino acid mass fractions
    aa_fracs = dict.fromkeys(AMINO_ACIDS_ONE, 0)
    g_per_mol = dict(zip(AMINO_ACIDS_ONE, MOLAR_MASSES))
    total_mass = 0
    with open(filename, 'r') as handle:
        for entry in SeqIO.parse(handle, 'fasta'):
            for aa in aa_fracs:
                aa_mass = entry.seq.count(aa) * g_per_mol[aa]
                total_mass += aa_mass
                aa_fracs[aa] += aa_mass

    # Calculate mass fractions from given protein amount
    [aa_fracs.update({aa: protein * (aa_fracs[aa] / total_mass)}) 
     for aa in aa_fracs]

    # Convert to molar fractions (also, mol to mmol)
    if not mass:
        [aa_fracs.update({aa: 1000 * (aa_fracs[aa] / g_per_mol[aa])}) 
         for aa in aa_fracs]
    return aa_fracs


if __name__ == '__main__':
    aa_fracs = get_aa_fractions(0.50, False)    # mmol/gDW amino acid
