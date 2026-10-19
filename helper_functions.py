import math

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    #check sequences are of equal length
    if len(reference_sequence) != len(distant_sequence):
        raise ValueError("Aligned sequences must have equal lengths.")

    #convert sequences to uppercase for consistency
    reference_sequence = reference_sequence.upper()
    distant_sequence = distant_sequence.upper()

    valid_positions = 0
    differences = 0
    #compare the two sequences and count valid positions and differences
    for a, b in zip(reference_sequence, distant_sequence):  
        if a != '-' and b != '-':
            valid_positions += 1
            if a != b:
                differences += 1 

    #calculate the proportion of differences
    if valid_positions == 0:
        raise ValueError("No valid positions to compare.")
    p = differences / valid_positions

    #apply the Jukes-Cantor correction
    if p >= 0.75:
        raise ValueError("Proportion of differences is too high for Jukes-Cantor correction.")
    
    d = -3/4 * math.log(1 - 4/3 * p)
    return d * valid_positions


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    raise NotImplementedError()
