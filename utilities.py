"""
Utility functions for DNA/RNA sequence analysis and PCR operations
"""

def transcribe(dna_seq):
    """Convert DNA to RNA (T -> U)"""
    return dna_seq.upper().replace('T', 'U')

def reverse_complement(seq):
    """Return reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join([complement.get(base, base) for base in seq.upper()][::-1])

def translate(rna_seq, start_codon=True):
    """Translate RNA to protein"""
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    rna_seq = rna_seq.upper()
    protein = []
    
    for i in range(0, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i+3]
        amino_acid = codon_table.get(codon, 'X')
        if amino_acid == '*':
            break
        protein.append(amino_acid)
    
    return ''.join(protein)

def gc_content(seq):
    """Calculate GC content percentage"""
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0

def find_orfs(dna_seq, min_length=100):
    """Find all Open Reading Frames (ORFs) in DNA sequence"""
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    
    dna_seq = dna_seq.upper()
    
    # Check all three reading frames
    for frame in range(3):
        i = frame
        while i < len(dna_seq) - 2:
            codon = dna_seq[i:i+3]
            if codon == start_codon:
                # Found start codon, look for stop
                for j in range(i+3, len(dna_seq) - 2, 3):
                    stop_codon = dna_seq[j:j+3]
                    if stop_codon in stop_codons:
                        orf_seq = dna_seq[i:j+3]
                        if len(orf_seq) >= min_length:
                            orfs.append({
                                'frame': frame + 1,
                                'start': i,
                                'end': j + 3,
                                'length': len(orf_seq),
                                'sequence': orf_seq
                            })
                        break
            i += 3
    
    return orfs

def calculate_tm_basic(primer):
    """Calculate melting temperature using basic formula (Wallace rule)"""
    primer = primer.upper()
    return 4 * (primer.count('G') + primer.count('C')) + 2 * (primer.count('A') + primer.count('T'))

def calculate_tm_advanced(primer):
    """Calculate melting temperature using nearest-neighbor method (simplified)"""
    primer = primer.upper()
    gc_count = primer.count('G') + primer.count('C')
    at_count = primer.count('A') + primer.count('T')
    
    # Simplified nearest-neighbor approximation
    if len(primer) < 14:
        return calculate_tm_basic(primer)
    else:
        return 64.9 + 41 * (gc_count - 16.4) / len(primer)

def validate_primer_pair(forward, reverse, template):
    """Validate primer pair for PCR"""
    template = template.upper()
    forward = forward.upper()
    reverse = reverse.upper()
    
    results = {
        'forward_found': forward in template,
        'reverse_found': reverse_complement(reverse) in template,
        'forward_tm': calculate_tm_advanced(forward),
        'reverse_tm': calculate_tm_advanced(reverse),
        'tm_difference': 0,
        'forward_gc': gc_content(forward),
        'reverse_gc': gc_content(reverse),
        'compatible': False
    }
    
    results['tm_difference'] = abs(results['forward_tm'] - results['reverse_tm'])
    
    # Check compatibility (Tm difference < 5°C is good)
    if results['tm_difference'] < 5 and results['forward_found'] and results['reverse_found']:
        results['compatible'] = True
    
    return results

def in_silico_pcr(forward, reverse, template):
    """Perform in silico PCR"""
    template = template.upper()
    forward = forward.upper()
    reverse_rc = reverse_complement(reverse.upper())
    
    # Find forward primer position
    forward_pos = template.find(forward)
    if forward_pos == -1:
        return None
    
    # Find reverse primer position (after forward)
    reverse_pos = template.find(reverse_rc, forward_pos)
    if reverse_pos == -1:
        return None
    
    # Extract amplicon
    amplicon = template[forward_pos:reverse_pos + len(reverse_rc)]
    
    return {
        'amplicon': amplicon,
        'length': len(amplicon),
        'forward_position': forward_pos,
        'reverse_position': reverse_pos,
        'gc_content': gc_content(amplicon)
    }

def codon_usage(dna_seq):
    """Calculate codon usage frequency"""
    dna_seq = dna_seq.upper()
    codon_count = {}
    
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        if len(codon) == 3:
            codon_count[codon] = codon_count.get(codon, 0) + 1
    
    return codon_count

def validate_sequence(seq, seq_type='DNA'):
    """Validate if sequence contains only valid nucleotides"""
    seq = seq.upper()
    if seq_type == 'DNA':
        valid_chars = set('ATGC')
    elif seq_type == 'RNA':
        valid_chars = set('AUGC')
    else:
        return False
    
    return all(char in valid_chars for char in seq)