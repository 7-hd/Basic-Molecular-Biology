"""
Digital Bioscience Suite (DBS)
Main application entry point
"""

import sys
from structure import SequenceManager, PrimerPair
from utilities import (
    transcribe, reverse_complement, translate, gc_content,
    find_orfs, calculate_tm_basic, calculate_tm_advanced,
    validate_primer_pair, in_silico_pcr, codon_usage,
    validate_sequence
)

# Global sequence manager
seq_manager = SequenceManager()

def print_header():
    """Print application header"""
    print("\n" + "="*60)
    print(" " * 15 + "Digital Bioscience Suite (DBS)")
    print("="*60 + "\n")

def main_menu():
    """Display main menu"""
    print("\n--- MAIN MENU ---")
    print("1. DNA/RNA Sequence Analysis")
    print("2. PCR & Primer Tools")
    print("3. Sequence Manager")
    print("4. Exit")
    return input("\nSelect option: ").strip()

def sequence_analysis_menu():
    """DNA/RNA sequence analysis submenu"""
    while True:
        print("\n--- SEQUENCE ANALYSIS ---")
        print("1. Transcribe DNA to RNA")
        print("2. Translate RNA to Protein")
        print("3. Reverse Complement")
        print("4. Calculate GC Content")
        print("5. Find ORFs (Open Reading Frames)")
        print("6. Codon Usage Analysis")
        print("7. Back to Main Menu")
        
        choice = input("\nSelect option: ").strip()
        
        if choice == '1':
            dna = input("Enter DNA sequence: ").strip()
            if validate_sequence(dna, 'DNA'):
                rna = transcribe(dna)
                print(f"\nRNA: {rna}")
            else:
                print("Error: Invalid DNA sequence!")
        
        elif choice == '2':
            rna = input("Enter RNA sequence: ").strip()
            if validate_sequence(rna, 'RNA'):
                protein = translate(rna)
                print(f"\nProtein: {protein}")
                print(f"Length: {len(protein)} amino acids")
            else:
                print("Error: Invalid RNA sequence!")
        
        elif choice == '3':
            dna = input("Enter DNA sequence: ").strip()
            if validate_sequence(dna, 'DNA'):
                rev_comp = reverse_complement(dna)
                print(f"\nReverse Complement: {rev_comp}")
            else:
                print("Error: Invalid DNA sequence!")
        
        elif choice == '4':
            seq = input("Enter DNA/RNA sequence: ").strip()
            gc = gc_content(seq)
            print(f"\nGC Content: {gc:.2f}%")
        
        elif choice == '5':
            dna = input("Enter DNA sequence: ").strip()
            min_len = input("Minimum ORF length (default 100): ").strip()
            min_len = int(min_len) if min_len else 100
            
            if validate_sequence(dna, 'DNA'):
                orfs = find_orfs(dna, min_len)
                print(f"\nFound {len(orfs)} ORF(s):")
                for idx, orf in enumerate(orfs, 1):
                    print(f"\nORF {idx}:")
                    print(f"  Frame: {orf['frame']}")
                    print(f"  Position: {orf['start']}-{orf['end']}")
                    print(f"  Length: {orf['length']} bp")
                    print(f"  Sequence: {orf['sequence'][:50]}..." if len(orf['sequence']) > 50 else f"  Sequence: {orf['sequence']}")
            else:
                print("Error: Invalid DNA sequence!")
        
        elif choice == '6':
            dna = input("Enter DNA sequence: ").strip()
            if validate_sequence(dna, 'DNA'):
                usage = codon_usage(dna)
                print("\nCodon Usage:")
                for codon, count in sorted(usage.items(), key=lambda x: x[1], reverse=True):
                    print(f"  {codon}: {count}")
            else:
                print("Error: Invalid DNA sequence!")
        
        elif choice == '7':
            break
        
        else:
            print("Invalid option!")

def pcr_primer_menu():
    """PCR and primer tools submenu"""
    while True:
        print("\n--- PCR & PRIMER TOOLS ---")
        print("1. Calculate Primer Tm (Melting Temperature)")
        print("2. Validate Primer Pair")
        print("3. In Silico PCR")
        print("4. Design Primer Pair (Manual)")
        print("5. Back to Main Menu")
        
        choice = input("\nSelect option: ").strip()
        
        if choice == '1':
            primer = input("Enter primer sequence: ").strip()
            if validate_sequence(primer, 'DNA'):
                tm_basic = calculate_tm_basic(primer)
                tm_advanced = calculate_tm_advanced(primer)
                gc = gc_content(primer)
                print(f"\nPrimer: {primer}")
                print(f"Length: {len(primer)} bp")
                print(f"GC Content: {gc:.2f}%")
                print(f"Tm (Wallace): {tm_basic:.1f}°C")
                print(f"Tm (Advanced): {tm_advanced:.1f}°C")
            else:
                print("Error: Invalid primer sequence!")
        
        elif choice == '2':
            forward = input("Enter forward primer: ").strip()
            reverse = input("Enter reverse primer: ").strip()
            template = input("Enter template DNA: ").strip()
            
            if all(validate_sequence(s, 'DNA') for s in [forward, reverse, template]):
                results = validate_primer_pair(forward, reverse, template)
                print("\n--- Primer Pair Validation ---")
                print(f"Forward primer found: {results['forward_found']}")
                print(f"Reverse primer found: {results['reverse_found']}")
                print(f"Forward Tm: {results['forward_tm']:.1f}°C")
                print(f"Reverse Tm: {results['reverse_tm']:.1f}°C")
                print(f"Tm Difference: {results['tm_difference']:.1f}°C")
                print(f"Forward GC: {results['forward_gc']:.1f}%")
                print(f"Reverse GC: {results['reverse_gc']:.1f}%")
                print(f"Compatible: {'YES' if results['compatible'] else 'NO'}")
            else:
                print("Error: Invalid sequences!")
        
        elif choice == '3':
            forward = input("Enter forward primer: ").strip()
            reverse = input("Enter reverse primer: ").strip()
            template = input("Enter template DNA: ").strip()
            
            if all(validate_sequence(s, 'DNA') for s in [forward, reverse, template]):
                result = in_silico_pcr(forward, reverse, template)
                if result:
                    print("\n--- PCR Product ---")
                    print(f"Amplicon length: {result['length']} bp")
                    print(f"GC content: {result['gc_content']:.1f}%")
                    print(f"Forward position: {result['forward_position']}")
                    print(f"Reverse position: {result['reverse_position']}")
                    print(f"\nAmplicon sequence:")
                    print(result['amplicon'][:100] + "..." if len(result['amplicon']) > 100 else result['amplicon'])
                else:
                    print("\nNo PCR product found! Check primer binding sites.")
            else:
                print("Error: Invalid sequences!")
        
        elif choice == '4':
            name = input("Enter primer pair name: ").strip()
            forward = input("Enter forward primer sequence: ").strip()
            reverse = input("Enter reverse primer sequence: ").strip()
            
            if validate_sequence(forward, 'DNA') and validate_sequence(reverse, 'DNA'):
                pair = PrimerPair(name, forward, reverse)
                info = pair.get_info()
                print(f"\n--- Primer Pair: {info['name']} ---")
                print(f"Forward: {info['forward_seq']}")
                print(f"  Tm: {info['forward_tm']:.1f}°C, GC: {info['forward_gc']:.1f}%")
                print(f"Reverse: {info['reverse_seq']}")
                print(f"  Tm: {info['reverse_tm']:.1f}°C, GC: {info['reverse_gc']:.1f}%")
            else:
                print("Error: Invalid primer sequences!")
        
        elif choice == '5':
            break
        
        else:
            print("Invalid option!")

def sequence_manager_menu():
    """Sequence manager submenu"""
    while True:
        print("\n--- SEQUENCE MANAGER ---")
        print("1. Add Sequence")
        print("2. View Sequence")
        print("3. List All Sequences")
        print("4. Delete Sequence")
        print("5. Back to Main Menu")
        
        choice = input("\nSelect option: ").strip()
        
        if choice == '1':
            name = input("Enter sequence name: ").strip()
            sequence = input("Enter sequence: ").strip()
            seq_type = input("Type (DNA/RNA) [DNA]: ").strip().upper() or 'DNA'
            
            if validate_sequence(sequence, seq_type):
                seq_manager.add_sequence(name, sequence, seq_type)
                print(f"\nSequence '{name}' added successfully!")
            else:
                print("Error: Invalid sequence!")
        
        elif choice == '2':
            name = input("Enter sequence name: ").strip()
            seq = seq_manager.get_sequence(name)
            if seq:
                info = seq.get_info()
                print(f"\n--- {info['name']} ---")
                print(f"Type: {info['type']}")
                print(f"Length: {info['length']} bp/nt")
                print(f"GC Content: {info['gc_content']:.2f}%")
                print(f"Sequence: {info['sequence'][:100]}..." if len(info['sequence']) > 100 else f"Sequence: {info['sequence']}")
            else:
                print(f"Sequence '{name}' not found!")
        
        elif choice == '3':
            sequences = seq_manager.list_sequences()
            if sequences:
                print(f"\nStored Sequences ({seq_manager.get_count()}):")
                for name in sequences:
                    seq = seq_manager.get_sequence(name)
                    print(f"  - {seq}")
            else:
                print("\nNo sequences stored.")
        
        elif choice == '4':
            name = input("Enter sequence name to delete: ").strip()
            if seq_manager.delete_sequence(name):
                print(f"\nSequence '{name}' deleted.")
            else:
                print(f"Sequence '{name}' not found!")
        
        elif choice == '5':
            break
        
        else:
            print("Invalid option!")

def main():
    """Main application loop"""
    print_header()
    
    while True:
        choice = main_menu()
        
        if choice == '1':
            sequence_analysis_menu()
        elif choice == '2':
            pcr_primer_menu()
        elif choice == '3':
            sequence_manager_menu()
        elif choice == '4':
            print("\nExiting DBS. Goodbye!\n")
            sys.exit(0)
        else:
            print("Invalid option! Please try again.")

if __name__ == "__main__":
    main()