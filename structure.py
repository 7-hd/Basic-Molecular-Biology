"""
Data structures and sequence storage for DBS
"""

class Sequence:
    """Class to store and manage biological sequences"""
    
    def __init__(self, name, sequence, seq_type='DNA'):
        self.name = name
        self.sequence = sequence.upper()
        self.seq_type = seq_type
    
    def __str__(self):
        return f"{self.name} ({self.seq_type}): {len(self.sequence)} bp/nt"
    
    def get_info(self):
        """Return detailed information about the sequence"""
        from utilities import gc_content
        return {
            'name': self.name,
            'type': self.seq_type,
            'length': len(self.sequence),
            'gc_content': gc_content(self.sequence),
            'sequence': self.sequence
        }

class SequenceManager:
    """Manager class to handle multiple sequences"""
    
    def __init__(self):
        self.sequences = {}
    
    def add_sequence(self, name, sequence, seq_type='DNA'):
        """Add a new sequence to the manager"""
        self.sequences[name] = Sequence(name, sequence, seq_type)
        return True
    
    def get_sequence(self, name):
        """Retrieve a sequence by name"""
        return self.sequences.get(name)
    
    def list_sequences(self):
        """List all stored sequences"""
        return list(self.sequences.keys())
    
    def delete_sequence(self, name):
        """Delete a sequence from the manager"""
        if name in self.sequences:
            del self.sequences[name]
            return True
        return False
    
    def get_count(self):
        """Get total number of stored sequences"""
        return len(self.sequences)

class Primer:
    """Class to store primer information"""
    
    def __init__(self, name, sequence, direction='forward'):
        self.name = name
        self.sequence = sequence.upper()
        self.direction = direction
    
    def __str__(self):
        return f"{self.name} ({self.direction}): {self.sequence}"

class PrimerPair:
    """Class to store and manage primer pairs"""
    
    def __init__(self, name, forward_primer, reverse_primer):
        self.name = name
        self.forward = Primer(f"{name}_F", forward_primer, 'forward')
        self.reverse = Primer(f"{name}_R", reverse_primer, 'reverse')
    
    def get_info(self):
        """Return primer pair information"""
        from utilities import calculate_tm_advanced, gc_content
        return {
            'name': self.name,
            'forward_seq': self.forward.sequence,
            'reverse_seq': self.reverse.sequence,
            'forward_tm': calculate_tm_advanced(self.forward.sequence),
            'reverse_tm': calculate_tm_advanced(self.reverse.sequence),
            'forward_gc': gc_content(self.forward.sequence),
            'reverse_gc': gc_content(self.reverse.sequence)
        }