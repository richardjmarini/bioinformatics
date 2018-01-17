class Nucleotide:

    _table= {
        'A': {
            'name': 'Adeine',
            'type': 'purine'
        },
        'C': {
            'name': 'Cytosine',
            'type': 'pyridmidine'
        },
        'T': {
            'name': 'Thymine',
            'type': 'pyridmidine'
        },
        'G': {
            'name': 'Guanine',
            'type': 'purine'
        },
        'U': {
            'name': 'Uracil',
            'type': 'pyridmidine'
        }
    }

    def __init__(self, code):
        self.code= code
        properties= self._table.get(self.code.upper())
        self.__dict__.update(properties.items())

    def __str__(self):

        return self.code

