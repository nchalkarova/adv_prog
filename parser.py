import pandas as pd
from Bio import SeqIO


class FastaParser:
    def __init__(self):
        self._df = None

    def parse(self, file:str, format:str = 'fasta'):
        try: 
            with open(file) as f:
                pass
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {file}")

        data = {}
        output = pd.DataFrame()

        records = list(SeqIO.parse(file, format))
        if not records:
            raise ValueError(f"No valid records found in file: {file}")
        for record in records:
            for key, value in record.__dict__.items():
                if key not in data:
                    data[key] = []
                if hasattr(value, '_data'):
                    value = str(value)
                data[key].append(value)

        for key, values in data.items():
            if values and len(values) == len(records):
                output[key.strip('_')] = values

        if 'seq' in output.columns:
            output['length'] = output['seq'].str.len()
        self._df = output
        return output

    def get_dataframe(self):
        if self._df is None:
            raise ValueError("No data: parse() not yet called")
        return self._df

# Usage 
if __name__ == "__main__":
    parser = FastaParser()
    df = parser.parse('synthetic_mtDNA_dataset.fasta')
    print(df.head())
    parser.get_dataframe().to_csv('out.csv')
