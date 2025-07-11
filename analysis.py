import matplotlib.pyplot as plt
import numpy as np
from parser import FastaParser
from models import MitochondrialDNA, MotifFinder, SequenceAlignment

class MitoAnalysisSystem:
    def __init__(self, fasta_file):
        parser = FastaParser()
        self.data = parser.parse(fasta_file, 'fasta')
        self.genomes = [
            MitochondrialDNA(seq=row['seq'], ID=row['id'], description=row['description'])
            for _, row in self.data.iterrows()
        ]

    def compare_two_species(self, idx1=0, idx2=1):
        g1 = self.genomes[idx1]
        g2 = self.genomes[idx2]
        aligner = SequenceAlignment(g1.seq, g2.seq)
        aligned_seq1, comparison_line, aligned_seq2 = aligner.align_sequences()
        self.visualize_differences_bar(aligned_seq1, aligned_seq2, g1.id, g2.id)
        return aligned_seq1, comparison_line, aligned_seq2

    def visualize_differences_bar(self, seq1, seq2, label1, label2):
        matches = mismatches = gaps = 0
        block_size = 60
        length = len(seq1)
        print(f"\nComparison between {label1} and {label2}:\n")
        for i in range(0, length, block_size):
            block1 = seq1[i:i+block_size]
            block2 = seq2[i:i+block_size]
            comp_line = ""
            for a, b in zip(block1, block2):
                if a == b and a != "-":
                    comp_line += "*"
                    matches += 1
                elif a == '-' or b == '-':
                    comp_line += " "
                    gaps += 1
                else:
                    comp_line += "|"
                    mismatches += 1
            label_width = max(len(label1), len(label2))
            print(f"{label1.ljust(label_width)}: {block1}")
            print(f"{' '.ljust(label_width)}  {comp_line}")
            print(f"{label2.ljust(label_width)}: {block2}\n")
        print("Summary:")
        print(f"  Matches   : {matches}")
        print(f"  Mismatches: {mismatches}")
        print(f"  Gaps      : {gaps}")
        print(f"  Total     : {length} positions")

    def motif_conservation_heatmap(self, motifs, output_path, threshold_score=4):
        conservation_matrix = []
        for motif in motifs:
            row = []
            for genome in self.genomes:
                aligner = SequenceAlignment(motif, genome.seq)
                score = aligner.get_alignment_scores(algo="local")
                row.append(1 if score >= threshold_score else 0)
            conservation_matrix.append(row)
        conservation_matrix = np.array(conservation_matrix)
        # Summary print
        for i, motif in enumerate(motifs):
            present_in = [self.genomes[j].id for j in range(len(self.genomes)) if conservation_matrix[i][j] == 1]
            print(f"Motif '{motif}' conserved in {len(present_in)} genome(s): {present_in}")
        # Plot
        fig, ax = plt.subplots(figsize=(14, 5))
        cax = ax.imshow(conservation_matrix, cmap="Greens", aspect='auto')
        ax.set_xticks(np.arange(len(self.genomes)))
        ax.set_yticks(np.arange(len(motifs)))
        ax.set_xticklabels([g.id for g in self.genomes], rotation=90, ha='right', fontsize=6)
        ax.set_yticklabels(motifs)
        plt.title("Conserved Motifs Across Species")
        plt.xlabel("Genomes")
        plt.ylabel("Motifs")
        plt.colorbar(cax, label='1 = conserved')
        plt.tight_layout()
        plt.savefig(output_path)
        # plt.show()  # Uncomment this if you want to see the plot interactively
        plt.close()

    def similarity_to_reference(self, ref_idx=0):
        reference = self.genomes[ref_idx]
        reference_id = reference.id
        print(f"\nReference genome: {reference_id}\n{'-'*50}")
        results = []
        for target in self.genomes[1:]:
            aligner = SequenceAlignment(reference.seq, target.seq)
            _, comp, _ = aligner.align_sequences()
            score = aligner.get_alignment_scores()
            matches = comp.count("*")
            total = len(comp)
            similarity = (matches / total) * 100 if total > 0 else 0
            results.append({
                "id": target.id,
                "score": score,
                "similarity": similarity
            })
            print(f"{target.id} | Score: {score} | Similarity: {similarity:.2f}%")
        # Bar chart
        plt.figure(figsize=(12, 5))
        ids = [r["id"] for r in results]
        similarities = [r["similarity"] for r in results]
        bars = plt.bar(ids, similarities, color="teal")
        plt.axhline(y=100, color='gray', linestyle='--', linewidth=0.5)
        plt.title(f"Similarity to Reference Genome {reference_id}")
        plt.ylabel("Similarity (%)")
        plt.xlabel("Genome ID")
        plt.xticks(rotation=90, ha='right', fontsize=6)
        plt.tight_layout()
        plt.show()
        
# # Usage example - Uncomment 
# if __name__ == "__main__":
#     mito_system = MitoAnalysisSystem('synthetic_mtDNA_dataset.fasta')
#     
#     # 1. Compare two species and show bar visualization
#     mito_system.compare_two_species(0, 1)
#     
#     # 2. Motif conservation across genomes
#     motifs = ["GATC", "TATA", "ATCG", "CGCG", "TTAA", "CTAG"]
#     mito_system.motif_conservation_heatmap(motifs)
#     
#     # 3. Similarity to reference genome (the first one)
#     mito_system.similarity_to_reference(ref_idx=0)
