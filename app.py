from flask import Flask, request, render_template, redirect, url_for, send_from_directory
import os, uuid
import tempfile
from analysis import MitoAnalysisSystem
from models import SequenceAlignment, MotifFinder

webapp = Flask(__name__)
webapp.config['UPLOAD_FOLDER'] = tempfile.gettempdir()

mito_system = None  # Global instance


@webapp.route('/')
def index():
    return render_template("index.html")


@webapp.route('/upload', methods=['POST'])
def upload():
    global mito_system
    fasta_file = request.files['fasta']
    if fasta_file:
        path = os.path.join(webapp.config['UPLOAD_FOLDER'], fasta_file.filename)
        fasta_file.save(path)
        mito_system = MitoAnalysisSystem(path)
        return redirect(url_for('menu'))  
    return "Upload failed."

@webapp.route('/menu')
def menu():
    if mito_system is None:
        return redirect(url_for('index'))
    return render_template("menu.html")



@webapp.route('/stats')
def stats():
    if mito_system is None:
        return redirect(url_for('index'))
    df = mito_system.data.copy()
    df["GC content (%)"] = df["seq"].apply(lambda x: (x.count("G") + x.count("C")) / len(x) * 100)
    return render_template("stats.html", rows=df.to_dict(orient='records'))


@webapp.route('/compare', methods=['GET'])
def compare_form():
    if mito_system is None:
        return redirect(url_for('index'))
    genomes = mito_system.genomes
    return render_template("compare_select.html", genomes=genomes)

@webapp.route('/compare/results', methods=['POST'])
def compare_results():
    if mito_system is None:
        return redirect(url_for('index'))

    id1 = request.form.get("id1")
    id2 = request.form.get("id2")

    if id1 == id2:
        return "Please choose two different genomes."

    g1 = next((g for g in mito_system.genomes if g.id == id1), None)
    g2 = next((g for g in mito_system.genomes if g.id == id2), None)

    if not g1 or not g2:
        return "Invalid genome IDs selected."

    aligner = SequenceAlignment(g1.seq, g2.seq)
    seq1, comp, seq2 = aligner.align_sequences()

    # Block-wise formatting 
    block_size = 60
    rows = []
    for i in range(0, len(seq1), block_size):
        block1 = seq1[i:i+block_size]
        comp_block = comp[i:i+block_size]
        block2 = seq2[i:i+block_size]
        id_width = max(len(g1.id), len(g2.id))
        rows.append((
            f"{g1.id.ljust(id_width)}: {block1}",
            " " * (id_width + 2) + comp_block,
            f"{g2.id.ljust(id_width)}: {block2}"
        ))

    return render_template("compare.html", g1=g1, g2=g2, rows=rows)


@webapp.route('/motifs', methods=['GET'])
def motifs_form():
    if mito_system is None:
        return redirect(url_for('index'))
    return render_template("motif_form.html")


@webapp.route('/motifs/results', methods=['POST'])
def motifs_results():
    if mito_system is None:
        return redirect(url_for('index'))

    motif_input = request.form.get("motifs") or ""
    motifs = [m.strip().upper() for m in motif_input.splitlines() if m.strip()]

    if not motifs:
        return "No motifs entered."

    counts = []
    for motif in motifs:
        finder = MotifFinder(motif)
        row = [finder.count_occurrences(genome.seq) for genome in mito_system.genomes]
        counts.append((motif, row))

    # Image path
    img_filename = f"motif_plot_{uuid.uuid4().hex}.png"
    img_path = os.path.join(webapp.config['UPLOAD_FOLDER'], img_filename)

    # Plot
    import matplotlib.pyplot as plt
    import numpy as np

    if len(motifs) == 1:
        # === BAR PLOT ===
        motif = motifs[0]
        values = counts[0][1]  # row of counts
        genome_ids = [genome.id for genome in mito_system.genomes]

        plt.figure(figsize=(12, 5))
        bars = plt.bar(genome_ids, values, color='seagreen')
        plt.title(f"Occurrences of Motif '{motif}' Across Genomes")
        plt.xlabel("Genome ID")
        plt.ylabel("Count")
        plt.xticks(rotation=90, fontsize=6)
        plt.tight_layout()
        plt.savefig(img_path)
        plt.close()
    else:
        # Heatmap
        mito_system.motif_conservation_heatmap(motifs, img_path)

    return render_template("motif_results.html", motifs=counts, genomes=mito_system.genomes, heatmap_img=img_filename)

@webapp.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(webapp.config['UPLOAD_FOLDER'], filename)


@webapp.route('/heatmaps/<filename>')
def heatmap_file(filename):
    return send_from_directory(webapp.config['UPLOAD_FOLDER'], filename)


@webapp.route('/similarity', methods=['GET'])
def similarity_form():
    if mito_system is None:
        return redirect(url_for('index'))
    genomes = mito_system.genomes
    return render_template("similarity_select.html", genomes=genomes)

@webapp.route('/similarity/results', methods=['POST'])
def similarity_results():
    if mito_system is None:
        return redirect(url_for('index'))
    ref_id = request.form.get("ref_id")
    # find the genome index for this id
    ref_idx = next((i for i, g in enumerate(mito_system.genomes) if g.id == ref_id), 0)
    reference = mito_system.genomes[ref_idx]
    results = []
    for i, target in enumerate(mito_system.genomes):
        if i == ref_idx:
            continue
        aligner = SequenceAlignment(reference.seq, target.seq)
        _, comp, _ = aligner.align_sequences()
        score = aligner.get_alignment_scores()
        matches = comp.count("*")
        similarity = (matches / len(comp)) * 100 if comp else 0
        results.append({"id": target.id, "score": score, "similarity": similarity})
    return render_template("similarity.html", reference=reference, results=results)


webapp.run(debug=True)
