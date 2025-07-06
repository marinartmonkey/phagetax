# phagetax.py
# Main script for the PhageTax tool, integrates all functionality for batch annotation, feature extraction, matrix generation and taxonomic prediction.

import os
import glob
import subprocess
import pandas as pd
import pickle
from tqdm import tqdm
import sys

# -------------------------
# Global constants and paths
# -------------------------
GENOMES_DIR = "entrada_genomas"
PHAROKKA_OUTPUT_DIR = "pharokka_out"
MODELS_DIR = "models"
PHROGS_INFO_DIR = "phrogs_info"

BINARY_MATRIX_PATH = "df_phrogs_binaria.tsv"
SUMMARY_OUTPUT_PATH = "matriz_phrogs_tax_gc_simple.tsv"
FINAL_MATRIX_PATH = "matriz_phrogs_tax_gc_final.tsv"
PREDICTIONS_OUTPUT_PATH = "predicciones_nuevos.tsv"
PHROGS_REPORT_PATH = "phrogs_presentes_por_genoma.tsv"

TAX_LEVELS = ["Kingdom", "Phylum", "Class", "Family", "Genus", "Species"]
THREADS = 4

# -------------------------
# Function 1: Batch annotation with Pharokka
# -------------------------
def annotate_pending_genomes():
    fasta_files = glob.glob(os.path.join(GENOMES_DIR, "*.fasta"))
    print(f"[INFO] Found {len(fasta_files)} .fasta files in {GENOMES_DIR}")

    os.makedirs(PHAROKKA_OUTPUT_DIR, exist_ok=True)

    for fasta in fasta_files:
        genome_name = os.path.splitext(os.path.basename(fasta))[0]
        genome_output_dir = os.path.join(PHAROKKA_OUTPUT_DIR, genome_name)

        if os.path.exists(genome_output_dir):
            print(f"[INFO] {genome_name} already annotated. Skipping...")
            continue

        print(f"[INFO] Annotating {genome_name}...")
        cmd = f"pharokka.py -i {fasta} -o {genome_output_dir} -t {THREADS} -f"
        subprocess.run(cmd, shell=True, check=True)

    print(f"[INFO] Annotation complete. Results saved in {PHAROKKA_OUTPUT_DIR}")

# -------------------------
# Function 2: Generate binary PHROG matrix
# -------------------------
def generate_phrog_matrix():
    genomes = sorted([
        d for d in os.listdir(PHAROKKA_OUTPUT_DIR)
        if os.path.isdir(os.path.join(PHAROKKA_OUTPUT_DIR, d)) and not d.startswith(".")
    ])

    phrog_dict = {}
    for genome_id in tqdm(genomes, desc="Processing annotated genomes"):
        phrog_file = os.path.join(PHAROKKA_OUTPUT_DIR, genome_id, "pharokka_cds_final_merged_output.tsv")
        try:
            df_phrog = pd.read_csv(phrog_file, sep="\t")
            if "phrog" not in df_phrog.columns:
                raise ValueError("Missing 'phrog' column.")

            phrogs = df_phrog["phrog"].dropna().unique()
            phrogs = [str(p).strip().replace("phrog_", "") for p in phrogs if isinstance(p, str) and p != "No_PHROG"]
            phrog_dict[genome_id] = phrogs

        except Exception as e:
            print(f"[ERROR] Reading PHROGs for {genome_id}: {e}")
            phrog_dict[genome_id] = []

    unique_phrogs = sorted(set(phrog for lst in phrog_dict.values() for phrog in lst))

    df_bin = pd.DataFrame(0, index=genomes, columns=[f"PHROG_{p}" for p in unique_phrogs])
    for genome_id, lst in phrog_dict.items():
        for phrog in lst:
            col = f"PHROG_{phrog}"
            if col in df_bin.columns:
                df_bin.at[genome_id, col] = 1

    df_bin.to_csv(BINARY_MATRIX_PATH, sep="\t")
    print(f"[INFO] Binary PHROG matrix saved to {BINARY_MATRIX_PATH}")

# -------------------------
# Function 3: Genome summary extraction
# -------------------------
def generate_genome_summary():
    genomes = sorted([
        d for d in os.listdir(PHAROKKA_OUTPUT_DIR)
        if os.path.isdir(os.path.join(PHAROKKA_OUTPUT_DIR, d)) and not d.startswith(".")
    ])

    info_data = []
    for genome in tqdm(genomes, desc="Extracting genome information"):
        resumen_path = os.path.join(PHAROKKA_OUTPUT_DIR, genome, "pharokka_length_gc_cds_density.tsv")
        phrog_path = os.path.join(PHAROKKA_OUTPUT_DIR, genome, "pharokka_cds_final_merged_output.tsv")
        try:
            df_resumen = pd.read_csv(resumen_path, sep="\t")
            length = df_resumen["length"].iloc[0]
            gc_perc = df_resumen["gc_perc"].iloc[0]
            num_genes = pd.read_csv(phrog_path, sep="\t").shape[0] if os.path.exists(phrog_path) else 0
            info_data.append({"genome_id": genome, "genome_length": length, "gc_percent": gc_perc, "num_genes": num_genes})

        except Exception as e:
            print(f"[ERROR] Processing {genome}: {e}")

    df_info = pd.DataFrame(info_data).set_index("genome_id")
    df_info.to_csv(SUMMARY_OUTPUT_PATH, sep="\t")
    print(f"[INFO] Genome summary saved to {SUMMARY_OUTPUT_PATH}")

# -------------------------
# Function 4: Final matrix construction with PHROG features
# -------------------------
def generate_final_matrix():
    df_bin = pd.read_csv(BINARY_MATRIX_PATH, sep="\t", index_col=0)
    df_pharokka = pd.read_csv(SUMMARY_OUTPUT_PATH, sep="\t", index_col=0)
    phrog_cols = [col for col in df_bin.columns if col.startswith("PHROG_")]

    for col in phrog_cols:
        if col in df_pharokka.columns:
            df_bin[col] = df_bin[col] | df_pharokka[col]

    df_tax_gc = df_pharokka.drop(columns=[col for col in df_pharokka.columns if col.startswith("PHROG_")], errors="ignore")
    df_final = df_tax_gc.join(df_bin, how="outer").fillna(0)

    phrogs_especificos = {nivel: pd.read_csv(f"{PHROGS_INFO_DIR}/phrogs_especificos_por_{nivel.lower()}.tsv", sep="\t") for nivel in TAX_LEVELS}
    phrogs_caracteristicos = {nivel: pd.read_csv(f"{PHROGS_INFO_DIR}/phrogs_caracteristicos_por_{nivel.lower()}.tsv", sep="\t") for nivel in TAX_LEVELS}

    phrog_cols_final = [col for col in df_final.columns if col.startswith("PHROG_")]
    df_extra = pd.DataFrame(index=df_final.index)

    for nivel in TAX_LEVELS:
        total_esp = [p for col in phrogs_especificos[nivel].columns for p in phrogs_especificos[nivel][col].dropna()]
        total_carac = [p for col in phrogs_caracteristicos[nivel].columns for p in phrogs_caracteristicos[nivel][col].dropna()]

        df_extra[f"especificos_{nivel.lower()}"] = [3 * len(set(row[phrog_cols_final][row[phrog_cols_final] == 1].index) & set(total_esp)) for _, row in df_final.iterrows()]
        df_extra[f"caracteristicos_{nivel.lower()}"] = [len(set(row[phrog_cols_final][row[phrog_cols_final] == 1].index) & set(total_carac)) for _, row in df_final.iterrows()]

    df_final = pd.concat([df_final, df_extra], axis=1)
    df_final.to_csv(FINAL_MATRIX_PATH, sep="\t")
    print(f"[INFO] Final matrix saved to {FINAL_MATRIX_PATH}")

# -------------------------
# Function 5: Taxonomic prediction
# -------------------------
def predict_from_matrix():
    df = pd.read_csv(FINAL_MATRIX_PATH, sep="\t", index_col=0)
    resultados = pd.DataFrame(index=df.index)

    for nivel in TAX_LEVELS:
        print(f"[INFO] Predicting taxonomic level: {nivel}")
        with open(os.path.join(MODELS_DIR, f"modelo_{nivel.lower()}.pkl"), "rb") as f:
            modelo = pickle.load(f)

        cols_usadas = modelo.feature_names_in_.tolist()
        cols_presentes = [col for col in cols_usadas if col in df.columns]
        cols_faltantes = [col for col in cols_usadas if col not in df.columns]

        df_completo = pd.concat([df[cols_presentes], pd.DataFrame(0, index=df.index, columns=cols_faltantes)], axis=1)[cols_usadas]
        predicciones = modelo.predict(df_completo)
        probs = modelo.predict_proba(df_completo)

        resultados[nivel] = predicciones
        resultados[f"Confianza_{nivel}"] = [probs[i, list(modelo.classes_).index(pred)] for i, pred in enumerate(predicciones)]

    resultados.to_csv(PREDICTIONS_OUTPUT_PATH, sep="\t")
    print(f"[INFO] Predictions saved to {PREDICTIONS_OUTPUT_PATH}")

# -------------------------
# Function 6: PHROGs presence report
# -------------------------
def generate_phrog_presence_report():
    df = pd.read_csv(FINAL_MATRIX_PATH, sep="\t", index_col=0)
    phrog_cols = [col for col in df.columns if col.startswith("PHROG_")]

    phrogs_especificos = {nivel: pd.read_csv(f"{PHROGS_INFO_DIR}/phrogs_especificos_por_{nivel.lower()}.tsv", sep="\t") for nivel in TAX_LEVELS}
    phrogs_caracteristicos = {nivel: pd.read_csv(f"{PHROGS_INFO_DIR}/phrogs_caracteristicos_por_{nivel.lower()}.tsv", sep="\t") for nivel in TAX_LEVELS}

    results = []
    for idx, row in df.iterrows():
        presentes = row[phrog_cols][row[phrog_cols] == 1].index.tolist()
        info = {"genome_id": idx}

        for nivel in TAX_LEVELS:
            total_esp = [p for col in phrogs_especificos[nivel].columns for p in phrogs_especificos[nivel][col].dropna()]
            total_carac = [p for col in phrogs_caracteristicos[nivel].columns for p in phrogs_caracteristicos[nivel][col].dropna()]
            info[f"PHROGs_especificos_{nivel}"] = ";".join(sorted(set(presentes) & set(total_esp))) if set(presentes) & set(total_esp) else "None"
            info[f"PHROGs_caracteristicos_{nivel}"] = ";".join(sorted(set(presentes) & set(total_carac))) if set(presentes) & set(total_carac) else "None"

        results.append(info)

    pd.DataFrame(results).set_index("genome_id").to_csv(PHROGS_REPORT_PATH, sep="\t")
    print(f"[INFO] PHROGs presence report saved to {PHROGS_REPORT_PATH}")

# -------------------------
# Main logic with command-line arguments
# -------------------------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python phagetax.py [annotate_batch | generate_matrix | genome_summary | final_matrix | predict | phrog_report]")
        sys.exit(1)

    cmd = sys.argv[1]

    if cmd == "annotate_batch":
        annotate_pending_genomes()
    elif cmd == "generate_matrix":
        generate_phrog_matrix()
    elif cmd == "genome_summary":
        generate_genome_summary()
    elif cmd == "final_matrix":
        generate_final_matrix()
    elif cmd == "predict":
        predict_from_matrix()
    elif cmd == "phrog_report":
        generate_phrog_presence_report()
    else:
        print(f"[ERROR] Unknown command: {cmd}")
        sys.exit(1)
