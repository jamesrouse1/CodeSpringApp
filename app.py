from __future__ import annotations

import ast
import base64
import json
import os
import re
import shutil
import socket
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components


def find_codespringlab_root() -> Path:
    env_root = os.environ.get("CSL_CODESPRINGLAB_ROOT", "").strip()
    candidates = []
    if env_root:
        candidates.append(Path(env_root).expanduser())
    app_path = Path(__file__).resolve()
    candidates.extend([
        app_path.parents[1],
        Path.cwd(),
        Path.cwd().parent,
        Path("~/CodeSpringLab").expanduser(),
        Path("~/CSH/CodeSpringLab").expanduser(),
        Path("/grid/bsr/home/rouse/CodeSpringLab"),
        Path("/Users/rouse/CSH/CodeSpringLab"),
    ])
    for candidate in candidates:
        if (candidate / "scripts_DoNotTouch").is_dir():
            return candidate
    return app_path.parents[1]


REPO_ROOT = find_codespringlab_root()
SCRIPTS = REPO_ROOT / "scripts_DoNotTouch"
APP_HOME = Path(os.environ.get("CSL_WEB_HOME", "~/.codespringlab_web")).expanduser()
PROJECTS_PATH = APP_HOME / "projects.json"
JOBS_PATH = APP_HOME / "jobs.json"
CONFIGS_DIR = APP_HOME / "configs"
FASTQ_SUFFIXES = [".fastq.gz", ".fq.gz", ".fastq", ".fq"]

GENOME_RESOURCES = {
    "mouse": {
        "star_index": "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/GRCm39_M29_gencode_starindex",
        "kallisto_index": "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.transcripts.idx",
        "gtf": "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation.gtf",
        "strand_bed": "/grid/bsr/data/data/utama/genome/GRCm39_M29_gencode/gencode.vM29.annotation_forStrandDetect_geneID.bed",
    },
    "human": {
        "star_index": "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/hg38_p13_gencode_rel42_all_starindex",
        "kallisto_index": "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v45.transcripts.idx",
        "gtf": "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation.gtf",
        "strand_bed": "/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation_forStrandDetect_geneID.bed",
    },
}


def now_stamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def ensure_home() -> None:
    APP_HOME.mkdir(parents=True, exist_ok=True)
    CONFIGS_DIR.mkdir(parents=True, exist_ok=True)


def read_json(path: Path, default):
    ensure_home()
    if not path.exists():
        return default
    try:
        return json.loads(path.read_text())
    except Exception:
        return default


def write_json(path: Path, data) -> None:
    ensure_home()
    path.write_text(json.dumps(data, indent=2, sort_keys=True))


def analysis_slug(value: str) -> str:
    raw = str(value or "RNA-seq").lower()
    if "atac" in raw:
        return "atac_seq"
    if "chip" in raw:
        return "chip_seq"
    return "rna_seq"


def project_id(project: dict) -> str:
    return f"{analysis_slug(project.get('analysis_type', 'RNA-seq'))}/{clean_name(project.get('name', 'project'), 'project')}"


def project_config_path(project: dict) -> Path:
    return CONFIGS_DIR / analysis_slug(project.get("analysis_type", "RNA-seq")) / f"{clean_name(project.get('name', 'project'), 'project')}.json"


def normalize_project(project: dict) -> dict:
    project = dict(project or {})
    project["name"] = clean_name(project.get("name", "project"), "project")
    project["analysis_type"] = project.get("analysis_type", "RNA-seq")
    project["project_id"] = project_id(project)
    return project


def load_projects() -> Dict[str, dict]:
    registry = read_json(PROJECTS_PATH, {})
    projects: Dict[str, dict] = {}
    if isinstance(registry, dict):
        for key, value in registry.items():
            if isinstance(value, dict):
                project = normalize_project(value)
                projects[project_id(project)] = project

    ensure_home()
    if CONFIGS_DIR.exists():
        for config_path in sorted(CONFIGS_DIR.glob("*/*.json")):
            try:
                project = normalize_project(json.loads(config_path.read_text()))
            except Exception:
                continue
            projects[project_id(project)] = project

    for project in discover_codespringlab_projects():
        pid = project_id(project)
        if pid not in projects:
            projects[pid] = project
    return projects


def save_project(project: dict) -> dict:
    projects = load_projects()
    project = normalize_project(project)
    project["updated_at"] = now_stamp()
    if not project.get("created_at"):
        project["created_at"] = project["updated_at"]
    projects[project_id(project)] = project
    write_json(PROJECTS_PATH, projects)
    config_path = project_config_path(project)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(project, indent=2, sort_keys=True))
    return project


ANALYSIS_DIR_TO_LABEL = {
    "rna": "RNA-seq",
    "rnaseq": "RNA-seq",
    "rna_seq": "RNA-seq",
    "atac": "ATAC-seq",
    "atacseq": "ATAC-seq",
    "atac_seq": "ATAC-seq",
    "chip": "ChIP-seq",
    "chipseq": "ChIP-seq",
    "chip_seq": "ChIP-seq",
}

ANALYSIS_NOTEBOOK_DIR = {
    "rna": "bulkRNAseq",
    "rnaseq": "bulkRNAseq",
    "rna_seq": "bulkRNAseq",
    "atac": "bulkATACseq",
    "atacseq": "bulkATACseq",
    "atac_seq": "bulkATACseq",
    "chip": "bulkChIPseq",
    "chipseq": "bulkChIPseq",
    "chip_seq": "bulkChIPseq",
}


def read_python_config(path: Path) -> dict:
    values = {}
    try:
        tree = ast.parse(path.read_text())
    except Exception:
        return values
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        try:
            value = ast.literal_eval(node.value)
        except Exception:
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and not target.id.startswith("_"):
                values[target.id] = value
    return values


def legacy_analysis_key(path: Path, values: dict) -> str:
    raw = str(values.get("analysis_type", "") or path.parent.name or "rna").lower()
    if "atac" in raw:
        return "atac"
    if "chip" in raw:
        return "chip"
    return "rna"


def legacy_base_dir(analysis_key: str) -> Path:
    return REPO_ROOT / ANALYSIS_NOTEBOOK_DIR.get(analysis_key, "bulkRNAseq")


def resolve_legacy_path(value, analysis_key: str, prefer_folder: bool = True) -> str:
    value = str(value or "").strip()
    if not value:
        return ""
    path = Path(value).expanduser()
    if path.is_absolute():
        resolved = path
    else:
        resolved = (legacy_base_dir(analysis_key) / path).resolve()
    if prefer_folder and value.endswith("/"):
        return str(resolved)
    return str(resolved)


def legacy_results_root(values: dict, analysis_key: str, project_name: str) -> str:
    visualizer = str(values.get("visualizer_data_dir", "")).strip()
    if visualizer:
        data_path = Path(resolve_legacy_path(visualizer, analysis_key)).expanduser()
        if data_path.name == "data" and data_path.parent.name == project_name:
            return str(data_path.parent.parent)
    raw = values.get("results_directory", "../../csl_results/")
    return resolve_legacy_path(raw, analysis_key)


def legacy_project_from_values(values: dict, config_path: Path) -> Optional[dict]:
    analysis_key = legacy_analysis_key(config_path, values)
    analysis_type = ANALYSIS_DIR_TO_LABEL.get(analysis_key, "RNA-seq")
    project_name = str(values.get("project_name", "") or config_path.stem).strip()
    if not project_name:
        return None

    results_root = legacy_results_root(values, analysis_key, project_name)
    inpath_design = resolve_legacy_path(values.get("inpath_design", ""), analysis_key)
    read_dest = resolve_legacy_path(values.get("read_path_destination", ""), analysis_key)
    read_orig = resolve_legacy_path(values.get("read_path_original", ""), analysis_key)
    visualizer = resolve_legacy_path(values.get("visualizer_data_dir", ""), analysis_key)
    pairing = str(values.get("pairing", values.get("paired_end", ""))).strip().lower()
    paired = False if pairing in ["n", "no", "false", "single", "single-end", "se"] else True

    project = {
        "name": clean_name(project_name, "project"),
        "analysis_type": analysis_type,
        "workflow_mode": "Visualize existing results",
        "genome": str(values.get("genome", "mouse") or "mouse").strip().lower(),
        "paired_end": paired,
        "results_root": results_root,
        "fastq_dir": read_dest or read_orig,
        "design_matrix_path": inpath_design,
        "data_dir_override": visualizer,
        "metadata_columns": ["treatment"],
        "source_config": str(config_path),
        "source": "CodeSpringLab config",
    }
    project = normalize_project(project)
    project = apply_project_inference(project)
    project["source_config"] = str(config_path)
    project["source"] = "CodeSpringLab config"
    return project


def legacy_config_candidates() -> List[Path]:
    candidates = []
    search_roots = [
        SCRIPTS / "project_configs",
        REPO_ROOT / "project_configs",
    ]
    for root in search_roots:
        if root.is_dir():
            candidates.extend(sorted(root.glob("*/*.py")))
    active_config = SCRIPTS / "config.py"
    if active_config.exists():
        candidates.append(active_config)
    seen = set()
    unique = []
    for path in candidates:
        key = str(path.resolve())
        if key not in seen:
            seen.add(key)
            unique.append(path)
    return unique


def discover_codespringlab_projects() -> List[dict]:
    projects = []
    for config_path in legacy_config_candidates():
        values = read_python_config(config_path)
        if not values and config_path.name != "config.py":
            continue
        project = legacy_project_from_values(values, config_path)
        if project:
            projects.append(project)
    return projects


def load_jobs() -> List[dict]:
    return read_json(JOBS_PATH, [])


def save_job(job: dict) -> None:
    jobs = load_jobs()
    jobs.append(job)
    write_json(JOBS_PATH, jobs)


def project_jobs(project) -> List[dict]:
    if isinstance(project, dict):
        pid = project_id(project)
        name = project.get("name")
        analysis = project.get("analysis_type")
        return [
            j for j in load_jobs()
            if j.get("project_id") == pid
            or (j.get("project") == name and j.get("analysis_type", analysis) == analysis)
        ]
    return [j for j in load_jobs() if j.get("project") == project]


def clean_name(value: str, fallback: str = "sample") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", str(value).strip()).strip("_")
    return cleaned or fallback


def parse_metadata_columns(value: str, fallback: Optional[List[str]] = None) -> List[str]:
    fallback = fallback or ["treatment"]
    cols = [clean_name(x, "condition") for x in str(value or "").split(",") if x.strip()]
    cols = [c for c in cols if c not in ["sample", "filename", "include", "detected_status"]]
    deduped = []
    for col in cols:
        if col not in deduped:
            deduped.append(col)
    return deduped or fallback


def sync_design_editor_columns(df: pd.DataFrame, metadata_cols: List[str]) -> pd.DataFrame:
    if df.empty:
        df = pd.DataFrame(columns=["include", "sample", "filename", "detected_status"])
    df = df.copy()
    for col in ["include", "sample", "filename", "detected_status"]:
        if col not in df.columns:
            df[col] = True if col == "include" else ""
    for col in metadata_cols:
        if col not in df.columns:
            df[col] = ""
    ordered = ["include", "sample"] + metadata_cols + ["filename", "detected_status"]
    extra = [c for c in df.columns if c not in ordered]
    return df[ordered + extra]


def project_root(project: dict) -> Path:
    return Path(project["results_root"]).expanduser() / project["name"]


def data_dir(project: dict) -> Path:
    override = str(project.get("data_dir_override", "")).strip()
    if override:
        return Path(override).expanduser()
    return project_root(project) / "data"


def log_dir(project: dict) -> Path:
    return project_root(project) / "log"


def manifest_dir(project: dict) -> Path:
    return data_dir(project) / "manifest"


def design_matrix_path(project: dict) -> Path:
    override = str(project.get("design_matrix_path", "")).strip()
    if override:
        path = Path(override).expanduser()
        return path if path.name == "design_matrix.txt" else path / "design_matrix.txt"
    inpath_design = str(project.get("inpath_design", "")).strip()
    if inpath_design:
        path = Path(inpath_design).expanduser()
        return path if path.name == "design_matrix.txt" else path / "design_matrix.txt"
    return manifest_dir(project) / "design_matrix.txt"


def candidate_design_matrix_path(project: dict) -> Path:
    candidates = [
        design_matrix_path(project),
        data_dir(project) / "manifest" / "design_matrix.txt",
        data_dir(project) / "design_matrix" / "design_matrix.txt",
        data_dir(project) / "design_matrix.txt",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def infer_fastq_dir(project: dict) -> str:
    current = str(project.get("fastq_dir", "")).strip()
    if current and Path(current).expanduser().is_dir():
        return str(Path(current).expanduser())
    candidates = [
        data_dir(project) / "fastq",
        data_dir(project) / "cutadapt",
    ]
    for candidate in candidates:
        if candidate.is_dir():
            return str(candidate)
    return current


def infer_metadata_columns_from_design(project: dict) -> List[str]:
    design_path = candidate_design_matrix_path(project)
    if not design_path.exists():
        return project.get("metadata_columns", ["treatment"])
    try:
        design = pd.read_table(design_path, nrows=5)
    except Exception:
        return project.get("metadata_columns", ["treatment"])
    cols = [c for c in design.columns if c not in ["sample", "filename"]]
    return cols or project.get("metadata_columns", ["treatment"])


def apply_project_inference(project: dict) -> dict:
    project["fastq_dir"] = infer_fastq_dir(project)
    inferred_design = candidate_design_matrix_path(project)
    if inferred_design.exists():
        project["design_matrix_path"] = str(inferred_design)
    project["metadata_columns"] = infer_metadata_columns_from_design(project)
    return project


def count_files(folder: Path, patterns: Iterable[str]) -> int:
    if not folder.exists():
        return 0
    total = 0
    for pattern in patterns:
        total += len(list(folder.rglob(pattern)))
    return total


def project_step_status(project: dict) -> pd.DataFrame:
    root = data_dir(project)
    design_path = candidate_design_matrix_path(project)
    fastq_dir = Path(str(project.get("fastq_dir", ""))).expanduser() if project.get("fastq_dir") else root / "fastq"
    rows = [
        {
            "step": "Setup",
            "status": "Complete" if project.get("name") and project.get("results_root") else "Needs attention",
            "evidence": str(project_root(project)),
            "count": "",
        },
        {
            "step": "Design matrix",
            "status": "Complete" if design_path.exists() else "Missing",
            "evidence": str(design_path),
            "count": "",
        },
        {
            "step": "FASTQ reads",
            "status": "Complete" if fastq_dir.is_dir() and len(fastq_files(str(fastq_dir))) > 0 else "Optional/missing",
            "evidence": str(fastq_dir),
            "count": len(fastq_files(str(fastq_dir))) if fastq_dir.is_dir() else 0,
        },
        {
            "step": "FastQC",
            "status": "Complete" if count_files(root / "fastqc", ["*.html"]) > 0 or count_files(root / "fastqc_cutadapt", ["*.html"]) > 0 else "Not found",
            "evidence": str(root / "fastqc"),
            "count": count_files(root / "fastqc", ["*.html"]) + count_files(root / "fastqc_cutadapt", ["*.html"]),
        },
        {
            "step": "Cutadapt",
            "status": "Complete" if count_files(root / "cutadapt", ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"]) > 0 else "Not found",
            "evidence": str(root / "cutadapt"),
            "count": count_files(root / "cutadapt", ["*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"]),
        },
        {
            "step": "STAR",
            "status": "Complete" if count_files(root / "star", ["*Aligned.sortedByCoord.out.bam"]) > 0 else "Not found",
            "evidence": str(root / "star"),
            "count": count_files(root / "star", ["*Aligned.sortedByCoord.out.bam"]),
        },
        {
            "step": "Kallisto",
            "status": "Complete" if count_files(root / "kallisto", ["abundance.tsv"]) > 0 else "Not found",
            "evidence": str(root / "kallisto"),
            "count": count_files(root / "kallisto", ["abundance.tsv"]),
        },
        {
            "step": "featureCounts",
            "status": "Complete" if count_files(root / "featurecounts", ["*_counts.txt"]) > 0 else "Not found",
            "evidence": str(root / "featurecounts"),
            "count": count_files(root / "featurecounts", ["*_counts.txt"]),
        },
        {
            "step": "Count matrix",
            "status": "Complete" if (root / "counts" / "count_matrix.txt").exists() else "Not found",
            "evidence": str(root / "counts" / "count_matrix.txt"),
            "count": "",
        },
        {
            "step": "DESeq2",
            "status": "Complete" if count_files(root / "deseq2", ["DEG*.txt", "*normalized*.txt"]) > 0 else "Not found",
            "evidence": str(root / "deseq2"),
            "count": count_files(root / "deseq2", ["DEG*.txt", "*normalized*.txt"]),
        },
        {
            "step": "Pathway analysis",
            "status": "Complete" if count_files(root / "gseapy", ["*.csv", "*.txt", "*.png", "*.pdf"]) > 0 else "Not found",
            "evidence": str(root / "gseapy"),
            "count": count_files(root / "gseapy", ["*.csv", "*.txt", "*.png", "*.pdf"]),
        },
    ]
    df = pd.DataFrame(rows)
    df["count"] = df["count"].astype(str)
    df.insert(0, "ready", df["status"].map(lambda x: "yes" if x == "Complete" else ""))
    return df


def next_recommended_step(project: dict) -> str:
    status = project_step_status(project).set_index("step")["status"].to_dict()
    if status.get("Design matrix") != "Complete":
        return "Design Matrix"
    if status.get("STAR") != "Complete" and status.get("Kallisto") != "Complete":
        if status.get("FastQC") != "Complete" and status.get("FASTQ reads") == "Complete":
            return "FastQC"
        return "STAR"
    if status.get("featureCounts") != "Complete" and status.get("Count matrix") != "Complete":
        return "featureCounts"
    if status.get("Count matrix") != "Complete":
        return "Count matrix"
    if status.get("DESeq2") != "Complete":
        return "DESeq2"
    return "Results Explorer"


def render_step_status(project: dict) -> None:
    status = project_step_status(project)
    st.dataframe(
        status,
        use_container_width=True,
        hide_index=True,
        column_config={
            "ready": st.column_config.TextColumn("ready", width="small"),
            "step": st.column_config.TextColumn("step", width="medium"),
            "status": st.column_config.TextColumn("status", width="small"),
            "count": st.column_config.TextColumn("count", width="small"),
            "evidence": st.column_config.TextColumn("path/evidence", width="large"),
        },
    )


def resume_step_details(project: dict, step: str) -> Tuple[str, List[str]]:
    root = data_dir(project)
    details = {
        "FastQC": (
            "Run read-level QC from raw or trimmed FASTQs.",
            [str(Path(project.get("fastq_dir", "")).expanduser()), str(root / "fastqc")],
        ),
        "Trim": (
            "Trim adapters with cutadapt and write cleaned reads.",
            [str(Path(project.get("fastq_dir", "")).expanduser()), str(root / "cutadapt")],
        ),
        "STAR": (
            "Align FASTQs or trimmed reads to the selected genome.",
            [str(candidate_design_matrix_path(project)), str(root / "star")],
        ),
        "Kallisto": (
            "Quantify transcript abundance from FASTQs or trimmed reads.",
            [str(candidate_design_matrix_path(project)), str(root / "kallisto")],
        ),
        "featureCounts": (
            "Count aligned STAR BAM files by gene_id or gene_name.",
            [str(root / "star"), str(root / "featurecounts")],
        ),
        "Count matrix": (
            "Merge featureCounts sample files into count_matrix.txt.",
            [str(root / "featurecounts"), str(root / "counts" / "count_matrix.txt")],
        ),
        "DESeq2": (
            "Run differential expression from count_matrix.txt and design_matrix.txt.",
            [str(root / "counts" / "count_matrix.txt"), str(candidate_design_matrix_path(project)), str(root / "deseq2")],
        ),
        "Results Explorer": (
            "Open the integrated Streamlit results viewer for completed outputs.",
            [str(root), str(candidate_design_matrix_path(project))],
        ),
    }
    return details.get(step, ("Resume from this step.", [str(root)]))


def render_resume_card(project: dict, step: str) -> None:
    description, paths = resume_step_details(project, step)
    st.markdown("**Resume guidance**")
    st.caption(description)
    st.code("\n".join(paths))


def split_fastq_suffix(filename: str) -> Tuple[str, str]:
    name = Path(str(filename).strip()).name
    lower = name.lower()
    for suffix in FASTQ_SUFFIXES:
        if lower.endswith(suffix):
            return name[:-len(suffix)], name[-len(suffix):]
    return name, ""


def mate_fastq_name(filename: str, mate: str) -> Optional[str]:
    stem, suffix = split_fastq_suffix(filename)
    if str(mate) == "2":
        replacements = [
            (r"([._-]R)1([._-]?\d*)$", r"\g<1>2\2"),
            (r"([._-])1$", r"\g<1>2"),
        ]
    else:
        replacements = [
            (r"([._-]R)2([._-]?\d*)$", r"\g<1>1\2"),
            (r"([._-])2$", r"\g<1>1"),
        ]
    for pattern, repl in replacements:
        new_stem, n = re.subn(pattern, repl, stem, flags=re.IGNORECASE)
        if n:
            return new_stem + suffix
    return None


def infer_sample_name(filename: str) -> str:
    stem, _suffix = split_fastq_suffix(filename)
    stem = re.sub(r"([._-]R)[12]([._-]?\d*)$", "", stem, flags=re.IGNORECASE)
    stem = re.sub(r"([._-])[12]$", "", stem)
    return clean_name(stem)


def fastq_files(folder: str) -> List[str]:
    path = Path(folder).expanduser()
    if not path.is_dir():
        return []
    return sorted([
        p.name for p in path.iterdir()
        if p.is_file() and p.name.lower().endswith(tuple(FASTQ_SUFFIXES))
    ])


def scan_fastqs(folder: str, paired: bool) -> pd.DataFrame:
    files = fastq_files(folder)
    file_set = set(files)
    rows = []
    used = set()
    if paired:
        for r1 in files:
            r2 = mate_fastq_name(r1, "2")
            if not r2:
                continue
            if r2 in file_set:
                rows.append({
                    "include": True,
                    "sample": infer_sample_name(r1),
                    "filename": f"{r1},{r2}",
                    "detected_status": "paired",
                })
                used.add(r1)
                used.add(r2)
            else:
                rows.append({
                    "include": False,
                    "sample": infer_sample_name(r1),
                    "filename": r1,
                    "detected_status": "missing R2",
                })
                used.add(r1)
    else:
        for name in files:
            mate1 = mate_fastq_name(name, "1")
            if mate1 and mate1 in file_set:
                continue
            rows.append({
                "include": True,
                "sample": infer_sample_name(name),
                "filename": name,
                "detected_status": "single",
            })
            used.add(name)
    return pd.DataFrame(rows)


def sample_fastq_pairs(project: dict, trimmed: bool = False) -> List[Tuple[str, Path, Path]]:
    design = read_design(project)
    base = data_dir(project) / "cutadapt" if trimmed else Path(project["fastq_dir"]).expanduser()
    pairs = []
    paired = bool(project.get("paired_end", True))
    for _, row in design.iterrows():
        parts = [x.strip() for x in str(row["filename"]).split(",") if x.strip()]
        if not parts:
            continue
        r1_name = Path(parts[0]).name if trimmed else parts[0]
        r1 = Path(r1_name).expanduser() if Path(r1_name).is_absolute() else base / r1_name
        if paired:
            if len(parts) > 1:
                r2_name = Path(parts[1]).name if trimmed else parts[1]
            else:
                inferred = mate_fastq_name(parts[0], "2")
                r2_name = Path(inferred).name if (trimmed and inferred) else inferred
            r2 = Path(r2_name).expanduser() if r2_name and Path(r2_name).is_absolute() else base / str(r2_name)
        else:
            r2 = r1
        pairs.append((str(row["sample"]), r1, r2))
    return pairs


def write_design(project: dict, edited: pd.DataFrame, metadata_cols: List[str]) -> Path:
    out = design_matrix_path(project)
    out.parent.mkdir(parents=True, exist_ok=True)
    keep = edited[edited["include"] == True].copy()
    if keep.empty:
        raise ValueError("No samples are included.")
    keep["sample"] = keep["sample"].map(clean_name)
    columns = ["sample"] + metadata_cols + ["filename"]
    for col in metadata_cols:
        if col not in keep.columns:
            keep[col] = "NA"
        keep[col] = keep[col].fillna("NA").astype(str).str.replace(r"\s+", "_", regex=True)
    keep[columns].to_csv(out, sep="\t", index=False)
    project["design_matrix_path"] = str(out)
    save_project(project)
    return out


def read_design(project: dict) -> pd.DataFrame:
    path = design_matrix_path(project)
    if path.exists():
        return pd.read_table(path)
    return pd.DataFrame()

def fastqc_html_name(read_name: str) -> str:
    stem, _suffix = split_fastq_suffix(Path(str(read_name)).name)
    return stem + "_fastqc.html"


def path_status(path: Path) -> str:
    return "ready" if path.exists() else "missing"


def sample_progress(project: dict) -> pd.DataFrame:
    design = read_design(project)
    if design.empty:
        return pd.DataFrame(columns=["sample", "FastQC", "Trim", "STAR", "Kallisto", "featureCounts"])
    root = data_dir(project)
    rows = []
    for _, row in design.iterrows():
        sample = str(row.get("sample", "")).strip()
        reads = [x.strip() for x in str(row.get("filename", "")).split(",") if x.strip()]
        read_names = [Path(x).name for x in reads]
        raw_fastqc = all((root / "fastqc" / fastqc_html_name(x)).exists() for x in read_names) if read_names else False
        trimmed_fastqc = all((root / "fastqc_cutadapt" / fastqc_html_name(x)).exists() for x in read_names) if read_names else False
        trimmed = all((root / "cutadapt" / Path(x).name).exists() for x in read_names) if read_names else False
        star_bam = root / "star" / sample / f"{sample}Aligned.sortedByCoord.out.bam"
        kallisto_abundance = root / "kallisto" / sample / "abundance.tsv"
        featurecounts_file = root / "featurecounts" / sample / f"{sample}_counts.txt"
        rows.append({
            "sample": sample,
            "FastQC": "ready" if raw_fastqc or trimmed_fastqc else "missing",
            "Trim": "ready" if trimmed else "missing",
            "STAR": path_status(star_bam),
            "Kallisto": path_status(kallisto_abundance),
            "featureCounts": path_status(featurecounts_file),
        })
    return pd.DataFrame(rows)


def run_selected_step(project: dict, step: str, use_trimmed: bool = False, feature: str = "gene_id", reference: str = "", comparison: str = "", redundant: str = "NoRedundant"):
    if step == "FastQC":
        return submit_fastqc(project, trimmed=use_trimmed)
    if step == "Trim":
        return submit_cutadapt(
            project,
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            "20",
        )
    if step == "STAR":
        return submit_star(project, use_trimmed=use_trimmed)
    if step == "Kallisto":
        return submit_kallisto(project, use_trimmed=use_trimmed)
    if step == "featureCounts":
        return submit_featurecounts(project, feature=feature)
    if step == "Count matrix":
        out = create_featurecounts_matrix(project)
        job = {
            "project": project["name"],
            "project_id": project_id(project),
            "analysis_type": project.get("analysis_type", "RNA-seq"),
            "step": "Count matrix",
            "command": "create_featurecounts_matrix",
            "stdout": "",
            "stderr": "",
            "submitted_at": now_stamp(),
            "job_id": None,
            "return_code": 0,
            "submit_output": f"Wrote {out}",
        }
        save_job(job)
        return job
    if step == "DESeq2":
        return submit_deseq2(project, reference, comparison, redundant)
    raise ValueError(f"Unsupported step: {step}")


def progress_tab(project: dict) -> None:
    st.subheader("Progress")
    st.caption("Select an analysis, open a project config, review completed outputs, then run or resume one step at a time.")
    c1, c2, c3 = st.columns(3)
    c1.metric("Analysis", project.get("analysis_type", "RNA-seq"))
    c2.metric("Project", project.get("name", "project"))
    c3.metric("Recommended next step", next_recommended_step(project))

    st.markdown("**Pipeline Status**")
    render_step_status(project)

    st.markdown("**Sample Progress**")
    sample_df = sample_progress(project)
    if sample_df.empty:
        st.info("No design matrix was found yet, so sample-level progress cannot be shown.")
    else:
        st.dataframe(sample_df, use_container_width=True, hide_index=True, height=360)
        st.download_button(
            "Download sample progress",
            data=sample_df.to_csv(index=False).encode(),
            file_name=f"{project.get('name', 'project')}_sample_progress.csv",
            key="download_sample_progress_"+project_id(project),
        )

    if project.get("analysis_type") != "RNA-seq":
        st.info("Run buttons are implemented for RNA-seq first. ATAC-seq and ChIP-seq configs still get status detection and output browsing.")
        return

    st.markdown("**Run One Step**")
    steps = ["FastQC", "Trim", "STAR", "Kallisto", "featureCounts", "Count matrix", "DESeq2", "Results Explorer"]
    recommended = next_recommended_step(project)
    default = steps.index(recommended) if recommended in steps else 0
    selected_step = st.selectbox("Step", steps, index=default, key="progress_step_"+project_id(project))
    if selected_step == "Results Explorer":
        integrated_results_explorer(project, key_prefix="progress_tab_"+project_id(project), compact=True)
        return

    use_trimmed = False
    feature = "gene_id"
    reference = ""
    comparison = ""
    redundant = "NoRedundant"
    if selected_step in ["FastQC", "STAR", "Kallisto"]:
        use_trimmed = st.toggle("Use trimmed reads", value=selected_step != "FastQC", key="progress_trimmed_"+project_id(project))
    if selected_step == "featureCounts":
        feature = st.selectbox("Feature attribute", ["gene_id", "gene_name"], key="progress_feature_"+project_id(project))
    if selected_step == "DESeq2":
        design = read_design(project)
        metadata_cols = [c for c in design.columns if c not in ["sample", "filename"]]
        if metadata_cols:
            design_col = st.selectbox("Design column", metadata_cols, key="progress_deseq_column_"+project_id(project))
            choices = sorted([x for x in design[design_col].dropna().astype(str).unique().tolist() if x])
            col1, col2, col3 = st.columns(3)
            with col1:
                reference = st.selectbox("Reference", choices or ["control"], key="progress_deseq_ref_"+project_id(project))
            with col2:
                comparison = st.selectbox("Comparison", choices or ["treated"], key="progress_deseq_comp_"+project_id(project))
            with col3:
                redundant = st.selectbox("Redundant covariate", ["NoRedundant"] + [c for c in metadata_cols if c != design_col], key="progress_deseq_redundant_"+project_id(project))
        else:
            st.warning("DESeq2 needs a design matrix with at least one metadata column.")
            return

    if st.button("Run selected step", type="primary", key="progress_run_"+project_id(project)):
        try:
            if selected_step == "DESeq2" and reference == comparison:
                st.error("Reference and comparison must be different.")
            else:
                job_submission_result(run_selected_step(project, selected_step, use_trimmed, feature, reference, comparison, redundant))
        except Exception as exc:
            st.error(str(exc))



def shell_quote(args: Iterable[object]) -> str:
    return " ".join(subprocess.list2cmdline([str(a)]) for a in args)


def command_exists(command: str) -> bool:
    return shutil.which(command) is not None


def scheduler_available() -> bool:
    return command_exists("sbatch")


def format_size(num_bytes: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    value = float(num_bytes)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024
    return str(num_bytes)


def run_command(args: List[object], cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    return subprocess.run(
        [str(x) for x in args],
        cwd=str(cwd or REPO_ROOT),
        text=True,
        capture_output=True,
        check=False,
    )


def parse_job_id(output: str) -> Optional[str]:
    match = re.search(r"Submitted batch job\s+(\d+)", output or "")
    if match:
        return match.group(1)
    return None


def submit_sbatch(project: dict, step: str, script: Path, script_args: List[object], log_name: str) -> dict:
    log_dir(project).mkdir(parents=True, exist_ok=True)
    stdout = log_dir(project) / f"output_{log_name}.txt"
    stderr = log_dir(project) / f"error_{log_name}.txt"
    cmd = ["sbatch", "-e", stderr, "-o", stdout, script] + script_args
    job = {
        "project": project["name"],
        "project_id": project_id(project),
        "analysis_type": project.get("analysis_type", "RNA-seq"),
        "step": step,
        "command": shell_quote(cmd),
        "stdout": str(stdout),
        "stderr": str(stderr),
        "submitted_at": now_stamp(),
        "job_id": None,
        "return_code": None,
        "submit_output": "",
    }
    if not command_exists("sbatch"):
        job["submit_output"] = "sbatch was not found on this machine. Run this app on the server to submit jobs."
        save_job(job)
        return job
    result = run_command(cmd, cwd=REPO_ROOT)
    output = (result.stdout or "") + (result.stderr or "")
    job["return_code"] = result.returncode
    job["submit_output"] = output.strip()
    job["job_id"] = parse_job_id(output)
    save_job(job)
    return job


def scheduler_status(job_id: Optional[str]) -> str:
    if not job_id:
        return "not submitted"
    if command_exists("squeue"):
        result = run_command(["squeue", "-j", job_id, "-h", "-o", "%T"])
        status = (result.stdout or "").strip()
        if status:
            return status
    if command_exists("sacct"):
        result = run_command(["sacct", "-j", job_id, "--format=State", "-n", "-P"])
        status = (result.stdout or "").strip().splitlines()
        if status:
            return status[0].split("|")[0]
    return "completed/unknown"


def read_tail(path: str, n: int = 120) -> str:
    p = Path(path)
    if not p.exists():
        return ""
    lines = p.read_text(errors="replace").splitlines()
    return "\n".join(lines[-n:])


def genome_resources(project: dict) -> dict:
    genome = str(project.get("genome", "mouse")).lower()
    return GENOME_RESOURCES.get(genome, GENOME_RESOURCES["mouse"])


def submit_fastqc(project: dict, trimmed: bool = False) -> List[dict]:
    outdir = data_dir(project) / ("fastqc_cutadapt" if trimmed else "fastqc")
    outdir.mkdir(parents=True, exist_ok=True)
    files = []
    for _sample, r1, r2 in sample_fastq_pairs(project, trimmed=trimmed):
        files.append(r1)
        if project.get("paired_end", True):
            files.append(r2)
    jobs = []
    for read in sorted(set(files)):
        jobs.append(submit_sbatch(
            project,
            "FastQC",
            SCRIPTS / "FastQC" / "qsub_fastqc.sh",
            [read, outdir, project["name"]],
            "fastQC",
        ))
    return jobs


def submit_cutadapt(project: dict, adapter1: str, adapter2: str, min_length: str) -> List[dict]:
    outdir = data_dir(project) / "cutadapt"
    outdir.mkdir(parents=True, exist_ok=True)
    paired = bool(project.get("paired_end", True))
    script = SCRIPTS / ("cutadapt_PE/qsub_cutadapt_PE.sh" if paired else "cutadapt_SE/qsub_cutadapt_SE.sh")
    jobs = []
    for _sample, r1, r2 in sample_fastq_pairs(project, trimmed=False):
        trimmed1 = outdir / r1.name
        trimmed2 = outdir / r2.name
        jobs.append(submit_sbatch(
            project,
            "Cutadapt",
            script,
            [min_length, adapter1, adapter2, trimmed1, trimmed2, r1, r2, project["name"]],
            "cutadapt",
        ))
    return jobs


def submit_star(project: dict, use_trimmed: bool = False) -> List[dict]:
    resources = genome_resources(project)
    outdir = data_dir(project) / "star"
    outdir.mkdir(parents=True, exist_ok=True)
    paired = bool(project.get("paired_end", True))
    script = SCRIPTS / ("STAR/qsub_star_PE.sh" if paired else "STAR/qsub_star_SE.sh")
    jobs = []
    for sample, r1, r2 in sample_fastq_pairs(project, trimmed=use_trimmed):
        sample_dir = outdir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        out_prefix = sample_dir / sample
        jobs.append(submit_sbatch(
            project,
            "STAR",
            script,
            [out_prefix, resources["star_index"], r1, r2, project["name"]],
            "star",
        ))
    return jobs


def submit_kallisto(project: dict, use_trimmed: bool = False) -> List[dict]:
    resources = genome_resources(project)
    outdir = data_dir(project) / "kallisto"
    outdir.mkdir(parents=True, exist_ok=True)
    paired = bool(project.get("paired_end", True))
    script = SCRIPTS / ("Kallisto/qsub_kallisto_PE.sh" if paired else "Kallisto/qsub_kallisto_SE.sh")
    jobs = []
    for sample, r1, r2 in sample_fastq_pairs(project, trimmed=use_trimmed):
        sample_dir = outdir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        jobs.append(submit_sbatch(
            project,
            "Kallisto",
            script,
            [sample_dir, resources["kallisto_index"], r1, r2, project["name"]],
            "kallisto",
        ))
    return jobs


def submit_featurecounts(project: dict, feature: str = "gene_id") -> List[dict]:
    resources = genome_resources(project)
    outdir = data_dir(project) / "featurecounts"
    outdir.mkdir(parents=True, exist_ok=True)
    paired = bool(project.get("paired_end", True))
    script = SCRIPTS / ("featureCounts/qsub_featurecounts_PE.sh" if paired else "featureCounts/qsub_featurecounts_SE.sh")
    jobs = []
    for sample in read_design(project)["sample"].astype(str).tolist():
        sample_dir = outdir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)
        bam = data_dir(project) / "star" / sample / f"{sample}Aligned.sortedByCoord.out.bam"
        count_prefix = sample_dir / sample
        jobs.append(submit_sbatch(
            project,
            "featureCounts",
            script,
            [bam, resources["gtf"], feature, count_prefix, resources["strand_bed"], project["name"]],
            "featurecounts",
        ))
    return jobs


def create_featurecounts_matrix(project: dict) -> Path:
    inpath = data_dir(project) / "featurecounts"
    outpath = data_dir(project) / "counts"
    outpath.mkdir(parents=True, exist_ok=True)
    matrices = []
    for sample_dir in sorted([p for p in inpath.iterdir() if p.is_dir()]):
        count_file = sample_dir / f"{sample_dir.name}_counts.txt"
        if not count_file.exists():
            continue
        df = pd.read_table(count_file, comment="#", index_col=0)
        drop_cols = [c for c in ["Chr", "Start", "End", "Strand", "Length"] if c in df.columns]
        df = df.drop(columns=drop_cols)
        if df.shape[1] > 0:
            df = df.rename(columns={df.columns[0]: sample_dir.name})
            matrices.append(df[[sample_dir.name]])
    if not matrices:
        raise FileNotFoundError("No featureCounts sample count files were found.")
    count_matrix = pd.concat(matrices, axis=1)
    out_file = outpath / "count_matrix.txt"
    count_matrix.to_csv(out_file, sep="\t")
    return out_file


def submit_deseq2(project: dict, reference: str, comparison: str, redundant: str = "NoRedundant") -> dict:
    outpath = data_dir(project) / "deseq2"
    outpath.mkdir(parents=True, exist_ok=True)
    count_matrix = data_dir(project) / "counts" / "count_matrix.txt"
    return submit_sbatch(
        project,
        "DESeq2",
        SCRIPTS / "DESeq2" / "qsub_deseq2.sh",
        [SCRIPTS / "DESeq2" / "DESeq2.R", count_matrix, design_matrix_path(project), outpath, reference, comparison, redundant or "NoRedundant", project["name"]],
        "deseq2",
    )


def available_port(start: int = 3838, end: int = 3900) -> int:
    for port in range(start, end + 1):
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            sock.bind(("127.0.0.1", port))
            return port
        except OSError:
            continue
        finally:
            sock.close()
    return start


def style() -> None:
    st.markdown(
        """
        <style>
        :root {
            --csl-ink:#17202f;
            --csl-muted:#5f6f85;
            --csl-line:#d8dee8;
            --csl-bg:#f6f8fb;
            --csl-blue:#1d4ed8;
            --csl-green:#0f766e;
        }
        .stApp { background: var(--csl-bg); color: var(--csl-ink); }
        h1, h2, h3 { letter-spacing: 0; color: var(--csl-ink); }
        h1 { font-size: 2.4rem !important; line-height: 1.08 !important; }
        h2 { font-size: 1.7rem !important; }
        h3 { font-size: 1.35rem !important; }
        [data-testid="stSidebar"] {
            background: #ffffff;
            border-right: 1px solid var(--csl-line);
            min-width: 300px !important;
            max-width: 360px !important;
        }
        [data-testid="stSidebar"] h1 { font-size: 1.75rem !important; }
        input,
        textarea,
        [data-baseweb="input"] input,
        [data-baseweb="select"] > div,
        [data-baseweb="textarea"] textarea {
            background: #ffffff !important;
            border-color: #c7d0dd !important;
            color: var(--csl-ink) !important;
        }
        [data-baseweb="select"] span,
        [data-baseweb="select"] div {
            color: var(--csl-ink) !important;
        }
        [data-testid="stWidgetLabel"] p,
        [data-testid="stSidebar"] label,
        [data-testid="stSidebar"] p {
            color: var(--csl-ink) !important;
        }
        [data-testid="stAlert"] p,
        [data-testid="stAlert"] div {
            color: var(--csl-ink) !important;
        }
        div[data-testid="stMetric"] {
            background: #ffffff;
            border: 1px solid var(--csl-line);
            border-radius: 8px;
            padding: 14px 16px;
        }
        [data-testid="stMetricLabel"] {
            color: #5b6678 !important;
            opacity: 1 !important;
        }
        [data-testid="stMetricValue"] {
            color: var(--csl-ink) !important;
            overflow-wrap: anywhere;
        }
        .csl-header {
            background: #ffffff;
            border: 1px solid var(--csl-line);
            border-radius: 8px;
            padding: 16px 18px;
            margin-bottom: 14px;
            overflow: hidden;
        }
        .csl-header h1 {
            font-size: 2.15rem !important;
            line-height: 1.05 !important;
            white-space: normal;
        }
        .csl-subtle { color: var(--csl-muted); font-size: 0.94rem; }
        .csl-badge {
            display:inline-block;
            border:1px solid var(--csl-line);
            border-radius:999px;
            padding: 4px 10px;
            margin-right: 6px;
            background:#ffffff;
            color:var(--csl-muted);
            font-size: 0.82rem;
        }
        .stButton>button {
            border-radius: 6px;
            border: 1px solid #b9c4d4;
            background: #ffffff;
            color: var(--csl-ink);
            font-weight: 600;
        }
        .stButton>button[kind="primary"] {
            background: var(--csl-blue);
            color: white;
            border-color: var(--csl-blue);
        }
        button[data-baseweb="tab"] p {
            color: var(--csl-muted) !important;
            font-weight: 650 !important;
        }
        button[data-baseweb="tab"][aria-selected="true"] p {
            color: var(--csl-blue) !important;
        }
        div[data-testid="stDataFrame"],
        div[data-testid="stDataEditor"] {
            border: 1px solid var(--csl-line);
            border-radius: 8px;
            overflow: hidden;
            background: #ffffff;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def header(project: Optional[dict]) -> None:
    if project:
        st.markdown(
            f"""
            <div class="csl-header">
              <h1 style="margin:0;">CodeSpringLab Control Center</h1>
              <div class="csl-subtle">Project: <b>{project["name"]}</b> · {project.get("analysis_type","RNA-seq")} · {project.get("genome","mouse")} · {"paired-end" if project.get("paired_end", True) else "single-end"}</div>
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            """
            <div class="csl-header">
              <h1 style="margin:0;">CodeSpringLab Control Center</h1>
              <div class="csl-subtle">Create or select a project to scan reads, build design matrices, submit jobs, and inspect outputs.</div>
            </div>
            """,
            unsafe_allow_html=True,
        )


def sidebar_project_selector() -> Optional[dict]:
    st.sidebar.title("Projects")
    analysis_options = ["RNA-seq", "ATAC-seq", "ChIP-seq", "All analyses"]
    selected_analysis = st.sidebar.selectbox("Analysis", analysis_options, key="sidebar_analysis")
    projects = load_projects()
    filtered = {
        pid: project for pid, project in projects.items()
        if selected_analysis == "All analyses" or project.get("analysis_type", "RNA-seq") == selected_analysis
    }
    project_ids = sorted(
        filtered,
        key=lambda pid: filtered[pid].get("updated_at", filtered[pid].get("created_at", "")),
        reverse=True,
    )
    options = ["New project"] + project_ids
    default_project = project_ids[0] if project_ids else "New project"
    selected = st.session_state.pop("project_to_select", st.session_state.get("selected_project", default_project))
    if selected not in options:
        selected = default_project

    def label_project(option: str) -> str:
        if option == "New project":
            return "New project"
        project = filtered[option]
        source = " · CSL config" if project.get("source_config") else ""
        return f"{project.get('name')} ({project.get('analysis_type', 'RNA-seq')}{source})"

    choice = st.sidebar.selectbox(
        "Project config",
        options,
        index=options.index(selected),
        key="selected_project",
        format_func=label_project,
    )
    st.sidebar.caption(f"CodeSpringLab root: {REPO_ROOT}")
    legacy_count = len([p for p in filtered.values() if p.get("source_config")])
    if legacy_count:
        st.sidebar.caption(f"Detected {legacy_count} CodeSpringLab project config(s).")
    if choice != "New project":
        return filtered[choice]

    with st.sidebar.form("new_project"):
        workflow_mode = st.selectbox(
            "Workflow",
            ["Start new analysis", "Resume existing analysis", "Visualize existing results"],
        )
        name = st.text_input("Project name", value="example_dataset")
        default_analysis = selected_analysis if selected_analysis != "All analyses" else "RNA-seq"
        analysis_type = st.selectbox(
            "Analysis type",
            ["RNA-seq", "ATAC-seq", "ChIP-seq"],
            index=["RNA-seq", "ATAC-seq", "ChIP-seq"].index(default_analysis),
        )
        genome = st.selectbox("Genome", ["mouse", "human"])
        paired = st.toggle("Paired-end reads", value=True)
        results_root = st.text_input("Results root", value=str(Path("~/csl_results").expanduser()))
        fastq_dir = st.text_input("FASTQ folder", value="")
        design_path_input = st.text_input("Design matrix path or folder", value="")
        submitted = st.form_submit_button("Create / import project", type="primary")
    if submitted:
        project = {
            "name": clean_name(name, "project"),
            "analysis_type": analysis_type,
            "workflow_mode": workflow_mode,
            "genome": genome,
            "paired_end": paired,
            "results_root": str(Path(results_root).expanduser()),
            "fastq_dir": str(Path(fastq_dir).expanduser()) if fastq_dir else "",
            "design_matrix_path": str(Path(design_path_input).expanduser()) if design_path_input else "",
            "metadata_columns": ["treatment"],
        }
        project = apply_project_inference(project)
        project = save_project(project)
        st.session_state["project_to_select"] = project_id(project)
        st.sidebar.success("Project saved.")
        st.rerun()
    return None


def setup_tab(project: dict) -> dict:
    st.subheader("Project Setup")
    st.caption("Use this page to start new analyses, import older projects, or configure a visualize-only project.")
    with st.form("project_setup_form"):
        col1, col2 = st.columns(2)
        with col1:
            modes = ["Start new analysis", "Resume existing analysis", "Visualize existing results"]
            project["workflow_mode"] = st.selectbox("Workflow", modes, index=modes.index(project.get("workflow_mode", modes[0])) if project.get("workflow_mode", modes[0]) in modes else 0)
            project["analysis_type"] = st.selectbox("Analysis type", ["RNA-seq", "ATAC-seq", "ChIP-seq"], index=["RNA-seq", "ATAC-seq", "ChIP-seq"].index(project.get("analysis_type", "RNA-seq")))
            project["genome"] = st.selectbox("Genome", ["mouse", "human"], index=["mouse", "human"].index(project.get("genome", "mouse")))
            project["paired_end"] = st.toggle("Paired-end reads", value=bool(project.get("paired_end", True)))
        with col2:
            project["results_root"] = st.text_input("Results root", value=str(project.get("results_root", Path("~/csl_results").expanduser())))
            project["fastq_dir"] = st.text_input("FASTQ folder", value=str(project.get("fastq_dir", "")))
            project["design_matrix_path"] = st.text_input("Design matrix path or folder", value=str(project.get("design_matrix_path", "")))
            if project.get("data_dir_override"):
                project["data_dir_override"] = st.text_input("Visualizer data folder", value=str(project.get("data_dir_override", "")))
            metadata = st.text_input("Design metadata columns", value=", ".join(project.get("metadata_columns", ["treatment"])))
        saved = st.form_submit_button("Save setup / re-detect paths", type="primary")
    if saved:
        project["metadata_columns"] = [clean_name(x, "condition") for x in metadata.split(",") if x.strip()] or ["treatment"]
        project = apply_project_inference(project)
        save_project(project)
        st.success("Project setup saved and paths re-detected.")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Design matrix", "ready" if candidate_design_matrix_path(project).exists() else "missing")
    c2.metric("FASTQ files", len(fastq_files(project.get("fastq_dir", ""))))
    c3.metric("Submitted jobs", len(project_jobs(project)))
    c4.metric("Next step", next_recommended_step(project))

    if project.get("source_config"):
        st.markdown("**Imported CodeSpringLab Config**")
        st.code(str(project.get("source_config")))
    st.markdown("**Detected Project State**")
    render_step_status(project)
    return project


def design_tab(project: dict) -> None:
    st.subheader("Design Matrix")
    st.caption("Scan the FASTQ folder to prefill filenames and inferred sample names, then type metadata directly into the table.")
    scan_key = "scan_df_"+project_id(project)
    cols_key = "metadata_cols_"+project_id(project)
    current_cols = project.get("metadata_columns", ["treatment"])

    col_a, col_b = st.columns([3, 1])
    with col_a:
        metadata_text = st.text_input(
            "Metadata columns",
            value=", ".join(st.session_state.get(cols_key, current_cols)),
            help="Comma-separated columns to add to design_matrix.txt, for example treatment, batch, replicate.",
            key="metadata_text_"+project_id(project),
        )
    with col_b:
        if st.button("Apply columns", key="apply_metadata_columns_"+project_id(project)):
            metadata_cols = parse_metadata_columns(metadata_text, current_cols)
            project["metadata_columns"] = metadata_cols
            save_project(project)
            st.session_state[cols_key] = metadata_cols
            if scan_key in st.session_state:
                st.session_state[scan_key] = sync_design_editor_columns(st.session_state[scan_key], metadata_cols)
            st.success("Metadata columns updated.")
            st.rerun()

    metadata_cols = st.session_state.get(cols_key, parse_metadata_columns(metadata_text, current_cols))
    project["metadata_columns"] = metadata_cols

    scan_col, empty_col = st.columns([1, 1])
    with scan_col:
        if st.button("Scan FASTQ folder / prefill filenames", type="primary"):
            scanned = scan_fastqs(project.get("fastq_dir", ""), bool(project.get("paired_end", True)))
            scanned = sync_design_editor_columns(scanned, metadata_cols)
            st.session_state[scan_key] = scanned
    with empty_col:
        if st.button("Start empty design table"):
            st.session_state[scan_key] = sync_design_editor_columns(pd.DataFrame(columns=["include", "sample", "filename", "detected_status"]), metadata_cols)

    existing = read_design(project)
    if scan_key not in st.session_state:
        if not existing.empty:
            display = existing.copy().fillna("")
            display.insert(0, "include", True)
            display["detected_status"] = "saved"
            st.session_state[scan_key] = sync_design_editor_columns(display, metadata_cols)
        else:
            st.session_state[scan_key] = sync_design_editor_columns(pd.DataFrame(), metadata_cols)
    else:
        st.session_state[scan_key] = sync_design_editor_columns(st.session_state[scan_key], metadata_cols)

    if st.session_state[scan_key].empty:
        st.info("Scan a FASTQ folder to prefill filenames, or start an empty design table and add rows manually.")

    display_cols = ["include", "sample"] + metadata_cols + ["filename", "detected_status"]
    column_config = {
        "include": st.column_config.CheckboxColumn("include", help="Uncheck to exclude a detected sample."),
        "sample": st.column_config.TextColumn("sample", width="medium", help="Edit sample names before saving."),
        "filename": st.column_config.TextColumn("FASTQ file(s)", width="large", help="Comma-separated R1,R2 filenames for paired-end projects."),
        "detected_status": st.column_config.TextColumn("status", width="small"),
    }
    for col in metadata_cols:
        column_config[col] = st.column_config.TextColumn(col, width="medium")

    edited = st.data_editor(
        st.session_state[scan_key],
        use_container_width=True,
        hide_index=True,
        num_rows="dynamic",
        column_order=[c for c in display_cols if c in st.session_state[scan_key].columns],
        disabled=["detected_status"],
        column_config=column_config,
    )
    st.session_state[scan_key] = sync_design_editor_columns(edited, metadata_cols)

    col1, col2 = st.columns([1, 3])
    with col1:
        if st.button("Save design_matrix.txt", type="primary"):
            try:
                project["metadata_columns"] = metadata_cols
                save_project(project)
                path = write_design(project, st.session_state[scan_key], metadata_cols)
                st.success(f"Saved {path}")
            except Exception as exc:
                st.error(str(exc))
    with col2:
        path = design_matrix_path(project)
        st.code(str(path))
    if path.exists():
        st.download_button(
            "Download design_matrix.txt",
            data=path.read_bytes(),
            file_name="design_matrix.txt",
            key="download_design_"+project_id(project),
        )


def job_submission_result(jobs) -> None:
    if not isinstance(jobs, list):
        jobs = [jobs]
    rows = []
    for job in jobs:
        rows.append({
            "step": job.get("step"),
            "job_id": job.get("job_id") or "",
            "status": scheduler_status(job.get("job_id")),
            "return_code": job.get("return_code"),
            "message": job.get("submit_output", "")[:220],
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    if any(not job.get("job_id") for job in jobs):
        st.caption("No scheduler job ID means the command was not submitted, usually because this is being tested away from the SLURM server.")


def run_tab(project: dict) -> None:
    st.subheader("Run Pipeline")
    if project.get("analysis_type") != "RNA-seq":
        st.info("This first web runner implements the RNA-seq run path. ATAC-seq and ChIP-seq projects can still use setup, design matrix, logs, and output browsing here.")
        return

    mode = project.get("workflow_mode", "Start new analysis")
    st.caption("Workflow: "+mode+" · recommended next step: "+next_recommended_step(project))
    with st.expander("Detected progress and resume guide", expanded=(mode != "Start new analysis")):
        render_step_status(project)
        resume_options = ["FastQC", "Trim", "STAR", "Kallisto", "featureCounts", "Count matrix", "DESeq2", "Results Explorer"]
        recommended = next_recommended_step(project)
        resume_default = resume_options.index(recommended) if recommended in resume_options else 0
        resume_step = st.selectbox(
            "I want to resume from",
            resume_options,
            index=resume_default,
        )
        render_resume_card(project, resume_step)
        st.info("Open the matching tab below and run only that step. Existing upstream outputs are read from the paths shown above.")

    if read_design(project).empty:
        st.warning("Save or import a design matrix before running analysis jobs. Visual outputs can still be browsed from the Outputs tab if files already exist.")
        if mode != "Visualize existing results":
            return
    if scheduler_available():
        st.success("SLURM scheduler detected. Run buttons will submit jobs with sbatch.")
    else:
        st.warning("Local preview mode: sbatch is not available here, so run buttons will record the command but will not submit jobs. Run this app on the server to execute the pipeline.")

    step_tabs = st.tabs(["FastQC", "Trim", "STAR", "Kallisto", "featureCounts", "DESeq2", "Results Explorer"])

    with step_tabs[0]:
        trimmed = st.toggle("Use trimmed reads", value=False, key="fastqc_trimmed")
        if st.button("Run FastQC", type="primary"):
            job_submission_result(submit_fastqc(project, trimmed=trimmed))

    with step_tabs[1]:
        a1 = st.text_input("R1 adapter", value="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")
        a2 = st.text_input("R2 adapter", value="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
        min_len = st.text_input("Minimum length", value="20")
        if st.button("Run cutadapt", type="primary"):
            job_submission_result(submit_cutadapt(project, a1, a2, min_len))

    with step_tabs[2]:
        use_trimmed = st.toggle("Use trimmed reads", value=False, key="star_trimmed")
        if st.button("Run STAR", type="primary"):
            job_submission_result(submit_star(project, use_trimmed=use_trimmed))

    with step_tabs[3]:
        use_trimmed = st.toggle("Use trimmed reads", value=False, key="kallisto_trimmed")
        if st.button("Run Kallisto", type="primary"):
            job_submission_result(submit_kallisto(project, use_trimmed=use_trimmed))

    with step_tabs[4]:
        feature = st.selectbox("Feature attribute", ["gene_id", "gene_name"])
        c1, c2 = st.columns(2)
        with c1:
            if st.button("Run featureCounts", type="primary"):
                job_submission_result(submit_featurecounts(project, feature=feature))
        with c2:
            if st.button("Build count_matrix.txt"):
                try:
                    out = create_featurecounts_matrix(project)
                    st.success(f"Wrote {out}")
                except Exception as exc:
                    st.error(str(exc))

    with step_tabs[5]:
        design = read_design(project)
        metadata_cols = [c for c in design.columns if c not in ["sample", "filename"]]
        if not metadata_cols:
            st.warning("DESeq2 needs at least one metadata column in the design matrix.")
        else:
            design_col = st.selectbox("Design column", metadata_cols, key="deseq_design_column_"+project_id(project))
            options = sorted([x for x in design[design_col].dropna().astype(str).unique().tolist() if x])
            col1, col2, col3 = st.columns(3)
            with col1:
                ref = st.selectbox("Reference", options or ["control"], key="deseq_reference_"+project_id(project))
            with col2:
                comp = st.selectbox("Comparison", options or ["treated"], key="deseq_comparison_"+project_id(project))
            with col3:
                redundant_options = ["NoRedundant"] + [c for c in metadata_cols if c != design_col]
                redundant = st.selectbox("Redundant covariate", redundant_options, key="deseq_redundant_"+project_id(project))
            if st.button("Run DESeq2", type="primary"):
                if ref == comp:
                    st.error("Reference and comparison must be different.")
                else:
                    job_submission_result(submit_deseq2(project, ref, comp, redundant))

    with step_tabs[6]:
        integrated_results_explorer(project, key_prefix="run_tab_"+project_id(project), compact=True)


def jobs_tab(project: dict) -> None:
    st.subheader("Jobs and Logs")
    jobs = project_jobs(project)
    if not jobs:
        st.info("No jobs submitted from the web app yet.")
        return
    rows = []
    for job in jobs:
        rows.append({
            "submitted": job.get("submitted_at"),
            "step": job.get("step"),
            "job_id": job.get("job_id") or "",
            "status": scheduler_status(job.get("job_id")),
            "return_code": job.get("return_code"),
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    selected = st.selectbox("Inspect job", list(range(len(jobs))), format_func=lambda i: f"{jobs[i].get('step')} {jobs[i].get('job_id') or ''} {jobs[i].get('submitted_at')}")
    job = jobs[selected]
    st.code(job.get("command", ""))
    col1, col2 = st.columns(2)
    with col1:
        st.caption(job.get("stdout", ""))
        st.text_area("stdout tail", value=read_tail(job.get("stdout", "")), height=320)
    with col2:
        st.caption(job.get("stderr", ""))
        st.text_area("stderr tail", value=read_tail(job.get("stderr", "")), height=320)


def safe_read_table(path: Path, preview_rows: Optional[int] = None) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    kwargs = {"sep": None, "engine": "python"}
    if preview_rows is not None:
        kwargs["nrows"] = preview_rows
    try:
        return pd.read_csv(path, **kwargs)
    except Exception:
        try:
            return pd.read_table(path, nrows=preview_rows)
        except Exception:
            return pd.DataFrame()


def download_df_button(df: pd.DataFrame, file_name: str, key: str) -> None:
    if df is None or df.empty:
        return
    st.download_button(
        "Download table",
        data=df.to_csv(index=False).encode(),
        file_name=file_name,
        key=key,
    )


def dataframe_view(path: Path, label: str, key: str, preview_rows: int = 1000, height: int = 420) -> None:
    st.markdown(f"**{label}**")
    if not path.exists():
        st.info(f"Not found: {path}")
        return
    df = safe_read_table(path, preview_rows=preview_rows)
    if df.empty:
        st.warning(f"Could not read table: {path}")
        st.download_button(
            f"Download {path.name}",
            data=path.read_bytes(),
            file_name=path.name,
            key=key+"_raw_download",
        )
        return
    st.caption(f"{path.name} · {format_size(path.stat().st_size)} · previewing up to {preview_rows} rows")
    st.dataframe(df, use_container_width=True, height=height)
    download_df_button(df, path.name.replace(".txt", ".csv"), key+"_df_download")


def matching_files(root: Path, patterns: Iterable[str], limit: int = 500) -> List[Path]:
    if not root.exists():
        return []
    found = []
    for pattern in patterns:
        found.extend(root.rglob(pattern))
    files = sorted([p for p in found if p.is_file()])
    return files[:limit]


def select_file(label: str, files: List[Path], key: str, root: Optional[Path] = None) -> Optional[Path]:
    if not files:
        st.info(f"No {label.lower()} found.")
        return None
    root = root or files[0].parent
    return st.selectbox(label, files, format_func=lambda p: p.relative_to(root).as_posix() if root in p.parents or p == root else p.name, key=key)


def render_file_inline(path: Path, key: str, height: int = 850) -> None:
    suffix = path.suffix.lower()
    if suffix == ".png":
        st.image(str(path), use_container_width=True)
    elif suffix in [".jpg", ".jpeg", ".webp"]:
        st.image(str(path), use_container_width=True)
    elif suffix == ".html":
        components.html(path.read_text(errors="replace"), height=height, scrolling=True)
    elif suffix == ".pdf":
        encoded = base64.b64encode(path.read_bytes()).decode()
        components.html(
            f'<iframe src="data:application/pdf;base64,{encoded}" width="100%" height="{height}" style="border:1px solid #d8dde8;border-radius:8px;"></iframe>',
            height=height + 20,
            scrolling=False,
        )
        st.download_button(f"Download {path.name}", data=path.read_bytes(), file_name=path.name, key=key+"_download_pdf")
    elif suffix in [".txt", ".csv", ".tsv"]:
        dataframe_view(path, path.name, key)
    else:
        st.download_button(f"Download {path.name}", data=path.read_bytes(), file_name=path.name, key=key+"_download")
        st.code(str(path))


def design_metadata(project: dict) -> pd.DataFrame:
    design = read_design(project)
    if design.empty:
        return design
    if "sample" not in design.columns:
        first = design.columns[0]
        design = design.rename(columns={first: "sample"})
    return design


def sorted_sample_options(project: dict, sort_col: Optional[str] = None) -> List[str]:
    design = design_metadata(project)
    if design.empty or "sample" not in design.columns:
        return []
    design = design.copy()
    design["sample"] = design["sample"].astype(str)
    if sort_col and sort_col in design.columns:
        design = design.sort_values([sort_col, "sample"], kind="stable")
    else:
        design = design.sort_values("sample", kind="stable")
    return design["sample"].tolist()


def fastqc_reports(project: dict) -> List[Path]:
    root = data_dir(project)
    return matching_files(root / "fastqc", ["*.html"]) + matching_files(root / "fastqc_cutadapt", ["*.html"])


def star_sample_files(project: dict, sample: str) -> List[Path]:
    root = data_dir(project) / "star"
    if not root.exists():
        return []
    candidates = []
    for sample_dir in [root / sample, root]:
        if sample_dir.exists():
            candidates.extend(matching_files(sample_dir, [f"*{sample}*Log.final.out", f"*{sample}*SJ.out.tab", f"*{sample}*.txt"], limit=50))
    return sorted(set(candidates))


def featurecounts_sample_files(project: dict, sample: str) -> List[Path]:
    root = data_dir(project) / "featurecounts"
    if not root.exists():
        return []
    return matching_files(root / sample, ["*.txt", "*.summary"], limit=50) + matching_files(root, [f"*{sample}*.txt", f"*{sample}*.summary"], limit=50)


def render_overview(project: dict, key_prefix: str) -> None:
    root = data_dir(project)
    rows = [
        {"resource": "Data folder", "status": "found" if root.exists() else "missing", "path": str(root)},
        {"resource": "Design matrix", "status": "found" if candidate_design_matrix_path(project).exists() else "missing", "path": str(candidate_design_matrix_path(project))},
        {"resource": "FastQC reports", "status": str(len(fastqc_reports(project))), "path": str(root / "fastqc")},
        {"resource": "STAR summary", "status": "found" if (root / "star_summary" / "summary_matrix.txt").exists() else "missing", "path": str(root / "star_summary" / "summary_matrix.txt")},
        {"resource": "featureCounts summary", "status": "found" if (root / "counts" / "featurecounts_summary.txt").exists() else "missing", "path": str(root / "counts" / "featurecounts_summary.txt")},
        {"resource": "Raw count matrix", "status": "found" if (root / "counts" / "count_matrix.txt").exists() else "missing", "path": str(root / "counts" / "count_matrix.txt")},
        {"resource": "DESeq2 files", "status": str(len(matching_files(root / "deseq2", ["*.txt", "*.csv", "*.png", "*.pdf"]))), "path": str(root / "deseq2")},
        {"resource": "GSEA files", "status": str(len(matching_files(root / "gseapy", ["*.txt", "*.csv", "*.png", "*.pdf"]))), "path": str(root / "gseapy")},
    ]
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True, height=310)
    design = design_metadata(project)
    if not design.empty:
        st.markdown("**Design Matrix**")
        st.dataframe(design, use_container_width=True, hide_index=True, height=260)
        download_df_button(design, f"{project.get('name', 'project')}_design_matrix.csv", key_prefix+"_download_design")


def render_qc_explorer(project: dict, key_prefix: str) -> None:
    root = data_dir(project)
    reports = fastqc_reports(project)
    selected = select_file("FastQC report", reports, key_prefix+"_fastqc", root)
    if selected:
        components.html(selected.read_text(errors="replace"), height=1050, scrolling=True)


def render_alignment_explorer(project: dict, key_prefix: str) -> None:
    root = data_dir(project)
    design = design_metadata(project)
    sort_cols = [c for c in design.columns if c not in ["sample", "filename"]] if not design.empty else []
    sort_col = st.selectbox("Sort samples by", ["sample"] + sort_cols, key=key_prefix+"_sort_samples") if sort_cols else "sample"
    samples = sorted_sample_options(project, None if sort_col == "sample" else sort_col)
    sample = st.selectbox("Sample", samples or [""], key=key_prefix+"_sample")
    col1, col2 = st.columns(2)
    with col1:
        dataframe_view(root / "star_summary" / "summary_matrix.txt", "STAR Summary Across Samples", key_prefix+"_star_summary")
        sample_files = star_sample_files(project, sample) if sample else []
        selected = select_file("Selected STAR sample file", sample_files, key_prefix+"_star_sample", root)
        if selected:
            render_file_inline(selected, key_prefix+"_star_sample_file", height=500)
    with col2:
        dataframe_view(root / "counts" / "featurecounts_summary.txt", "FeatureCounts Summary Across Samples", key_prefix+"_fc_summary")
        fc_files = featurecounts_sample_files(project, sample) if sample else []
        selected_fc = select_file("Selected featureCounts sample file", fc_files, key_prefix+"_fc_sample", root)
        if selected_fc:
            render_file_inline(selected_fc, key_prefix+"_fc_sample_file", height=500)


def render_counts_explorer(project: dict, key_prefix: str) -> None:
    root = data_dir(project)
    tabs = st.tabs(["Raw Counts", "RSEM", "Kallisto", "DESeq2 Normalized"])
    with tabs[0]:
        dataframe_view(root / "counts" / "count_matrix.txt", "Raw Count Matrix", key_prefix+"_raw_counts", preview_rows=5000, height=560)
    with tabs[1]:
        files = matching_files(root / "rsem", ["*.genes.results", "*.isoforms.results", "*matrix*.txt", "*.txt", "*.csv"], limit=200)
        selected = select_file("RSEM table", files, key_prefix+"_rsem", root)
        if selected:
            dataframe_view(selected, selected.relative_to(root).as_posix(), key_prefix+"_rsem_table", preview_rows=5000, height=560)
    with tabs[2]:
        files = matching_files(root / "kallisto", ["abundance.tsv", "*matrix*.txt", "*.tsv", "*.csv"], limit=300)
        selected = select_file("Kallisto table", files, key_prefix+"_kallisto", root)
        if selected:
            dataframe_view(selected, selected.relative_to(root).as_posix(), key_prefix+"_kallisto_table", preview_rows=5000, height=560)
    with tabs[3]:
        files = matching_files(root / "deseq2", ["*normalized*.txt", "*normalized*.csv"], limit=100)
        selected = select_file("DESeq2 normalized counts", files, key_prefix+"_norm_counts", root)
        if selected:
            dataframe_view(selected, selected.relative_to(root).as_posix(), key_prefix+"_norm_counts_table", preview_rows=5000, height=560)


def render_deseq2_explorer(project: dict, key_prefix: str) -> None:
    root = data_dir(project)
    deseq2 = root / "deseq2"
    tabs = st.tabs(["Tables", "Plots"])
    with tabs[0]:
        files = matching_files(deseq2, ["DEG*.txt", "DEG*.csv", "*.txt", "*.csv"], limit=300)
        selected = select_file("DESeq2 table", files, key_prefix+"_deseq_table", root)
        if selected:
            dataframe_view(selected, selected.relative_to(root).as_posix(), key_prefix+"_deseq_table_view", preview_rows=5000, height=620)
    with tabs[1]:
        files = matching_files(deseq2, ["*.png", "*.pdf"], limit=300)
        selected = select_file("DESeq2 plot", files, key_prefix+"_deseq_plot", root)
        if selected:
            render_file_inline(selected, key_prefix+"_deseq_plot_view", height=900)


def render_gsea_explorer(project: dict, key_prefix: str) -> None:
    root = data_dir(project)
    gseapy = root / "gseapy"
    tabs = st.tabs(["Tables", "Summary Plots", "Individual Pathway"])
    with tabs[0]:
        files = matching_files(gseapy, ["*.csv", "*.txt"], limit=500)
        selected = select_file("GSEA table", files, key_prefix+"_gsea_table", root)
        if selected:
            dataframe_view(selected, selected.relative_to(root).as_posix(), key_prefix+"_gsea_table_view", preview_rows=5000, height=620)
    with tabs[1]:
        files = matching_files(gseapy, ["*.png", "*.pdf"], limit=500)
        files = [p for p in files if "gseapy" in p.name.lower() or "summary" in p.name.lower() or p.parent == gseapy]
        selected = select_file("GSEA summary plot", files, key_prefix+"_gsea_summary", root)
        if selected:
            render_file_inline(selected, key_prefix+"_gsea_summary_view", height=900)
    with tabs[2]:
        files = matching_files(gseapy, ["*.png", "*.pdf"], limit=1000)
        selected = select_file("Pathway plot", files, key_prefix+"_gsea_pathway", root)
        if selected:
            render_file_inline(selected, key_prefix+"_gsea_pathway_view", height=950)


def render_file_browser(project: dict, key_prefix: str) -> None:
    root = data_dir(project)
    files = result_files(project)
    selected = select_file("Result file", files, key_prefix+"_all_files", root)
    if selected:
        st.caption(str(selected))
        render_file_inline(selected, key_prefix+"_all_files_view", height=900)


def integrated_results_explorer(project: dict, key_prefix: str = "results_explorer", compact: bool = False) -> None:
    root = data_dir(project)
    if not compact:
        st.subheader("RNA-seq Results Explorer")
    st.caption("Integrated Streamlit viewer. No separate Shiny app, no second port.")
    st.code(str(root))
    if not root.exists():
        st.warning("The project data folder was not found. Check the project config, results root, or visualizer data folder.")
        return
    tabs = st.tabs(["Overview", "QC", "Alignment QC", "Counts", "DESeq2", "GSEA", "Files"])
    with tabs[0]:
        render_overview(project, key_prefix+"_overview")
    with tabs[1]:
        render_qc_explorer(project, key_prefix+"_qc")
    with tabs[2]:
        render_alignment_explorer(project, key_prefix+"_alignment")
    with tabs[3]:
        render_counts_explorer(project, key_prefix+"_counts")
    with tabs[4]:
        render_deseq2_explorer(project, key_prefix+"_deseq2")
    with tabs[5]:
        render_gsea_explorer(project, key_prefix+"_gsea")
    with tabs[6]:
        render_file_browser(project, key_prefix+"_files")


def table_preview(path: Path, label: str, preview_rows: int = 500) -> None:
    if not path.exists():
        return
    size = path.stat().st_size
    st.markdown(f"**{label}**")
    st.caption(f"{path.name} · {format_size(size)} · previewing up to {preview_rows} rows")
    try:
        df = pd.read_csv(path, sep=None, engine="python", nrows=preview_rows)
        st.dataframe(df, use_container_width=True, height=360)
        st.download_button(
            f"Download {path.name}",
            data=path.read_bytes(),
            file_name=path.name,
            key="download_table_"+str(abs(hash(str(path)))),
        )
    except Exception:
        preview = path.read_text(errors="replace")[:20000]
        st.text_area(path.name, value=preview, height=260)


def result_files(project: dict, limit: int = 1000) -> List[Path]:
    root = data_dir(project)
    if not root.exists():
        return []
    keep = []
    for path in root.rglob("*"):
        if path.is_file() and path.suffix.lower() in [".txt", ".csv", ".tsv", ".html", ".png", ".pdf"]:
            keep.append(path)
            if len(keep) >= limit:
                break
    return sorted(keep)


def results_tab(project: dict) -> None:
    st.subheader("Outputs")
    col1, col2 = st.columns([3, 1])
    with col1:
        st.code(str(data_dir(project)))
    with col2:
        if st.button("Re-detect project paths"):
            project = apply_project_inference(project)
            save_project(project)
            st.success("Paths refreshed.")
    with st.expander("Detected project state", expanded=False):
        render_step_status(project)

    if project.get("analysis_type") == "RNA-seq":
        integrated_results_explorer(project, key_prefix="outputs_tab_"+project_id(project))
    else:
        render_file_browser(project, key_prefix="outputs_tab_"+project_id(project))


def main() -> None:
    st.set_page_config(page_title="CodeSpringLab Control Center", layout="wide")
    style()
    project = sidebar_project_selector()
    header(project)
    if not project:
        st.info("Create a project from the sidebar or select an existing one.")
        return

    tabs = st.tabs(["Setup", "Progress", "Design Matrix", "Run Pipeline", "Jobs & Logs", "Outputs"])
    with tabs[0]:
        project = setup_tab(project)
    with tabs[1]:
        progress_tab(project)
    with tabs[2]:
        design_tab(project)
    with tabs[3]:
        run_tab(project)
    with tabs[4]:
        jobs_tab(project)
    with tabs[5]:
        results_tab(project)


if __name__ == "__main__":
    main()
