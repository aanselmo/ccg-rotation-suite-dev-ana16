#!/usr/bin/env python3
"""
export_nexus_ogt.py - run `cnvkit.py export nexus-ogt` for every sample in a
sarek project, then plot each result with plot_ogt_v3.py. Both steps run
inside a Docker image by default (see the Dockerfile); --no-docker runs the
tools installed on this machine instead.

For a project name given on the command line, searches (including their
subfolders) only
    <PROJECTS_ROOT>/<project>/variant_calling/cnvkit/   for
        {sample}.cnr                              (export)
        {sample}.call.cns                         (plot segments)
    <PROJECTS_ROOT>/<project>/annotation/mutect2/       for
        {sample}.mutect2.filtered_VEP.ann.vcf.gz  (export)
checks that every sample has all three, warns about any gaps, and writes
    <PROJECTS_ROOT>/<project>/nexus-ogt-files/{sample}.nexus-ogt
    <PROJECTS_ROOT>/<project>/cnv_plots/{sample}.nexus-ogt.png

A sample missing its .cnr or VCF is skipped. A sample missing only its
.call.cns is still exported and plotted without segments (with a warning).
Use --strict to skip any sample that is missing any of the three.
Every sample with a .nexus-ogt file (new or already present) is plotted.

Docker: each command runs in `docker run --rm` as your own user and groups,
so output files belong to you. The project directory is mounted at the same
path inside the container (plus the real location of any input file that is
a symlink pointing outside it), so all paths are identical inside and out.

Examples
    python export_nexus_ogt.py SO31252a
    python export_nexus_ogt.py SO31252a --dry-run
    python export_nexus_ogt.py SO31252a --image cnv-ogt:0.9.14 --force   # another build
    python export_nexus_ogt.py SO31252a --cytoband /ref/hg38_cytoBand.txt.gz
    python export_nexus_ogt.py SO31252a --no-docker
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Hardcoded settings - edit these for your system
# ---------------------------------------------------------------------------
PROJECTS_ROOT = Path("/mnt/storage/centers/ccgd/Illumina/pipeline_sarek")
OGT_SUBDIR = "nexus-ogt-files"
PLOT_SUBDIR = "cnv_plots"
GAIN_CUTOFF = 0.2
LOSS_CUTOFF = -0.25

# Docker mode (default): tool paths are paths INSIDE the image
DOCKER_BIN = "docker"
DOCKER_IMAGE = "us-central1-docker.pkg.dev/ccgd-174313/cnv-ogt:0.9.13"
IMAGE_CNVKIT = "cnvkit.py"
IMAGE_PYTHON = "python"
IMAGE_PLOT_SCRIPT = "/opt/cnv/plot_ogt_v3.py"
IMAGE_CYTOBAND = "/opt/cnv/cytoBand.txt.gz"

# --no-docker mode: tool paths on this machine
LOCAL_CNVKIT = "/home/anselmo/tonyccg/environments/singlecell/bin/cnvkit.py"
LOCAL_PLOT_SCRIPT = Path(__file__).resolve().parent / "plot_ogt_v3.py"
LOCAL_CYTOBAND = Path(__file__).resolve().parent / "cytoBand.txt.gz"

# File type -> filename suffix
SUFFIXES = {
    "cnr": ".cnr",
    "call": ".call.cns",
    "vcf": ".mutect2.filtered_VEP.ann.vcf.gz",
}
REQUIRED = tuple(SUFFIXES)
# File type -> directory searched (relative to the project; subfolders included)
SEARCH_DIRS = {
    "cnr": "variant_calling/cnvkit",
    "call": "variant_calling/cnvkit",
    "vcf": "annotation/mutect2",
}
# Files from sarek's germline / tumor-only CNVkit runs carry one of these tags
# before the suffix (e.g. S1.germline.call.cns, S1.tumor_only.call.cns) and
# are ignored; only the tumor/normal (somatic) files are used.
EXCLUDE_TAGS = ("germline", "tumor_only")


def warn(msg):
    print(f"WARNING: {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Finding and checking input files
# ---------------------------------------------------------------------------
def gather(project_dir):
    """Return {file_type: {sample: Path}} for all file types, searching only
    each type's directory in SEARCH_DIRS."""
    found = {k: {} for k in SUFFIXES}
    dups = []
    excluded = 0
    for d in sorted(set(SEARCH_DIRS.values())):
        if not (project_dir / d).is_dir():
            warn(f"search directory not found: {project_dir / d}")
    for ftype, suffix in SUFFIXES.items():
        search_dir = project_dir / SEARCH_DIRS[ftype]
        if not search_dir.is_dir():
            continue
        for p in sorted(search_dir.rglob(f"*{suffix}")):
            if not p.is_file():
                continue
            sample = p.name[: -len(suffix)]
            if not sample:
                continue
            if any(tag in sample.split(".")[1:] for tag in EXCLUDE_TAGS):
                excluded += 1
                continue
            if sample in found[ftype]:
                dups.append((ftype, sample, found[ftype][sample], p))
                continue
            found[ftype][sample] = p
    for ftype, sample, kept, extra in dups:
        warn(f"duplicate {ftype} file for {sample}: using {kept}, ignoring {extra}")
    if excluded:
        print(f"Ignored {excluded} {'/'.join(EXCLUDE_TAGS)} file(s).")
    return found


def check(found):
    """Print a completeness report; return (all samples, {sample: [missing]})."""
    samples = sorted(set().union(*(found[t].keys() for t in REQUIRED)))
    missing = {s: [t for t in REQUIRED if s not in found[t]] for s in samples}
    missing = {s: m for s, m in missing.items() if m}

    print(f"Found {len(samples)} sample(s): "
          f"{len(found['cnr'])} .cnr, {len(found['call'])} .call.cns, "
          f"{len(found['vcf'])} VEP VCF")
    if not samples:
        return samples, missing
    if missing:
        warn(f"{len(missing)} sample(s) do not have all three files:")
        w = max(len(s) for s in missing)
        for s, m in missing.items():
            print(f"    {s:<{w}}  missing: {', '.join(SUFFIXES[t] for t in m)}",
                  file=sys.stderr)
    else:
        print("All samples have .cnr, .call.cns and VEP VCF.")
    return samples, missing


# ---------------------------------------------------------------------------
# Running commands, natively or in Docker
# ---------------------------------------------------------------------------
class Runner:
    """Builds and runs commands either directly or via `docker run`."""

    def __init__(self, image=None, docker_bin=DOCKER_BIN):
        self.image = image
        self.docker_bin = docker_bin
        self.rw = set()   # directories mounted read-write
        self.ro = set()   # directories mounted read-only

    def mount(self, path, writable=False):
        """Make a host directory (or a file's directory) visible in the
        container at the same path."""
        if not self.image:
            return
        p = Path(path).resolve()
        d = p if p.is_dir() else p.parent
        (self.rw if writable else self.ro).add(d)

    def _mounts(self):
        """Mount list with nested directories dropped (a read-write parent
        covers its children; a read-only child under a read-write parent is
        redundant)."""
        rw = sorted(self.rw)
        rw = [d for d in rw if not any(o != d and o in d.parents for o in rw)]
        ro = sorted(d for d in self.ro
                    if not any(o == d or o in d.parents for o in rw))
        ro = [d for d in ro if not any(o != d and o in d.parents for o in ro)]
        return [(d, "rw") for d in rw] + [(d, "ro") for d in ro]

    def wrap(self, cmd):
        if not self.image:
            return cmd
        dk = [self.docker_bin, "run", "--rm",
              "--user", f"{os.getuid()}:{os.getgid()}"]
        for g in sorted(set(os.getgroups()) - {os.getgid()}):
            dk += ["--group-add", str(g)]   # keep group write access
        dk += ["-e", "HOME=/tmp", "-e", "MPLCONFIGDIR=/tmp/matplotlib",
               "-w", "/tmp"]
        for d, mode in self._mounts():
            dk += ["-v", f"{d}:{d}:{mode}"]
        return dk + [self.image] + cmd

    def run(self, cmd, dry_run):
        """Run a command; return True on success. Prints the tail of the
        error output on failure."""
        full = self.wrap(cmd)
        if dry_run:
            print("   " + " ".join(full))
            return True
        res = subprocess.run(full, capture_output=True, text=True)
        if res.returncode != 0:
            print(f"   FAILED (exit {res.returncode})", file=sys.stderr)
            for line in (res.stderr or res.stdout).strip().splitlines()[-5:]:
                print(f"   | {line}", file=sys.stderr)
            return False
        return True


def check_docker(docker_bin, image):
    if not shutil.which(docker_bin):
        sys.exit(f"'{docker_bin}' not found on PATH (or use --no-docker)")
    res = subprocess.run([docker_bin, "image", "inspect", image],
                         capture_output=True, text=True)
    if res.returncode != 0:
        msg = (res.stderr or "").strip().splitlines()
        msg = msg[-1] if msg else ""
        if "permission denied" in msg.lower():
            sys.exit(f"Can't talk to Docker: {msg}\n"
                     "Your user may need to be in the 'docker' group.")
        sys.exit(f"Docker image '{image}' not available: {msg}\n"
                 f"Build it with:  docker build -t {image} .")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project", help="project name, e.g. SO31252a")
    ap.add_argument("--image", default=DOCKER_IMAGE,
                    help=f"Docker image with cnvkit + plot script [{DOCKER_IMAGE}]")
    ap.add_argument("--docker-bin", default=DOCKER_BIN,
                    help="docker-compatible CLI, e.g. podman [docker]")
    ap.add_argument("--no-docker", action="store_true",
                    help="run cnvkit and the plot script installed on this machine")
    ap.add_argument("--cnvkit",
                    help=f"cnvkit.py [image: {IMAGE_CNVKIT}; local: {LOCAL_CNVKIT}]")
    ap.add_argument("--plot-script",
                    help=f"plot_ogt_v3.py on this machine; in Docker mode it is mounted "
                         f"and used instead of the copy in the image "
                         f"[image: {IMAGE_PLOT_SCRIPT}; local: {LOCAL_PLOT_SCRIPT}]")
    ap.add_argument("--cytoband",
                    help=f"UCSC cytoBand.txt(.gz) on this machine; in Docker mode it is "
                         f"mounted and used instead of the copy in the image "
                         f"[image: {IMAGE_CYTOBAND}; local: {LOCAL_CYTOBAND}]")
    ap.add_argument("--gain-cutoff", type=float, default=GAIN_CUTOFF,
                    help=f"gain cutoff for plots [{GAIN_CUTOFF}]")
    ap.add_argument("--loss-cutoff", type=float, default=LOSS_CUTOFF,
                    help=f"loss cutoff for plots [{LOSS_CUTOFF}]")
    ap.add_argument("--strict", action="store_true",
                    help="skip samples missing any of the three files")
    ap.add_argument("--force", action="store_true",
                    help="overwrite existing .nexus-ogt and .png files")
    ap.add_argument("--no-plot", action="store_true", help="export only, no plots")
    ap.add_argument("--dry-run", action="store_true",
                    help="check files and print commands without running them")
    a = ap.parse_args()

    use_docker = not a.no_docker
    do_plot = not a.no_plot

    project_dir = (PROJECTS_ROOT / a.project).resolve()
    if not project_dir.is_dir():
        sys.exit(f"Project directory not found: {project_dir}")
    ogt_dir = project_dir / OGT_SUBDIR
    plot_dir = project_dir / PLOT_SUBDIR

    # ---- Tools: resolve and check up front --------------------------------
    runner = Runner(a.image if use_docker else None, a.docker_bin)
    if use_docker:
        if not a.dry_run:
            check_docker(a.docker_bin, a.image)
        cnvkit = a.cnvkit or IMAGE_CNVKIT
        python = IMAGE_PYTHON
        # Host files given on the command line are mounted; defaults live in the image
        plot_script, cytoband = IMAGE_PLOT_SCRIPT, IMAGE_CYTOBAND
        for opt, val in (("--plot-script", a.plot_script), ("--cytoband", a.cytoband)):
            if val:
                if not Path(val).is_file():
                    sys.exit(f"{opt} not found: {val}")
                runner.mount(val)
        if a.plot_script:
            plot_script = str(Path(a.plot_script).resolve())
        if a.cytoband:
            cytoband = str(Path(a.cytoband).resolve())
        runner.mount(project_dir, writable=True)
    else:
        c = a.cnvkit or LOCAL_CNVKIT
        cnvkit = c if Path(c).is_file() else shutil.which(c)
        if not cnvkit and not a.dry_run:
            sys.exit(f"cnvkit.py not found: {c} (use --cnvkit)")
        cnvkit = cnvkit or c
        python = sys.executable
        plot_script = str(a.plot_script or LOCAL_PLOT_SCRIPT)
        cytoband = str(a.cytoband or LOCAL_CYTOBAND)
        if do_plot and not a.dry_run:
            if not Path(plot_script).is_file():
                sys.exit(f"plot script not found: {plot_script} (use --plot-script)")
            if not Path(cytoband).is_file():
                sys.exit(f"cytoband file not found: {cytoband} (use --cytoband)")

    # ---- Inputs -----------------------------------------------------------
    print(f"Searching {project_dir} in: "
          + ", ".join(sorted(set(SEARCH_DIRS.values()))))
    found = gather(project_dir)
    samples, missing = check(found)
    if not samples:
        sys.exit("No .cnr, .call.cns or VEP VCF files found.")

    # sarek output can be symlinks into a work directory elsewhere; use real
    # paths and mount their directories so the container can read them
    for ftype in found:
        for s, p in found[ftype].items():
            found[ftype][s] = p.resolve()
            runner.mount(found[ftype][s])

    if not a.dry_run:
        ogt_dir.mkdir(parents=True, exist_ok=True)
        if do_plot:
            plot_dir.mkdir(parents=True, exist_ok=True)

    print(f"Running with: {'Docker image ' + a.image if use_docker else 'local tools'}")

    # ---- Step 1: export ---------------------------------------------------
    print("\n== Exporting .nexus-ogt ==")
    exp_ok, exp_present, exp_missing, exp_fail = [], [], [], []
    have_ogt = []   # samples with a .nexus-ogt to plot (new or existing)
    for s in samples:
        m = missing.get(s, [])
        out = ogt_dir / f"{s}.nexus-ogt"
        if "cnr" in m or "vcf" in m or (a.strict and m):
            print(f">> {s}  SKIP (missing {', '.join(SUFFIXES[t] for t in m)})")
            exp_missing.append(s)
            continue
        if out.exists() and not a.force:
            print(f">> {s}  exists, not re-exported (use --force to overwrite)")
            exp_present.append(s)
            have_ogt.append(s)
            continue

        print(f">> {s}")
        cmd = [cnvkit, "export", "nexus-ogt",
               str(found["cnr"][s]), str(found["vcf"][s]), "-o", str(out)]
        if runner.run(cmd, a.dry_run) and (a.dry_run or out.exists()):
            if not a.dry_run:
                print(f"   OK -> {out}")
            exp_ok.append(s)
            have_ogt.append(s)
        else:
            exp_fail.append(s)

    # ---- Step 2: plot -----------------------------------------------------
    plt_ok, plt_present, plt_fail = [], [], []
    if do_plot:
        print("\n== Plotting ==")
        for s in have_ogt:
            ogt = ogt_dir / f"{s}.nexus-ogt"
            png = plot_dir / f"{s}.nexus-ogt.png"
            if png.exists() and not a.force and s not in exp_ok:
                print(f">> {s}  exists, not re-plotted (use --force to overwrite)")
                plt_present.append(s)
                continue

            seg = found["call"].get(s)
            print(f">> {s}" + ("" if seg else "  (no .call.cns; plotting without segments)"))
            cmd = [python, plot_script, str(ogt), "--cytoband", cytoband]
            if seg:
                cmd += ["--segments", str(seg),
                        "--gain-cutoff", str(a.gain_cutoff),
                        "--loss-cutoff", str(a.loss_cutoff)]
            cmd += ["--out", str(png)]
            if runner.run(cmd, a.dry_run) and (a.dry_run or png.exists()):
                if not a.dry_run:
                    print(f"   OK -> {png}")
                plt_ok.append(s)
            else:
                plt_fail.append(s)

    # ---- Summary ----------------------------------------------------------
    print(f"\nExport: {len(exp_ok)} done, {len(exp_present)} already present, "
          f"{len(exp_missing)} missing inputs, {len(exp_fail)} failed.  "
          f"Output: {ogt_dir}")
    if do_plot:
        print(f"Plots:  {len(plt_ok)} done, {len(plt_present)} already present, "
              f"{len(plt_fail)} failed.  Output: {plot_dir}")
    if exp_missing:
        print("Missing inputs (not exported): " + ", ".join(exp_missing))
    if exp_fail:
        print("Export failed: " + ", ".join(exp_fail), file=sys.stderr)
    if plt_fail:
        print("Plot failed: " + ", ".join(plt_fail), file=sys.stderr)
    if exp_fail or plt_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
