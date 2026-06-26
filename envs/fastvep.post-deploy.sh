#!/usr/bin/env bash
# Snakemake conda post-deploy hook for the fastvep env.
# Runs automatically after the env is created, with $CONDA_PREFIX set to the env.
# Builds + installs the fastvep Rust CLI into the env so the `fastvep` binary
# lands on the env PATH ($CONDA_PREFIX/bin/fastvep).
set -euo pipefail

# Pinned for reproducibility: Huang-lab/fastVEP, v0.2.0.
FASTVEP_GIT="https://github.com/Huang-lab/fastVEP.git"
FASTVEP_REV="785922ebcaacd3f646d5f1edf374f40f1a39efe5"

# Use the env's own cargo/rustc (installed via fastvep.yaml).
export PATH="${CONDA_PREFIX}/bin:${PATH}"

# Keep cargo's build/registry state inside the env so it is self-contained and
# does not collide with the user's ~/.cargo while building.
export CARGO_HOME="${CONDA_PREFIX}/.cargo"

# cargo install places the binary at $CONDA_PREFIX/bin/fastvep.
# Package = fastvep-cli (workspace member), binary = fastvep.
cargo install \
    --git "${FASTVEP_GIT}" \
    --rev "${FASTVEP_REV}" \
    --bin fastvep \
    --root "${CONDA_PREFIX}" \
    --locked \
    fastvep-cli

# Sanity check that the binary is on the env PATH.
"${CONDA_PREFIX}/bin/fastvep" --version
