# list all commands
default:
    @just --list

# generate fixture outputs by running pdbe-arpeggio against fixtures
fixtures:
    #!/usr/bin/env bash
    source oracle_env/bin/activate
    python tests/tools/generate_fixtures.py

#format files
fmt:
    cargo fmt

#run with no changes , just pass/fail
fmt-check:
    cargo fmt -- --check

#run clippy for linting
clippy:
    cargo clippy

#hard fail linting in ci 
clippy-ci:
    cargo clippy --all-targets -- -D warnings

# run tests
test:
    cargo test

# run local build
build:
    cargo build

# run release build
build-release:
    cargo build --release

#runs full local dev check
check: fmt clippy test

#runs ci check with hard failures
ci: fmt-check clippy-ci test

