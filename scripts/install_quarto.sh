#!/usr/bin/env bash
set -euo pipefail

version="${QUARTO_VERSION:-1.11.4}"
install_root="${QUARTO_INSTALL_ROOT:-$HOME/.local/share/quarto}"
bin_dir="${QUARTO_BIN_DIR:-$HOME/.local/bin}"

case "$(uname -m)" in
  x86_64) architecture="amd64" ;;
  aarch64|arm64) architecture="arm64" ;;
  *)
    echo "Unsupported architecture: $(uname -m)" >&2
    exit 1
    ;;
esac

target="$install_root/$version"
launcher="$bin_dir/quarto"
if [[ -e "$launcher" && ! -L "$launcher" ]]; then
  echo "The path $launcher exists and is not a symbolic link." >&2
  exit 1
fi

temporary_dir="$(mktemp -d)"
trap 'rm -r "$temporary_dir"' EXIT
archive="$temporary_dir/quarto.tar.gz"
checksums="$temporary_dir/checksums.txt"
base_url="https://github.com/quarto-dev/quarto-cli/releases/download/v${version}"

curl -fsSL -o "$archive" "$base_url/quarto-${version}-linux-${architecture}.tar.gz"
curl -fsSL -o "$checksums" "$base_url/quarto-${version}-checksums.txt"
expected="$(awk -v file="quarto-${version}-linux-${architecture}.tar.gz" '$2 == file {print $1}' "$checksums")"
actual="$(sha256sum "$archive" | awk '{print $1}')"
if [[ -z "$expected" || "$actual" != "$expected" ]]; then
  echo "The Quarto archive checksum does not match." >&2
  exit 1
fi

mkdir -p "$target" "$bin_dir"
tar -xzf "$archive" -C "$target" --strip-components=1
ln -sfn "$target/bin/quarto" "$launcher"
"$launcher" --version

if [[ ":$PATH:" != *":$bin_dir:"* ]]; then
  echo "Add $bin_dir to PATH before you start the notebook kernel."
fi
