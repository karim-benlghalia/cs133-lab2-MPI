#!/bin/sh
if [ "$(uname -s)" == 'Darwin' ]; then
  xcode-select --install
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  /usr/local/bin/brew install llvm
else
  sudo apt update
  sudo env DEBIAN_FRONTEND=noninteractive apt install build-essential -y
fi
