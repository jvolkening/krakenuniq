#!/bin/bash

# Portions (c) 2017-2018, Florian Breitwieser <fbreitwieser@jhu.edu>
# Copyright 2013-2015, Derrick Wood <dwood@cs.jhu.edu>
#
# Kraken is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Kraken is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Kraken.  If not, see <http://www.gnu.org/licenses/>.

set -e

DIR=$(dirname $0)
VERSION=`cat $(dirname $0)/VERSION`
INSTALL_JELLYFISH=1
MAKE_ARGS=()
MAKE_CLEAN="clean"
ADD_DEBUG_INFO=0
STATIC_LINK=0
LINK_DIR=

USAGE="Usage: $(basename $0) [OPTIONS] INSTALL_DIR

OPTIONS:
    -l BIN_DIR  Link KrakenUniq executables to BIN_DIR, e.g /usr/local/bin or ~/bin.
    -s          Skip building the bundled Jellyfish fork (installed as 'krakenfish').
    -c BIN      Use compiler BIN instead of g++.
    -S          Statically link compiled executables (requires static system libraries).
    -g          Add debug info
    -h          This help message

On MacOS, if you experience the error \"clang: fatal error: unsupported option '-fopenmp'\" on OSX, try installing g++ with brew, and using the option \"-c g++-7\".
"


while getopts "ChsSc:gl:" OPTION; do
    case $OPTION in
    c) MAKE_ARGS+=("CXX=$OPTARG") ;;
    C) MAKE_CLEAN="" ;;
    g) ADD_DEBUG_INFO=1 ;;
    s) INSTALL_JELLYFISH=0 ;;
    S) STATIC_LINK=1 ;;
    l) LINK_DIR="$OPTARG" ;;
    h) echo "$USAGE"; exit 0 ;;
    *) echo "Incorrect options provided. $USAGE"
       exit 1 ;;
    esac
done
shift $((OPTIND -1))

if [ -z "$1" ]
then
  echo "$USAGE"
  exit 0
fi

if [ "$1" = "KRAKEN_DIR" ]
then
  echo "Please replace \"KRAKEN_DIR\" with the name of the directory"
  echo "that you want to install Kraken in."
  exit 1
fi


# Perl cmd used to canonicalize dirname - "readlink -f" doesn't work
# on OS X.
export KRAKEN_DIR=$(perl -MCwd=abs_path -le 'print abs_path(shift)' "$1")

mkdir -p "$KRAKEN_DIR"
if [ "$INSTALL_JELLYFISH" == "1" ]; then
  if [[ ! -x "$KRAKEN_DIR/krakenfish" || "$STATIC_LINK" == 1 ]]; then
    echo "Installing krakenfish (bundled Jellyfish fork)"
    JF_SRC=$(perl -MCwd=abs_path -le 'print abs_path(shift)' "$DIR/external/jellyfish")
    JF_BUILD=$(mktemp -d)
    JF_PREFIX="$JF_BUILD/prefix"
    JF_CONFIGURE_ARGS=()
    JF_LDFLAGS=${LDFLAGS:-}
    if [[ "$STATIC_LINK" == 1 ]]; then
      JF_CONFIGURE_ARGS+=(--enable-all-static)
      JF_LDFLAGS="-static${JF_LDFLAGS:+ $JF_LDFLAGS}"
    fi
    cp -a "$JF_SRC/." "$JF_BUILD/src" && rm -rf "$JF_BUILD/src/autom4te.cache"
    (cd "$JF_BUILD/src" \
      && { [ -x ./configure ] || autoreconf -i; } \
      && LDFLAGS="$JF_LDFLAGS" ./configure --prefix="$JF_PREFIX" --disable-shared "${JF_CONFIGURE_ARGS[@]}" CXXFLAGS="-O2 -std=gnu++11 -Wno-error" --program-transform-name='s/^jellyfish$/krakenfish/' \
      && make -j && make install \
      && { if [ -x "$JF_PREFIX/bin/krakenfish" ]; then cp "$JF_PREFIX/bin/krakenfish" "$KRAKEN_DIR/krakenfish";
           else cp "$JF_PREFIX/bin/jellyfish" "$KRAKEN_DIR/krakenfish"; fi; } \
      && chmod +x "$KRAKEN_DIR/krakenfish"
    ) || {
      rm -rf "$JF_BUILD"
      if [[ "$STATIC_LINK" == 1 ]]; then
        echo "Static krakenfish build failed; check that static system libraries are installed." >&2
        exit 1
      fi
      echo "krakenfish installation failed, building new database may not work properly"
    }
    rm -rf "$JF_BUILD"
  fi
fi

[[ "$ADD_DEBUG_INFO" == 1 ]] && MAKE_ARGS+=("NDEBUG=-g")
if [[ "$STATIC_LINK" == 1 ]]; then
  MAKE_ARGS+=("LDFLAGS=-static${LDFLAGS:+ $LDFLAGS}")
fi
MAKE_GOALS=(install)
[[ -n "$MAKE_CLEAN" ]] && MAKE_GOALS=(clean install)
echo make -C "$DIR/src" "${MAKE_GOALS[@]}" "${MAKE_ARGS[@]}"
make -C "$DIR/src" "${MAKE_GOALS[@]}" "${MAKE_ARGS[@]}" || { echo "Error building KrakenUniq. See $(basename "$0") -h for options." >&2; exit 1; }
for file in $DIR/scripts/*
do
  [[ -f $file ]] || continue;
  perl -pl -e 'BEGIN { while (@ARGV) { $_ = shift; ($k,$v) = split /=/, $_, 2; $H{$k} = $v } }'\
           -e 's/#####=(\w+)=#####/$H{$1}/g' \
           "KRAKEN_DIR=$KRAKEN_DIR" "VERSION=$VERSION" \
           < "$file" > "$KRAKEN_DIR/$(basename $file)"
  if [ -x "$file" ]
  then
    chmod +x "$KRAKEN_DIR/$(basename $file)"
  fi
done
mkdir -p $KRAKEN_DIR/File
cp -r $DIR/scripts/File/* $KRAKEN_DIR/File

echo "Installed KrakenUniq in $KRAKEN_DIR."

if [[ "$LINK_DIR" != "" ]]; then
    [[ -d "$LINK_DIR" ]] || mkdir -p $LINK_DIR
    for file in $KRAKEN_DIR/krakenuniq*; do
        [ -x "$file" ] && ln -sf $file $LINK_DIR/`basename $file`
    done
    echo "Linked KrakenUniq files to $LINK_DIR."
else
    echo "You may want to link KrakenUniq scripts to /usr/bin or another directory for executables, or add the install directory to the PATH environment variable. Main executables: "

    for file in $KRAKEN_DIR/krakenuniq*; do
        [ -x "$file" ] && echo $file 
    done
fi

echo "KrakenUniq installation complete."
exit 0
