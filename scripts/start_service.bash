#!/bin/bash


archive_path="$HOME/data/spectroscape/spectral_archives_latest"
bin_path="$HOME/code/SpectralArchive/build/bin/spectroscape"
wwwroot_path="`dirname $bin_path`/../www"
cd $archive_path && spawn-fcgi -p 8710 -n -- $bin_path --run --wwwroot $wwwroot_path