version=1.1.7

wget https://github.com/wulongict/SpectralArchive/releases/download/v${version}/Spectroscape_CPU-${version}.deb

# install
# remove
sudo apt remove -y spectroscape_cpu
sudo apt install -y ./Spectroscape_CPU-${version}.deb
rm ./Spectroscape_CPU-${version}.deb

# test.
(cd ~/data/spectroscape/spectral_archives && spectroscape --add --datasearchpath ../new_data)

# remove package.
