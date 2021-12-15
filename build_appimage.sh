
# After a new release 
# Download the lastet release of IGoR and then
./configure --prefix=/usr --datarootdir=/usr/share/applications

# modify dataroot?
make DESTDIR=$(pwd)/AppDir install

mkdir -p AppDir/usr/share/icons
mkdir -p AppDir/usr/share/icons/hicolor
mkdir -p AppDir/usr/share/icons/hicolor/256x256
mkdir -p AppDir/usr/share/icons/hicolor/256x256/apps/

cp igoricon.png AppDir/usr/share/icons/hicolor/256x256/apps/
cp igor.desktop AppDir/usr/share/applications/


./linuxdeploy-x86_64.AppImage --appdir AppDir --output appimage
