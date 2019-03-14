

pkgbuild --version 1.6.2.Alpha --identifier com.cbica.captk --install-location /Applications --component ${TRAVIS_BUILD_DIR}/bin/_CPack_Packages/OSX/DragNDrop/CaPTk_1.6.2.Alpha/CaPTk_1.6.2.Alpha.app/  ${TRAVIS_BUILD_DIR}/bin/CaPTk_1.6.2.Alpha.pkg
productbuild --synthesize --package CaPTk_1.6.2.Alpha.pkg ${TRAVIS_BUILD_DIR}/bin/distribution.xml

echo '<?xml version="1.0" encoding="utf-8"?>
<installer-gui-script minSpecVersion="1">
    <title>CaPTk_1.6.2.Alpha</title>
    <license file="Combined.txt"></license>
    <pkg-ref id="com.cbica.captk"/>
    <options customize="never" require-scripts="false"/>
    <choices-outline>
        <line choice="default">
            <line choice="com.cbica.captk"/>
        </line>
    </choices-outline>
    <choice id="default"/>
    <choice id="com.cbica.captk" visible="false">
        <pkg-ref id="com.cbica.captk"/>
    </choice>
    <pkg-ref id="com.cbica.captk" version="1.6.2.Alpha" onConclusion="none">CaPTk_1.6.2.Alpha.pkg</pkg-ref>
</installer-gui-script>' > "${TRAVIS_BUILD_DIR}/bin/distribution.xml"

productbuild --distribution ${TRAVIS_BUILD_DIR}/bin/distribution.xml --resources ${TRAVIS_BUILD_DIR}/bin/_CPack_Packages/OSX/DragNDrop/CaPTk_1.6.2.Alpha/CaPTk_1.6.2.Alpha.app/Contents/Resources/license/ --package-path ${TRAVIS_BUILD_DIR}/bin ${TRAVIS_BUILD_DIR}/bin/CaPTk_1.6.2.Alpha_Installer.pkg