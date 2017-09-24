This a minimal project using proxTV as a third-party.
To compile it, it requires that proxTV is installed somewhere in your system.
The installation folder can be modified when configuring proxTV with the option

    -DCMAKE_INSTALL_PREFIX:PATH="/proxTV_install_folder"

The default is usually the system folder `/usr`. Please read https://github.com/albarji/proxTV/README.md#c-interface for more info.

Once installed in your system, if it is not in your system path, provide the folder where the `proxTVConfig.cmake` file is located with the option `-proxTV_DIR`:

    mkdir ~/dummy_project ; cd ~/dummy_project
    cmake /path/proxTV/test/project_using_proxTV/ -DproxTV_DIR:PATH="/proxTV_install_folder/lib/cmake/proxTV"
