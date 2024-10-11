#!/bin/bash

echo "Running python setup.py sdist bdist_wheel..."
python setup.py sdist bdist_wheel

if [ $? -ne 0 ]; then
    echo "Failed to run python setup.py sdist bdist_wheel"
    exit 1
fi

echo "Installing .whl files from ./dist..."
for whl_file in ./dist/*.whl; do
    if [ -f "$whl_file" ]; then
        echo "Installing $whl_file..."
        pip install "$whl_file"

        if [ $? -ne 0 ]; then
            echo "Failed to install $whl_file"
            exit 1
        fi
    else
        echo "No .whl files found in ./dist"
    fi
done

echo "All .whl files installed successfully."
