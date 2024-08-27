#!/usr/bin/env python3
import sys, os

sys.path.append("./compiled_softwares/OrthoFinder_source/")

from scripts_of.__main__ import main

if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
