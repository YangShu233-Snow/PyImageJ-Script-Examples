import imagej
import os
import glob
import pandas as pd
from skimage import filters, measure

FIJI_PATH = "~/Fiji/Fiji.app"

def main():
    print("Initializing ImageJ2 (Headless)...")
    ij = imagej.init(FIJI_PATH, mode=imagej.Mode.HEADLESS)
    print(f"ImageJ2 initialized. Version: {ij.getVersion()}")

if __name__ == "__main__":
    main()