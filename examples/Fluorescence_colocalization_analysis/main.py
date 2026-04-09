import imagej
import numpy as np
import pandas as pd
import os
import pandas as pd
from pathlib import Path
from skimage import filters, measure
import scyjava as sj
import xarray

FIJI_PATH = "sc.fiji:fiji:2.14.0"

def process_czi_imagej(ij: imagej.GatewayAddons, input_czi: Path|str, output_csv: Path|str)->bool:
    print("Setting Bio-Formats Options...")
    bf = sj.jimport('loci.plugins.BF')
    options = sj.jimport('loci.plugins.in.ImporterOptions')()
    options.setSplitChannels(True)
    options.setOpenAllSeries(True)
    options.setVirtual(True)
    options.setId(str(input_czi))

    print("Applying Bio-Formats Options...")
    imps = bf.openImagePlus(options)

    cal = imps[0].getCalibration()
    pixel_width = cal.pixelWidth
    pixel_height = cal.pixelHeight
    spatial_unit = cal.getUnit() # 通常是 'micron' 或 'µm'
    
    # 单个像素代表的真实物理面积
    pixel_physical_area = pixel_width * pixel_height
    print(f"Spatial Calibration: 1 pixel = {pixel_width:.4f} x {pixel_height:.4f} {spatial_unit}")

    print("Loading Image...")
    fluor_imp_xr: xarray.Dataset = ij.py.from_java(imps[0])
    fluor_imp = fluor_imp_xr.values

    dapi_imp_xr: xarray.Dataset = ij.py.from_java(imps[1])
    dapi_imp: np.ndarray = dapi_imp_xr.values

    print("Analysis Image...")
    # transfer to 8-bit (linear)
    dapi_min, dapi_max = (dapi_imp.min(), dapi_imp.max())
    dapi_img_8bit: np.ndarray = ((dapi_imp - dapi_min) / (dapi_max - dapi_min) * 255).astype(np.uint8)

    # thresholding by mean
    threshold_val: float = filters.threshold_mean(dapi_img_8bit)
    dapi_mask = dapi_img_8bit > threshold_val

    # as same as fiji analysis particle
    cell_labels = measure.label(dapi_mask)

    # Measure fluorsence
    cell_props = measure.regionprops(cell_labels, fluor_imp)

    # Background Brightness (5% in not covered area)
    background_mask = ~dapi_mask
    background_pixels = np.percentile(fluor_imp[background_mask], 10)

    bg_brighrness_mean = np.mean(background_pixels)
    bg_brighrness_median = np.median(background_pixels)
    bg_brighrness_std = np.std(background_pixels)
    
    # select cell area and calculate fluorsence intensity
    results = []
    for prop in cell_props:
        physical_area = prop.area * pixel_physical_area
        
        if physical_area >= 50:
            results.append({
                'Cell_ID': prop.label,
                'Area_Pixels': prop.area,
                f'Area_{spatial_unit}2': physical_area,
                'Mean_Fluor_Intensity': prop.intensity_mean,
                'Max_Fluor_Intensity': prop.intensity_max,
                'Integrated_Density': prop.intensity_mean * prop.area,
                'Integrated_Density': prop.intensity_mean * physical_area, 
                'Raw_Integrated_Density': prop.intensity_mean * prop.area,
                f'Centroid_X_{spatial_unit}': prop.centroid[1] * pixel_width,
                f'Centroid_Y_{spatial_unit}': prop.centroid[0] * pixel_height,
                'Centroid_X': prop.centroid[1],
                'Centroid_Y': prop.centroid[0],
                'Background_Mean': bg_brighrness_mean,
                'Background_Median': bg_brighrness_median,
                'Backgournd_Std': bg_brighrness_std,
                'Corrected_Fluorescence_Intensity': prop.intensity_mean - bg_brighrness_mean
            })

    # save data
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"Processed {os.path.basename(input_czi)} -> Found {len(df)} cells.")


def main():
    print("Initializing ImageJ2 (Headless)...")
    sj.config.add_options('-Xmx6g')
    ij: imagej.GatewayAddons = imagej.init(FIJI_PATH, mode=imagej.Mode.HEADLESS)
    print(f"ImageJ2 initialized. Version: {ij.getVersion()}")

    for path in Path('./data').glob("*"):
        targer_dir = Path(f"./output/{path.name}")
        Path.mkdir(targer_dir, exist_ok=True)

        for file_path in Path(path).glob("*"):
            print(file_path.name[:-4])
            process_czi_imagej(ij, file_path, f"{targer_dir}/{file_path.name[:-4]}.csv")


if __name__ == "__main__":
    main()
