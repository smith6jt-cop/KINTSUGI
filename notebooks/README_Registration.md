# Image Registration and Merging Workflow

This directory contains tools for registering and merging PANCREAS and THYMUS images using the KINTSUGI pipeline.

## Files

- `Image_Registration_Workflow.ipynb` - Interactive Jupyter notebook for the registration workflow
- `registration_workflow.py` - Standalone Python script for batch processing
- `config_example.json` - Example configuration file
- `sample_channel_names.txt` - Example channel names file with 62 channels

## Quick Start

### Method 1: Using the Jupyter Notebook

1. Open `Image_Registration_Workflow.ipynb` in Jupyter/VS Code
2. Ensure your image files are in the `data/` directory
3. Run all cells to process both PANCREAS and THYMUS images

### Method 2: Using the Python Script

```bash
# Process both datasets
python registration_workflow.py

# Process only THYMUS
python registration_workflow.py --dataset thymus

# Process only PANCREAS  
python registration_workflow.py --dataset pancreas

# Use custom configuration
python registration_workflow.py --config config_example.json

# Use custom directories
python registration_workflow.py --data-dir /path/to/data --output-dir /path/to/outputs
```

## Required Files

### For PANCREAS Registration
Place in `data/` directory:
- At least 2 TIFF files with "pancreas" in the filename (fixed and moving images)
- Optional: `pancreas_channel_names.txt` with channel names

### For THYMUS Registration  
Place in `data/` directory:
- `ThymusDAPI_Xenium_downsampled2x_EDF.ome.tiff` (fixed image)
- `THYMUS.ome.tiff` (moving image)
- Optional: `thymus_channel_names.txt` with channel names

## Output

The workflow generates:
- `PANCREAS_registered_62ch_ImageJ.tiff` - Registered PANCREAS image with 62 channels
- `THYMUS_registered_62ch_ImageJ.tiff` - Registered THYMUS image with 62 channels

Both files are:
- Pyramidal TIFF format for efficient viewing
- OME metadata compliant
- Compatible with QuPath
- LZW compressed
- Exactly 62 channels (padded or cropped as needed)

## Workflow Steps

1. **Load channel names** from text files or use defaults
2. **Initialize VALIS registration** with optimized parameters
3. **Perform rigid and non-rigid registration** between images
4. **Warp and merge slides** with proper channel mapping
5. **Process arrays** to ensure exactly 62 channels
6. **Create OME-XML metadata** for QuPath compatibility
7. **Save as pyramidal TIFF** with compression

## Configuration

The `config_example.json` file shows the structure for custom configurations:

```json
{
  "pancreas": {
    "name": "PANCREAS",
    "fixed_image": "data/pancreas_fixed.ome.tiff",
    "moving_image": "data/pancreas_moving.ome.tiff",
    "channel_names_file": "data/pancreas_channel_names.txt",
    "output_file": "outputs/PANCREAS_registered_62ch_ImageJ.tiff"
  },
  "thymus": {
    "name": "THYMUS", 
    "fixed_image": "data/ThymusDAPI_Xenium_downsampled2x_EDF.ome.tiff",
    "moving_image": "data/THYMUS.ome.tiff",
    "channel_names_file": "data/thymus_channel_names.txt",
    "output_file": "outputs/THYMUS_registered_62ch_ImageJ.tiff"
  }
}
```

## Error Handling

The workflow includes comprehensive error handling:
- File existence checks
- Safe array processing
- Automatic cleanup of temporary files
- JVM cleanup after registration
- Graceful degradation with default values

## Dependencies

- KINTSUGI environment with all required packages
- Kreg module (included in KINTSUGI)
- VALIS registration library
- NumPy, tifffile, PIL

## Troubleshooting

1. **Import errors**: Ensure you're running in the KINTSUGI conda environment
2. **File not found**: Check paths in configuration and ensure files exist
3. **Memory errors**: Adjust `max_processed_image_dim_px` parameter for large images
4. **Registration failures**: Try different parameter values or check image quality
5. **JVM errors**: Restart Python session if Java Virtual Machine issues occur

## Contact

For issues or questions, refer to the main KINTSUGI documentation or repository.