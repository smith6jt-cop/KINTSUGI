
//This is needed for the segmentation notebook that uses the Mesmer package.  Also reduces the file size, but also impacts the pixel range.
run("8-bit");

// Check for border pixel distortion, if crop is used, all images must be cropped to the same size.
makeRectangle(40, 48, 9880, 9392);

run("Crop");

// Another method of fixing border pixel distortion without changing image size is to turn the pixels black.
run("Make Inverse");
run("Set...", "value=0");
// This collapses the histogram to where the bulk of the pixels are for visualization only, to make permanent (alter pixels) use run("Apply LUT");
run("Enhance Contrast", "saturated=0.35");

// This will typically smooth any remaining bright spots caused by AF removal.
run("Median...", "radius=2 slice");
run("Remove Outliers...", "radius=2 threshold=50 which=Bright slice");

// This will remove any "haze"
run("Subtract...", "value=2000 slice");

// If there are still bright/dim regions in the image that do not represent biological variation, CLAHE will adjust contrast based on small regional samples determined by "blocksize." This alters the relative pixel values.
run("Enhance Local Contrast (CLAHE)", "blocksize=65 histogram=125 maximum=5 mask=*None*");

// Run this to alter the pixel values to a compressed range to increase contrast.
run("Apply LUT");

