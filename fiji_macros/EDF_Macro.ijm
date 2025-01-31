// (17, 9484, 9962) image 
// NVIDIA RTX A4500: 42.1s in pyimagej macro split 2x2
// Intel(R) UHD Graphics 770: 1m59.4s in pyimagej macro split 2x2


// run("CLIJ2 Macro Extensions", "cl_device=[13th Gen Intel(R) Core(TM) i9-13900] automatic_output_naming=false");
run("CLIJ2 Macro Extensions", "cl_device=[NVIDIA RTX A4500] automatic_output_naming=false");

Ext.CLIJ2_clear();
setBatchMode(true);
GPU=true

radius_x = 5.0;
radius_y = 5.0;
sigma = 20.0;
xTiles = 3;
yTiles = 3;
z_start = 3;
z_end = 15;

Stack.getDimensions(width, height, channels, slices, frames);

if (slices != z_end - z_start + 1) {
    run("Slice Keeper", "first=z_start last=z_end increment=1");
    input = "deconvolved kept stack";
    Stack.getDimensions(width, height, channels, slices, frames);
    } else {
    input = "deconvolved";
}

if (GPU == true) {
	output = "output";
	newImage(output, "16-bit black", width, height, 1);
    } else {
    Ext.CLIJ2_push(input);
	Ext.CLIJ2_create2D(output, width, height, 16);
}

numTilesZ = 1;
numTilesX = xTiles;
numTilesY = yTiles;

tileDepth = round((slices/numTilesZ));
tileWidth = round((width/numTilesX));
tileHeight = round((height/numTilesY));


for (x = 0; x < numTilesX; x++) {
	for (y = 0; y < numTilesY; y++) {
		for (z = 0; z < numTilesZ; z++) {

			Ext.CLIJ2_pushTile(input, x, y, z, tileWidth, tileHeight, tileDepth, 0, 0, 0);		
			Ext.CLIJ2_extendedDepthOfFocusVarianceProjection(input, output, radius_x, radius_y, sigma);
            Ext.CLIJ2_release(input);
            if (GPU == true) {
				Ext.CLIJ2_pullTile(output, x, y, z, tileWidth, tileHeight, 1, 0, 0, 0);
				Ext.CLIJ2_release(output);
			}            
		}
	}
}

if (GPU == false) {
	Ext.CLIJ2_pull(output);
    Ext.CLIJ2_release(output);
}
selectImage(output);
setBatchMode("exit and display");
//saveAs("Tiff", out_folder + File.separator + file_name);
//
//Ext.CLIJ2_clear();
//run("Close All");
