// (17, 9484, 9962) image 
// NVIDIA RTX A4500: 42.1s in pyimagej macro split 2x2
// Intel(R) UHD Graphics 770: 1m59.4s in pyimagej macro split 2x2


run("CLIJ2 Macro Extensions", "cl_device=[13th Gen Intel(R) Core(TM) i9-13900] automatic_output_naming=false");
//run("CLIJ2 Macro Extensions", "cl_device=[NVIDIA RTX A4500] automatic_output_naming=false");

Ext.CLIJ2_clear();
GPU=false
radius_x = 5.0;
radius_y = 5.0;
sigma = 20.0;

input = getTitle();
Stack.getDimensions(width, height, channels, slices, frames);

if (GPU == true) {
	output = "output";
	newImage(output, "16-bit black", width, height, 1);
    } 
    else {
	Ext.CLIJ2_push(input);
    //run("Close All");
	Ext.CLIJ2_create2D(output, width, height, 16);
    }

numTilesZ = 1;
numTilesX = 5;
numTilesY = 5;

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
// Ext.CLIJ2_clear();
//selectImage(output);
