%estimate psf for every available bead image and write to appropriate file
bead_file_struct = dir('/Volumes/shared/HCSS1/Shared299/JHOS/PSF BEAD IMAGES_correct/*.ome.tif');
for j = 1:length(bead_file_struct)
    pathToMovie = bead_file_struct(j).name;
    jonathanEstimateSigma(pathToMovie,1,1);
end
