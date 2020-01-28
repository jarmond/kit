%estimate psf for every available bead image and write to appropriate file
bead_file_struct = dir('/run/user/1001/gvfs/smb-share:server=ads.warwick.ac.uk,share=shared/HCSS1/Shared299/JHOS/PSF BEAD IMAGES_correct/*.ome.tif');
for j = 1:length(bead_file_struct)
    pathToMovie = fullfile(bead_file_struct(j).folder,bead_file_struct(j).name);
    jonathanEstimateSigma(pathToMovie,0,1);
end
