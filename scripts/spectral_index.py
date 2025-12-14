from ska_img import SKAImage
import astropy.units as u
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate spectral index between  images")
    parser.add_argument("images", metavar='IMAGE', type=str, nargs='+')
    parser.add_argument("--map", action='store_true', help="Output spectral index map instead of average value", default=False)


    args = parser.parse_args()
    if len(args.images) < 2:
        raise Exception("Please provide exactly two images to calculate the spectral index")
    
    images = []
    for img_path in args.images:
        img = SKAImage(img_path)
        images.append(img)
    images = sorted(images, key=lambda x: x.restfreq())
    if args.map:
        si_map = SKAImage.spectral_index(images, map_mode=True)
        outname = "spectral_index_map.fits"
        si_map.tofits(outname, overwrite=True)
        print (f"Spectral index map written to {outname}")
    else:
        si = SKAImage.spectral_index(images)
        for i, img in enumerate(images):
            # Extract filename (without path)
            name = img.path.split("/")[-1]
            print (f"{name}: {img.janskys()} {img.restfreq().to(u.MHz)}")
        print (f"Spectral Index: {si}")
