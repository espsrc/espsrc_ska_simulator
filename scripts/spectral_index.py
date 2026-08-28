import argparse

import astropy.units as u
from ska_img import SKAImage

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate spectral index between  images"
    )
    parser.add_argument("images", metavar="IMAGE", type=str, nargs="+")
    parser.add_argument(
        "--map",
        action="store_true",
        help="Output spectral index map instead of average value",
        default=False,
    )

    args = parser.parse_args()
    if len(args.images) < 2:
        raise Exception(
            "Please provide exactly two images to calculate the spectral index"
        )

    images = []
    for img_path in args.images:
        img = SKAImage(img_path)
        img.header["BUNIT"] = "Jy/beam"
        img.data_to_pixel()
        images.append(img)
    images = sorted(images, key=lambda x: x.restfreq())

    freqs = [img.restfreq().to(u.Hz) for img in images]

    # center_pixels = [img.data2d[500,500] for img in images]
    # #Fit using log-log with numpy polyfit
    # log_freqs = np.log10([f.value for f in freqs])
    # log_values = np.log10(center_pixels)
    # slope, intercept = np.polyfit(log_freqs, log_values, 1)
    # print (slope)
    # sys.exit()
    if args.map:
        si_map = SKAImage.spectral_index(images, map_mode=True)
        outname = "spectral_index_map.fits"
        si_map.tofits(outname, overwrite=True)
        print(f"Spectral index map written to {outname}")
    else:
        si = SKAImage.spectral_index(images)
        for i, img in enumerate(images):
            # Extract filename (without path)
            name = img.path.split("/")[-1]
            print(f"{name}: {img.janskys()} {img.restfreq().to(u.MHz)}")
        print(f"Spectral Index: {si}")
