from ska_img import SKAImage
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import argparse

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description='Display FITS images with WCS support.')
   parser.add_argument('fits', metavar='FITS', type=str, nargs='+', help='FITS file(s) to display')
   parser.add_argument('--exp', type=float, default=0.3, help='Exponent for scaling (default: 0.3)')
   args = parser.parse_args()

   fits = args.fits
   images = [SKAImage(f) for f in fits]

   nRows = (len(images) - 1) // 2 + 1
   nCols = 2 if len(images) > 1 else 1

   fig = plt.figure(figsize=(16, 8 * nRows))
   plt.subplots_adjust(hspace=0.5, wspace=0.5)
   for idx, img in enumerate(images):
      name = img.path.split('/')[-2]
      name = (' ').join(name.split('_')[2:])
      img.name = name
      zoom = None
      zoom =["10:00:27.42","+02:20:59.28","69arcsec"]
      zoomImg = img.zoom(zoom)
      frame = fig.add_subplot(nRows, nCols, idx + 1, projection=WCS(zoomImg.header))
      zoomImg.draw(img=frame, exp=args.exp, cmapstr='inferno', scale='power', colorbar=True, plot=False, title=img.name)
      print (f"{zoomImg.janskys():.2e}")
      # Rotate x-axis labels for better readability
      frame.coords['ra'].set_ticklabel(rotation=0, va='bottom', ha='left', pad=0)

   plt.tight_layout()
   plt.show()
