import warnings
warnings.simplefilter('ignore')
import sys, os, time
import argparse, json
sys.path.append('%s/../src' % (os.path.dirname(os.path.realpath(__file__))))
from dendrosources import *
from commons import show_exc
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000
import astropy.units as u
mpl.rcParams['figure.figsize']=(10,10)
from scipy.ndimage.measurements import center_of_mass
from astrodendro.analysis import PPStatistic
from astrodendro import Dendrogram
import multiprocessing as mp

def dominated_by(spec, key):
    try:
        if spec is None:
            return "-"
        if key not in spec.keys():
            return "-"
        item = spec[key]
        freefree= float(item[4])
        dust = float(item[5])
        unknown = float(item[6])

        if freefree + dust + unknown == 0:      #ONLY in B6
            return "d"

        if freefree / max(dust, unknown, 1)  > 3:
            return "f"
        if dust / max(freefree, unknown, 1) > 3:
            return "d"
        return "-"
    except Exception as e:
        print (show_exc(e))
        return "E"

def zoom(img):
    (y,x) = img.data.shape
    rows = []
    cols = []

    for i in range(y):
        if (not (np.isnan(img.data[i,:]).all())):
            rows.append(i)
    for j in range(x):
        if (not (np.isnan(img.data[:,j]).all())):
            cols.append(j)

    img.header['CRPIX1'] -= min(cols)
    img.header['CRPIX2'] -= min(rows)
    img.data = img.data[min(rows):max(rows),min(cols):max(cols)]
    return (img)

def search_pbfile(project, region):
    folders={"bgps":"../data/ALMA-IMF-B6","mgps": "../data/ALMA-IMF-B3"}

    try:
        aux = region.name.split(' ')[0]
        files = glob.glob(os.path.join(folders[project], "{}**.fits".format(aux)))
        for f in files:
            if "pb.tt0.fits" in f:
                pb = f
                return (Image(pb))
            if "image.tt0.fits" in f:
                image = f
            if "tt0.pbcor.fits" in f:
                pbcor = f

        image = Image(image)
        pbcor = Image(pbcor)
        pb_model = Image(data=1/(pbcor.data/image.data), header=pbcor.header)
        return (pb_model)
    except Exception as e:
        print (show_exc(e))
        return None

def clump_flux(node, orig, flux_parent=0*u.Jy, corrected=True):
    try: 
        import copy
        img = copy.deepcopy(orig)
        if img is not None and node is not None:
            if not corrected:
                return img.janskys()

            if node.level == 0:
                img.data[~(node.get_mask())] = np.nan

                for i in node.descendants:
                    img.data[(i.get_mask())] = np.nan

                mean = np.nanmean(img.data)
                for i in node.descendants:
                    img.data[(i.get_mask())] = mean

                flux = img.janskys()
                return flux
            else:
                img.data[~(node.get_mask())] = np.nan
                return img.janskys() - (flux_parent / node.ancestor.get_npix() * node.get_npix())
        else:
            return (0*u.Jy)
    except Exception as e:
        mylog(show_exc(e))
        return (0*u.Jy)

def process_img(params):
    try:
        f = params[0]
        args = params[1]
        img = Image(f)
        region_name = img.name.split(' ')[0]
        try:
            distance = u.Quantity(args.distance)
        except:
            distance = 1 * u.kpc

        if args.posfix != "":
            posfix = "_{}".format(args.posfix)
        else:
            posfix = ""

        if args.pblim < 1.:
            trunk_color = 'white'

            mad = img.mad() # Calculate mad before crop the data, from full combination
            if "BOLOCAM" in f:
                if 'combination' in f:
                    pb = search_pbfile("bgps", img)
                    img.data[(pb.data < args.pblim)] = np.nan
                project = "B6"
            else:
                if 'combination' in f:
                    pb = search_pbfile("mgps", img)
                    img.data[(pb.data < args.pblim)] = np.nan
                project = "B3"
        else:
            mad = 3.55384e-4
            project = 'ALPHA'
            noise = np.random.normal(mad,0.1, img.data.shape)
            #img.data[(np.isnan(img.data))] = noise[(np.isnan(img.data))]
            img.data[(img.data == 0.5)] = noise[(img.data == 0.5)]
            img.data[(img.data == 0.)] = 8.*mad
            img.data[(img.data == 1.)] = 8.*mad
            trunk_color = 'blue'

        mylog(region_name, project)
        n_pixels = int((img.omega_beam()/img.omega_pix()).value) * 1.5
        d = Dendrogram.compute(img.data, min_value= args.min * mad, min_delta=mad * args.delta , min_npix=n_pixels)
        d.save_to("DDRG_{}-{}-dend_delta_{}_min_{}{}.fits".format(region_name, project, args.delta, args.min,posfix))
        leaves = d.leaves
        items = sorted(d.all_structures, key =lambda x:x.level)
        max_level = items[-1].level + 1
        cmap = matplotlib.cm.get_cmap('PuOr', max_level)
        mylog("{} ({}) --> {} sources".format(region_name, project, len(leaves)))
        coords = []
        mylist = []
        metadata = {}
        metadata['data_unit'] = u.Jy / u.beam
        metadata['spatial_scale'] = img.pixelarea**(0.5)
        metadata['beam_major'] = img.get_beam().major
        metadata['beam_minor'] = img.get_beam().minor
        metadata['wcs'] = WCS(img.header)
        if args.png:
#             img = zoom(img)
            fig = plt.figure()
            axis = img.draw(fig, exp=args.exp, title = '{} {} min={:.2f} delta={:.2f} mad={:.4E}'.format(img.name, project, args.min, args.delta, mad))
            p = d.plotter()
            ax = fig.get_axes()[0]

            if args.catalogues:
                for path in args.catalogues:
                    if region_name.split('.')[0] in path:
                        f = open(path, 'r')
                        lines = [line for line in f.readlines() if (not line.startswith('#') and not line.startswith('!'))]
                        for line in lines:
                            try:
                                fields = ' '.join(line.split(' ')).split()
                                if len(fields) > 6:
                                    coord = SkyCoord(fields[6], fields[7], frame=args.frame, unit=(u.deg, u.deg))
                                    coords.append(coord)
                            except Exception as e:
                                mylog(line)
                                mylog(show_exc(e))

                        counter = 0
                        out_of_region = 0
                        for coord in coords:
                            (x,y) = img.coords2pix(coord)
                            if x < img.data.shape[0] and y < img.data.shape[1] and x >= 0 and y >= 0:
                                counter += 1
                                mylog("{:.2f}%\r".format(counter/(len(coords))*100), end="")
                                axis.show_markers(coord.ra, coord.dec, edgecolor="black", facecolor="black", marker="s", alpha=0.5)
                                if img.data[int(x),int(y)] == np.nan:
                                    out_of_region += 1

                mylog("\n\n{} ({}) sources in {} ({:.2f}%)".format(counter, counter-out_of_region, img.name, counter/(len(coords)) * 100.))

            counter = 0
            for node in items:
                try:
                    stat = PPStatistic(node, metadata=metadata)
                    if node.level == 0:
                        if node.is_leaf:
                            counter+=1
                            p.plot_contour(ax, structure=node, color="orange")
                            coord_cm = img.center_of_mass(node.get_mask())
                            if args.numbers:
                                axis.add_label(coord_cm.ra.to(u.deg).value, coord_cm.dec.to(u.deg).value, f"{node.idx}",color="blue")
                            else:
                                axis.show_markers(coord_cm.ra, coord_cm.dec, edgecolor="black", facecolor="black", marker='x')

                        else:
                            if not args.onlyleaf:
                                p.plot_contour(ax, structure=node, color=trunk_color)
                                coord_cm = img.center_of_mass(node.get_mask())
                                if args.numbers:
                                    axis.add_label(coord_cm.ra.to(u.deg).value, coord_cm.dec.to(u.deg).value, f"{node.idx}",color="blue")
                                else:
                                    axis.show_markers(coord_cm.ra, coord_cm.dec, edgecolor="black", facecolor="black", marker='o')
                        mylist.append(node)
                    else:
                        if node.is_leaf and args.pblim < 1.:
                            counter+=1
                            p.plot_contour(ax, structure=node, color="red")
                            coord_cm = img.center_of_mass(node.get_mask())
                            if args.numbers:
                                axis.add_label(coord_cm.ra.to(u.deg).value, coord_cm.dec.to(u.deg).value, f"{node.idx}",color="black")
                            else:
                                axis.show_markers(coord_cm.ra, coord_cm.dec, edgecolor="black", facecolor="black", marker='x')
                            mylist.append(node)
                except Exception as e:
                    mylog(show_exc(e))

            #plt.plot(0,0,'x', label='baricenter', color='black')
            if args.catalogues:
                plt.plot(0,0,'s', label='catalogue', color='black')
            #plt.legend(loc='upper center', borderaxespad=0)
            if posfix != "":
                plt.savefig("{}-{}-dend{}.png".format(region_name, project, posfix))
            else:
                plt.savefig("{}-{}-dend_mad_{}_pblim_{}_delta_{}_min_{}{}.png".format(region_name, project, args.mad, args.pblim, args.delta, args.min,  posfix))
            mylist = items
            #img.tofits(f'zoom_{region_name}_{project}.fits')
        else:
            mylist = items


        cores = 0
        clumps = 0

        clump_list = []
        clump_list_flux = []

        for item in mylist:
            if item.level == 0:
                clump_list.append(item)
                clump_list_flux.append(clump_flux(item, img))

        table1 = []
        table2 = []
        sizes_pc = []
        spec_array = []
        for idx,node in enumerate(mylist):
            coord_cm = img.center_of_mass(node.get_mask())
            stat = PPStatistic(node, metadata=metadata)
            area = (node.get_npix() * img.pixelarea.to(u.arcsec**2))
            radius = (area.value / np.pi)**(0.5)
            radius_pc = arcsec2pc(radius * u.arcsec, distance).to(u.pc).value
            ellip = 1 - (stat.minor_sigma/stat.major_sigma)


            if node.level == 0:
                if node.is_leaf:
                    cores += 1
                    table1.append("{}-{}\#TRXX-LV{:03d}  &   {:s} &   {:s} & {:.2f} & {:.2f} & {:.2f}, {:.3f} & {:.2f} & {:.2f} & __DUST__ \\\\ \n".format(region_name,project,node.idx, coord_cm.ra.to_string(unit=u.hour, sep=("hms"), precision=2), coord_cm.dec.to_string(unit=u.degree, sep=("dms"), precision=2), area.value, ellip.value, radius, radius_pc, clump_flux(node, img).to(u.mJy).value,stat.flux.to(u.mJy).value))
                else:
                    table1.append("{}-{}\#TR{:02d}  &   {:s} &   {:s} &   {:.2f} & - & {:.2f}, {:.3f} & {:.2f} & {:.2f} & __DUST__ \\\\ \n".format(region_name,project,node.idx, coord_cm.ra.to_string(unit=u.hour, sep=("hms"), precision=2), coord_cm.dec.to_string(unit=u.degree, sep=("dms"), precision=2), area.value, radius, radius_pc, clump_flux(node, img).to(u.mJy).value,stat.flux.to(u.mJy).value))
            else:
                if node.is_leaf:
                    cores += 1
                    ((x,y),val) = node.get_peak()
                    (a,b) = img.coords2pix(coord_cm)
                    dist = ( ((b-x)**2 + (a-y)**2) ** (0.5))

                    table1.append("{}-{}\#TR{:02d}-LV{:03d}  &   {:s} &   {:s} &   {:.2f} &   {:.2f} &   {:.2f}, {:.3f}   & {:.2f} & {:.2f} & __DUST__ \\\\ \n".format(region_name,project,node.ancestor.idx, node.idx, coord_cm.ra.to_string(unit=u.hour, sep=("hms"), precision=2), coord_cm.dec.to_string(unit=u.degree, sep=("dms"), precision=2), area.value, ellip.value, radius, radius_pc, clump_flux(node,img,clump_list_flux[clump_list.index(node.ancestor)]).to(u.mJy).value,stat.flux.to(u.mJy).value))
                    table2.append("{}-{}\#TR{:02d}-LV{:03d}  &   {:s} &   {:s} &   {:.1f} & {:.3f} \\\\ \n".format(region_name,project,node.ancestor.idx, node.idx, coord_cm.ra.to_string(unit=u.hour, sep=("hms"), precision=2), coord_cm.dec.to_string(unit=u.degree, sep=("dms"), precision=2), dist, (img.pixelarea.to(u.arcsec**2)**(0.5)).value))

            if node.is_leaf:
                sizes_pc.append(radius_pc * 2)

        fig, axs = plt.subplots(1,1)
        N, bins, patches = axs.hist(sizes_pc)
        plt.title('{} {} min={:.2f} delta={:.2f}'.format(img.name, project, args.min, args.delta))
        axs.axvline(args.coresize, color='r', linestyle='dashed')
        plt.xlabel('[pc]')
        plt.ylabel('Number of cores')
        plt.savefig("{}-{}-dend_delta_{}_min_{}{}_hist.png".format(region_name, project, args.delta, args.min,posfix))



        return [table1, table2]
    except Exception as e:
        mylog(show_exc(e))
        return [[],[]]

def write_tex(fname, tables, headers, rows_per_pages=70, spec = {}):
    f = open(fname, "w")
    f.write("\\documentclass[a4paper,12pt,oneside]{article}\n")
    f.write("\\usepackage{geometry}\n")
    f.write("\\geometry{ a4paper, total={170mm,257mm}, left=20mm, top=10mm}\n")
    f.write("\\usepackage{tabularx,booktabs}\n")
    f.write("\\begin{document}\n")
    f.write("The ALMA-IMF data have been cropped to PB=0.4\n")

    for lines in tables:
        lines = sorted(lines)
        pages = int((len(lines) / rows_per_pages)+1 )
        for idx in range(pages):
            if (len(lines) > 0):
                region_name = lines[0].split('#')[0]
                band = region_name.split('-')[-1].replace('\\','')
                for header in headers:
                    f.write(header.replace('__IDREGION__', region_name))
                f.write("\\hline \n")
                f.write("\\hline \n")
                for idy in range(rows_per_pages):
                    try:
                        if len(lines) > 0:
                            line = lines.pop(0)
                            if "LV" not in line or "TRXX" in line:
                                f.write("\\hline \n")

                            if '__DUST__' in line:
                                fields = line.split('&')
                                if (len(fields) > 1):
                                    fields=[field.strip() for field in fields]
                                    key = '{}{}{}'.format(band,fields[1],fields[2])
                                    dom = dominated_by(spec, key)
                                    line=line.replace('__DUST__', dom)
                                else:
                                    line=line.replace('__DUST__','')
                            f.write(line.split('#')[1])
                        else:
                            break
                    except Exception as e:
                        print (line)
                        print (show_exc(e))
                f.write("\\hline \n")
                f.write("\\end{tabular} \\end{table} \n")
                f.write("\\newpage\n")
    f.write("\\end{document}\n")
    f.close()

def write_csv(fname, tables):
    region_name = ""
    f = None
    for lines in tables:
        lines = sorted(lines)
        for idx,line in enumerate(lines):
            current_region= line.split('#')[0].replace('\\','')

            if (current_region != region_name):
                if f:
                    f.close()
                f = open(f'{fname}_{current_region}.csv','a')
                region_name = current_region


            fields = line.split('#')[1].replace('\\','').split('&')
            fields=[field.strip() for field in fields]
            f.write('{}\n'.format(';'.join([current_region.split('-')[-1]]+fields+[region_name])))
    if f:
        f.close()

class Params():
    def __init__(self, args):
        self.png = args.png
        self.delta = args.delta
        self.pblim = args.pblim
        self.exp = args.exp
        self.min = args.min
        self.catalogues = args.catalogues
        self.frame = args.frame
        self.posfix = args.posfix
        self.onlyleaf = args.onlyleaf
        self.cpus = args.cpus
        self.numbers = args.numbers
        self.histsizes = args.histsizes
        self.coresize = args.coresize
        self.mad = args.mad

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Dendrogram')
    parser.add_argument('file', type=str, nargs="+")
    parser.add_argument('--png', action='store_true', default=False)
    parser.add_argument('--delta', type=float, default=5.0, help='delta * MAD => delta between trunks')
    parser.add_argument('--pblim', type=float, default=0.5, help='')
    parser.add_argument('--exp', type=float, default=0.5)
    parser.add_argument('--min', type=float, default=5.0, help='min * MAD => min value')
    parser.add_argument('--mad', type=float, default=1.0, help='MADn => MADn')
    parser.add_argument('--catalogues', type=str, nargs="*", help='Catalogues')
    parser.add_argument('--frame', type=str, default="icrs", help="Frame")
    parser.add_argument('--posfix', type=str, default="", help="Prefix")
    parser.add_argument('--distance', type=str, default="1kpc", help="Distance with units in astropy.units format (default='1kpc')")
    parser.add_argument('--onlyleaf', action='store_true', default=False)
    parser.add_argument('--cpus', type=int, default=1)
    parser.add_argument('--numbers', action='store_true', default=False)
    parser.add_argument('--histsizes', action='store_false', default=True)
    parser.add_argument('--coresize', type=float, default=0.1, help="Max Size (in pc) for cores. Default = 0.1")
    parser.add_argument('--spec_file', type=str, default=None, help="Text File with spectral index info")

    args = parser.parse_args()
    
    try:
        params = []
        for f in args.file:
            if f.endswith(".json"):
                lines = open(f,'r').readlines()
                lines = [item for item in lines if not item.startswith('#')]
                for line in lines:
                    jsonstr = {}
                    jsonstr = json.loads(line)
                    par = Params(args)
                    par.pblim = float(jsonstr['pblim'])
                    par.min = float(jsonstr['min'])
                    par.delta= float(jsonstr['delta'])
                    par.distance=jsonstr['distance']
                    try:
                        par.mad = jsonstr['mad']
                    except:
                        pass
                    if 'spec_index' in jsonstr.keys():
                        par.spec_index = jsonstr['spec_index']
                    params.append([jsonstr["file"],par])
            else:
                params.append([f,Params(args)])


        cpus = min(args.cpus, int(0.8 * mp.cpu_count()))
        with mp.Pool(cpus) as p:
            if args.posfix != "":
                posfix = "_{}".format(args.posfix)
            else:
                posfix = ""

            results = p.map(process_img,params)

            spec = {}
            if args.spec_file is not None:
                spec_file = open(args.spec_file, 'r')
                lines = spec_file.readlines()
                for line in lines:
                    fields = line.replace('\n','').split(' ')
                    spec['{}{}{}'.format(fields[0].split('-')[-1],fields[1],fields[2])] = fields
            lines = []
            lines2 = []

            headers = []
            headers.append("\\begin{table}[!ht] \\scriptsize \\centering \\begin{tabular}{l l l r r r r r r} \\hline \n" )
            for item in results:
                lines.append(item[0])
                lines2.append(item[1])

            ROW_PER_PAGES = 70
            headers.append("ID              & R.A. & Dec   & Area                 & Ellip & Radius       & ST\_flux & T\_flux & F/D \\\\ \n")
            headers.append(" __IDREGION__   &      &       & [arcsec$^2$]  &       & [arcsec, pc] & [mJy]    & [mJy] & \\\\ \n")
            write_tex("table{}.tex".format(posfix), lines, headers, ROW_PER_PAGES, spec=spec)

            #lines = sorted(lines2)
            headers = []
            headers.append(" \\begin{table}[!ht] \\scriptsize \\centering \\begin{tabular}{l l l l l } \\hline \n" )
            headers.append("ID & R.A. & Dec & Dist & Pix Size \\\\ \n")
            headers.append("   &      &     & [pix] & [arcsec]  \\\\ \n")

            write_tex("table_priv{}.tex".format(posfix), lines2, headers, ROW_PER_PAGES)

            for item in results:
                lines.append(item[0])
                lines2.append(item[1])

            write_csv("DDR_table{}".format(posfix), lines)


    except Exception as e:
        print (show_exc(e))

