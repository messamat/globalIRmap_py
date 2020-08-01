from utility_functions import *

#Directory structure
rootdir = os.path.dirname(os.path.abspath(__file__)).split('\\src')[0]
datdir = os.path.join(rootdir, 'data')

gsimdir = os.path.join(datdir, 'GSIM')

url_metadata = "http://store.pangaea.de/Publications/DoH-etal_2018/GSIM_metadata.zip"
url_indices = "http://store.pangaea.de/Publications/GudmundssonL-etal_2018/GSIM_indices.zip"

dlfile(url_metadata, outpath=gsimdir, outfile = os.path.split(url_metadata)[1], ignore_downloadable=True)
dlfile(url_indices, outpath=gsimdir, outfile = os.path.split(url_indices)[1], ignore_downloadable=True)