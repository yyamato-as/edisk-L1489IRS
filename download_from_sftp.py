import getpass
import fnmatch
import itertools
import os
import datetime
import time
import sys
import argparse
import tarfile
import subprocess

try:
    import paramiko
except ModuleNotFoundError:
    print("paramiko module not found. Please install as e.g., 'pip install paramiko'.")

default_ssh_key_filename = '/home/yamato/.ssh/id_rsa'

########## Authentic information ##########
hostname = "almaftp2.asiaa.sinica.edu.tw"
port = "2562"
username = "almauser"
# sftp_command = f"sftp -P {port} {username}@{hostname}"
key_filename = default_ssh_key_filename

########## What you want to download ##########

# parse the args
parser = argparse.ArgumentParser(description='Script to download data from eDisk sftp server')
parser.add_argument("source", help="Source name")
parser.add_argument("combine", choices=["SB", "SBLB"], help="Data set combination, SB or SBLB") 
parser.add_argument("data", choices=["image", "ms"], help="Data type, 'image' or 'ms'")
parser.add_argument("-d", "--desc", default="*", nargs="*", help="Data description, 'continuum' or line names. Multiple specification supported.")
parser.add_argument("-b", "--block", default="*", nargs="*", help="Execution block index for ms specification. Multiple specification supported.")
parser.add_argument("-p", "--imparam", default="*", nargs="*", help="Imaging parameter description. Multiple specification supported.")
parser.add_argument("-t", "--imtype", default="*", nargs="*", help="Type of image. Multiple specification supported.")
parser.add_argument("-f", "--filename", default=None, nargs="*", help="Direct specificatino of filenames. Multiple specification supported.")
parser.add_argument("--overwrite", action="store_true", help="Overwrite if the file exists already in the path.")
parser.add_argument("--dryrun", action="store_true", help="Dryrun to check infomation about the data that will be downloaded.")
parser.add_argument("--untar", action="store_true", help="Untar the compressed files if any.")
parser.add_argument("--path", default="./", help="Path to download the data.")
parser.add_argument("--timeout", default=60, help="Timeout for sftp connection.")

args = parser.parse_args()

########## setup the ssh client and got connection ##########
# get the input password for ssh private key
password = getpass.getpass(f"Enter passphrase for 'key {key_filename}':")

# set up client
client = paramiko.SSHClient()
client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
client.connect(hostname, port=port, username=username, password=password, key_filename=key_filename, timeout=args.timeout)
sftp = client.open_sftp()
print(f"Connected to the sftp server {hostname}.")

########## Download the files ##########
# a few useful functions
def get_remotepath(source, combine, data):
    combine = "ShortBaseline" if combine == "SB" else "LongBaseline"
    data = "eDisk_image_products" if data == "image" else "eDisk_calibrated_data"

    return f"/data/{source}/{combine}/{data}"

def get_filename_wc(source, combine, data, desc=["*"], block=["*"], imparam=["*"], imtype=["*"], filename=None):
    if filename is not None:
        return filename

    if data == "ms":
        return [f"{source}_{b}_{d}.ms.tgz" for b, d in itertools.product(block, desc)]

    elif data == "image":
        return [f"{source}_{combine}_{st}_{ip}.{it}.*" for st, ip, it in itertools.product(desc, imparam, imtype)]

remotepath = get_remotepath(args.source, args.combine, args.data)
sftp.chdir(remotepath)
filename_wc_list = get_filename_wc(args.source, args.combine, args.data, desc=args.desc, block=args.block, imparam=args.imparam, imtype=args.imtype, filename=args.filename)

to_download = []
for filename_wc in filename_wc_list:
    for filename in sftp.listdir():
        if fnmatch.fnmatch(filename, filename_wc):
            to_download.append(filename)

if args.dryrun:
    total_datavol = 0.0
    print(f"This script will download {len(to_download)} files below into {args.path}:")
    for filename in to_download:
        fstat = sftp.stat(filename)
        size = fstat.st_size*1e-9
        total_datavol += size
        mtime = datetime.datetime.fromtimestamp(fstat.st_mtime).strftime('%x %X')
        print("{:.3f} GB \t {:s} \t {:s}".format(size, mtime, filename))
    print("Total data volume: {:.2g} GB".format(total_datavol))
    out = input("Run actual download? (yes/no):")
    if out != "yes":
        print("Dryrun done.")
        sys.exit()

print(f"Starting to download {len(to_download)} files into {args.path}...")

start = time.time()

downloaded = []
for i, filename in enumerate(to_download):
    if not args.overwrite and os.path.exists(args.path + filename):
        print("[{:d}/{:d}] Skipped {:s} as it already exists.".format(i+1, len(to_download), filename))
        continue
    else:
        print("[{:d}/{:d}] Downloading {:s}...".format(i+1, len(to_download), filename))
        sftp.get(filename, args.path+filename)
        downloaded.append(filename)

end = time.time()
print("Download completed successfully in {}s.".format(end-start))

# close the sftp session
sftp.close()
client.close()

# check if any compressed files
to_untar = [filename for filename in downloaded if filename.endswith((".tgz", ".gz"))]
if len(to_untar) != 0 and args.untar:
    print(f"Will unzip {len(to_untar)} compressed files...")
    for i, filename in enumerate(to_untar):
        print("[{:d}/{:d}] Unzipping {:s}...".format(i+1, len(to_untar), filename))
        try:
            with tarfile.open(args.path+filename, mode="r") as tar:
                tar.extractall(path=args.path)
                os.remove(args.path+filename)
        except tarfile.ReadError:
            subprocess.run(["gunzip", args.path+filename])
    print("Unzip done.")


