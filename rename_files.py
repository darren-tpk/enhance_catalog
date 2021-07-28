import os
import subprocess

directory = '/home/ptan/enhance_catalog/data/mammoth/'

old_filenames = os.listdir(directory)

for old_filename in old_filenames:

    old_filepath = directory + old_filename

    pieces = old_filename.split('.')
    new_filename = '.'.join([pieces[1],pieces[0],pieces[2]])
    new_filepath = directory + new_filename

    rename_command = 'mv ' + old_filepath + ' ' + new_filepath
    subprocess.call(rename_command, shell=True)