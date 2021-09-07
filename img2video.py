
import cv2
import os
from os.path import isfile, join

pathIn = '/Users/darrentpk/Desktop/frames/'
pathOut = '/Users/darrentpk/Desktop/mammoth3.avi'
fps = 12
frame_array = []
files = [f for f in os.listdir(pathIn) if isfile(join(pathIn, f))]
files.sort()
files.pop(0)
frame_array = []
for i in range(len(files)):
    filename = pathIn + files[i]
    # reading each files
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width, height)

    # inserting the frames into an image array
    frame_array.append(img)
out = cv2.VideoWriter(pathOut, cv2.VideoWriter_fourcc(*'DIVX'), fps, size)
for i in range(len(frame_array)):
    # writing to a image array
    out.write(frame_array[i])
out.release()