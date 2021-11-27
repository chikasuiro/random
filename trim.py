import glob
import cv2

lst = glob.glob('*.PNG')

for filename in lst:
    img = cv2.imread(filename)
    img_trim = img[top:bottom, left:right]
    outfilename = filename.replace('.PNG', '_trim.png')
    cv2.imwrite('trim/'+outfilename, img_trim)
