import cv2

fourcc = cv2.VideoWriter_fourcc('m','p','4','v')
video = cv2.VideoWriter('transient.mp4', fourcc, 10.0, (1280, 720))

for i in range(1, 90):
    img_src = cv2.imread('transient-{}.ppm'.format(i))
    #img = cv2.resize(img, (1280,720))
    #img_src = img[119+50:839+50,:].copy()
    #img = cv2.resize(img_src, (1280,720))
    video.write(img_src)

for i in range(8):
    video.write(cv2.imread("main.ppm"))

video.release()
